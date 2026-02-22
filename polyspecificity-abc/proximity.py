#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import json
import os
import gzip
import argparse
from scipy.spatial.distance import cdist
import logging
logging.basicConfig(level=logging.INFO)

# get instances (asym_ids) and find nearby other instances (within threshold distance of any atom)

parser = argparse.ArgumentParser()
parser.add_argument('pdb_id', nargs='*', help='PDB ID, e.g. 2kgp 4KWX')
parser.add_argument('--in-tsv', help='input TSV file containing a pdb_id column')
parser.add_argument('--excel', help='output in Excel format')
parser.add_argument('--csv', help='output in CSV format')
parser.add_argument('--tsv', help='output in TSV format')
parser.add_argument('--stdout', action='store_true', default=False, help='output in TSV format to stdout')
parser.add_argument('--threshold', type=float, default=5.0, help='Distance threshold')
parser.add_argument('--pdb-repo', default='test-pdb', help='PDB repository path')

args = parser.parse_args()

PDB_REPO = args.pdb_repo

if (not args.pdb_id) and (not args.in_tsv):
    # quit if we have no structures
    parser.print_help(sys.stderr)
    sys.exit(1)

all_pdb_ids = list(args.pdb_id)
if args.in_tsv:
    structures = pd.read_csv(args.in_tsv, dialect='excel-tab')
    all_pdb_ids.extend(list(sorted(structures['pdb_id'].unique())))

results = []
for pdb_id in all_pdb_ids:
    pdb_id = pdb_id.lower()
    logging.info(f'**** {pdb_id}')

    assy_fns = [(0, f'{PDB_REPO}/mmjson/{pdb_id[1:3]}/{pdb_id}.json.gz')]
    assy_id = 1
    fn = f'{PDB_REPO}/biounit-mmjson/{pdb_id[1:3]}/{pdb_id}-assembly{assy_id:d}.json.gz'
    while os.path.exists(fn):
        assy_fns.append((assy_id, fn))
        assy_id += 1
        fn = f'{PDB_REPO}/biounit-mmjson/{pdb_id[1:3]}/{pdb_id}-assembly{assy_id:d}.json.gz'

    for assy_id, json_fn in assy_fns:
        try:
            with gzip.open(json_fn, 'rt', encoding='utf-8') as f:
                j = json.load(f)
                j = j[list(j.keys())[0]]
        except json.decoder.JSONDecodeError as e:
            logging.error(f'decode error in {pdb_id}')
            continue
        except FileNotFoundError as e:
            logging.error(f'structure mmjson file not found: {json_fn}')
            continue

        struct_asym = pd.DataFrame(j['struct_asym'])
        atom_site = pd.DataFrame(j['atom_site'])

        for asym_id in struct_asym['id']:
            query_coords = atom_site.loc[(atom_site['label_asym_id'] == asym_id) & (atom_site['type_symbol'] != 'H'), ['Cartn_x', 'Cartn_y', 'Cartn_z']].to_numpy()
            search_atoms = atom_site.loc[(atom_site['label_asym_id'] != asym_id) & (atom_site['type_symbol'] != 'H')]
            search_coords = search_atoms[['Cartn_x', 'Cartn_y', 'Cartn_z']].to_numpy()
            distances = cdist(search_coords, query_coords, 'euclidean')
            distances = np.min(distances, axis=1)
            distances = pd.Series(distances, index=search_atoms.index)
            idx = (distances < args.threshold)
            found_atoms = search_atoms.loc[idx].copy()
            found_atoms['distance'] = distances.loc[idx]
            found_atoms = found_atoms.sort_values('distance')
            found_asym_ids = found_atoms[['label_asym_id', 'distance']].drop_duplicates('label_asym_id', ignore_index=True)
            found_asym_ids['threshold'] = args.threshold
            found_asym_ids.insert(0, 'pdb_id', pdb_id)
            found_asym_ids.insert(1, 'asym_id', asym_id)
            found_asym_ids['assembly_id'] = assy_id
            results.append(found_asym_ids)
results = pd.concat(results)
results = results.groupby(['pdb_id', 'asym_id', 'label_asym_id']).agg({
    'distance': 'min',
    'threshold': 'first',
    'assembly_id': lambda a: ','.join(map(str, sorted(a))),
}).reset_index()

if args.excel:
    results.to_excel(args.excel, sheet_name='sheet', index=False)

if args.csv:
    results.to_csv(args.csv, index=False)

if args.tsv:
    results.to_csv(args.tsv, index=False, sep='\t')

if args.stdout:
    print(results.to_csv(args.csv, index=False, sep='\t'))

