#!/usr/bin/env python3

import logging
logging.basicConfig(level=logging.INFO)

import os
import sys

import pandas as pd
import argparse

import gzip
import json
import csv
import numpy as np
from scipy.spatial.distance import cdist

import itertools as it

from collections import OrderedDict, Counter

parser = argparse.ArgumentParser()
parser.add_argument('ligands_fn', help='Table with pdb_id, mon_id, asym_id, can be TSV or Excel, or a ligand_key')
parser.add_argument('proximity_tsv')
parser.add_argument('--mon-id', help='Only consider ligands from ligands_fn with specified mon_id')
parser.add_argument('--pdb-repo', default='test-pdb', help='PDB repository path')

args = parser.parse_args()

PDB_REPO = args.pdb_repo

# simple cache object for mmcif structure caching
class Cache:
    def __init__(self, max_size):
        self.max_size = max_size
        self.cache = OrderedDict()

    def get(self, key):
        return self.cache.get(key)

    def put(self, key, value):
        # Remove the key if it already exists to update its position
        if key in self.cache:
            self.cache.pop(key)
        # If we're at max capacity, pop the oldest item
        elif len(self.cache) >= self.max_size:
            self.cache.popitem(last=False)
        # Insert the new item; it will appear at the end (most recent)
        self.cache[key] = value

# cache last 10 files
mmjson_cache = Cache(10)

# string justification helper function for mmcif file writing
def ljust(s):
    s = s.astype(str)
    # fields containing newlines will be in a new line anyway
    return s.str.ljust(s[~s.str.contains('\n')].str.len().max())
    
if args.ligands_fn.endswith('.tsv'):
    ligands = pd.read_csv(args.ligands_fn, dialect='excel-tab')
elif args.ligands_fn.endswith('.xlsx'):
    ligands = pd.read_excel(args.ligands_fn)
else:
    # pretend that we have a ligand_key
    ligands = pd.DataFrame([args.ligands_fn.split('_')], columns=['pdb_id', 'mon_id', 'asym_id'])
    logging.warning('interpreting as ligand_key: %s', args.ligands_fn)

if args.mon_id:
    ligands = ligands[ligands['mon_id'] == args.mon_id]

proximity = pd.read_csv(args.proximity_tsv, dialect='excel-tab')

# start from ligand asym_id, and select: close polymer or ligand instances, and other ligand instances close to these polymers
# we choose an assembly that contains the most components around the ligand
# then we choose neighboring components in that assembly
# generate pymol code to align on ligand atoms that are present
# the align command should work with ligands: https://www.compchems.com/align-structures-in-pymol/#align-based-on-the-protein-backbone
pymol_selectors = []
for ligand in ligands.itertuples():
    p = proximity[(proximity['pdb_id'] == ligand.pdb_id) & (proximity['asym_id'] == ligand.asym_id)].copy()
    p['assembly_id'] = p['assembly_id'].str.split(',')
    p = p.explode('assembly_id')
    p['assembly_id'] = p['assembly_id'].astype('int')
    #if ligand.pdb_id == '3edl':
    #    # we blacklist this assembly because it contains a symmetry operation very close to identity
    #    p = p[p['assembly_id'] != 1]
    #if ligand.pdb_id == '9eok':
    #    # skip this structure because there are very serious clashes and close contacts
    #    continue

    cnt = Counter(p['assembly_id'])
    try:
        assembly_id = cnt.most_common(1)[0][0]
        logging.info(f'**** {ligand.pdb_id} {ligand.mon_id} {ligand.asym_id}: using assembly {assembly_id}')
    except IndexError as e:
        # seems nothing is close by
        logging.warning(f'**** {ligand.pdb_id} {ligand.mon_id} {ligand.asym_id}: no nearby chains found')
        assembly_id = None
        continue
    
    # load cif of biounit assembly
    if assembly_id == 0:
        fn = f'{PDB_REPO}/mmjson/{ligand.pdb_id[1:3]}/{ligand.pdb_id}.json.gz'
    else:
        fn = f'{PDB_REPO}/biounit-mmjson/{ligand.pdb_id[1:3]}/{ligand.pdb_id}-assembly{assembly_id}.json.gz'
    if not os.path.exists(fn):
        logging.error(f'file {fn} not found!')
        continue
    j = mmjson_cache.get(fn)
    if j is None:
        with gzip.open(fn, 'rt', encoding='utf-8') as f:
            j = json.load(f)
            j = j[next(iter(j.keys()))]
        mmjson_cache.put(fn, j)

    p = p[p['assembly_id'] == assembly_id]

    included_asym_ids = set([ligand.asym_id] + p['label_asym_id'].unique().tolist())
    logging.info(f'**** {ligand.pdb_id} {ligand.mon_id} {ligand.asym_id}: initial close chains: {included_asym_ids!r}')
    # TODO: asym_ids *might not* have the form "letter-dash-symmetry operation",
    # as there is a pdbx_chain_remapping structure in assembly files that maps
    # to the original asym_id and the symmetry operation.
    # For the proximity tables, maybe the assembly asym_id should be mapped to
    # the original asym_id + dash + symmetry op for standardization.
    # Then, we would need to seek out the corresponding asym_id here.
    # Example: 6b0l with helical symmetries

    # I've commented this out now because entry 6b0l has symmetry operation 27 as identity
    #if 'pdbx_chain_mapping' in j:
    #    mapping = pd.DataFrame(j['pdbx_chain_mapping'])
    #    mapping['std_asym_id'] = (mapping['orig_label_asym_id'] + '-' + mapping['applied_operations'].astype('str')).str.replace('-1', '')
    #    # map IDs
    #    included_asym_ids = set(mapping.loc[mapping['std_asym_id'].isin(included_asym_ids), 'label_asym_id'])
    ## if mapping table is not present, just go with the asym_ids we have

    # now look at non-polymer chains that are close to polymers chains in the initial seletion
    entities = pd.DataFrame(j['entity'])
    struct_asym = pd.DataFrame(j['struct_asym'])
    struct_asym = struct_asym.merge(entities[['id', 'type']].rename(columns={'type': 'entity_type'}).set_index('id'),
                                               how='left', left_on='entity_id', right_index=True)
    included_poly_entities = included_asym_ids & set(struct_asym.loc[struct_asym['entity_type'] == 'polymer', 'id'])
    nonpoly_entities = set(struct_asym.loc[struct_asym['entity_type'] != 'polymer', 'id'])
    included_asym_ids |= set(proximity.loc[
        (
            (proximity['pdb_id'] == ligand.pdb_id) & 
            (proximity['asym_id'].isin(included_poly_entities)) & 
            (proximity['label_asym_id'].isin(nonpoly_entities))
        ),
        'label_asym_id'])
    logging.info(f'**** {ligand.pdb_id} {ligand.mon_id} {ligand.asym_id}: final close chains: {included_asym_ids!r}')

    # now we need to write every chain from included_asym_ids
    atom_site = pd.DataFrame(j['atom_site'])

    # here we create a mapping for the single-letter pdb chain IDs
    # TODO: is there a more pandas-friendly way to do this?

    # we will write pdb files that need single-letter chain IDs
    chainid_map = dict()
    new_chain_ids = set()
    if 'pdbx_chain_remapping' in j:
        chain_remapping = pd.DataFrame(j['pdbx_chain_remapping'])
        for r in chain_remapping.itertuples():
            new_chain_id = str(r.orig_auth_asym_id)[0]
            if new_chain_id in new_chain_ids:
                # try lower case
                new_chain_id = new_chain_id.lower()
            #new_chain_id = r['orig_label_asym_id']
            #new_chain_id = r['label_asym_id']
            #if new_chain_id.endswith('-2'):
            #    new_chain_id = r['orig_label_asym_id'].lower()
            chainid_map[r.label_asym_id] = new_chain_id
            new_chain_ids.add(new_chain_id)
    else:
        for r in atom_site[['label_asym_id', 'auth_asym_id']].drop_duplicates().itertuples():
            # some structures contains numbers as asym_ids, and those are retained as ints in json
            new_chain_id = str(r.auth_asym_id)[0]
            if new_chain_id in new_chain_ids:
                new_chain_id = new_chain_id.lower()
            chainid_map[r.label_asym_id] = new_chain_id
            new_chain_ids.add(new_chain_id)

    # we filter H atoms here and only output alt_loc == A atoms
    atoms = atom_site[atom_site['label_asym_id'].isin(included_asym_ids) & (atom_site['type_symbol'] != 'H')].copy()
    atoms = atoms[atoms['label_alt_id'].isna() | (atoms['label_alt_id'] == 'A')]
    atoms['pdb_chain_id'] = atoms['label_asym_id'].apply(lambda a: chainid_map.get(a, ''))

    # select only model 1 (for NMR structures)
    atoms = atoms[atoms['pdbx_PDB_model_num'] == 1]

    # calculate and store closest atomic distance to ligand
    dist = cdist(atoms[['Cartn_x', 'Cartn_y', 'Cartn_z']].to_numpy(), atoms.loc[atoms['label_asym_id'] == ligand.asym_id, ['Cartn_x', 'Cartn_y', 'Cartn_z']].to_numpy())
    dist = np.min(dist, axis=1)
    # and store in occupancy field
    atoms['occupancy'] = dist

    # treat 5-letter comp_id codes for PDB
    # pdb comp_id codes 01-99 are reserved (https://www.wwpdb.org/news/news?year=2023#63ff72ccc031758bf1c30ff7)
    # we put a tilde at the front in addition
    legacy_comp_id_map = dict(zip(
        atoms.loc[atoms['label_comp_id'].str.len() > 3, 'label_comp_id'].unique(),
        (f'~{x:02d}' for x in it.count(1)), 
    ))
    atoms['legacy_comp_id'] = atoms['label_comp_id'].apply(lambda s: legacy_comp_id_map.get(s, s))
    
    fn = f'complex-strucs/{ligand.pdb_id}_{ligand.mon_id}_{ligand.asym_id}.complex.pdb'
    with open(fn, 'wt', encoding='utf-8') as f:
        for atom in atoms.itertuples():
            f.write(('{pdbrecname:<6s}{pdbserial:>5s} {pdbname:<4s}{pdbaltloc}{pdbresname:>3s} {chainid:1s}{pdbresid:>4s}    '
                     '{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}{occupancy:6.2f}{tempfactor:6.2f}      '
                     '{segid:<4s}{element:>2s}{charge:2s}\n').format(
                         pdbrecname=atom.group_PDB,
                         pdbserial=str(atom.id)[-5:],
                         pdbname=' '+atom.label_atom_id if (len(atom.type_symbol) == 1) and (len(atom.label_atom_id) <= 3) else atom.label_atom_id,
                         pdbaltloc=atom.label_alt_id if atom.label_alt_id is not None else ' ',
                         #pdbresname=atom.label_comp_id,
                         pdbresname=atom.legacy_comp_id,
                         chainid=atom.pdb_chain_id,
                         pdbresid=str(atom.auth_seq_id),
                         coords=[atom.Cartn_x, atom.Cartn_y, atom.Cartn_z],
                         occupancy=atom.occupancy,
                         tempfactor=atom.B_iso_or_equiv,
                         segid='',
                         element=atom.type_symbol,
                         charge='',
                     ))
    
    # now we write the same in mmcif
    # we need to write the chem_comp block for pdbe-arpeggio
    chem_comp = pd.DataFrame(j['chem_comp'])
    cif_chem_comp = chem_comp.copy()
    # convert numbers to string
    cif_chem_comp['formula_weight'] = cif_chem_comp['formula_weight'].apply(lambda a: '{:<.3f}'.format(a))
    # make sure everything is string
    cif_chem_comp = cif_chem_comp.astype('string')
    # fill missing values
    cif_chem_comp = cif_chem_comp.fillna({
        'mon_nstd_flag': '.',
        'pdbx_synonyms': '?',
    }).fillna('.')
    # rudimentary quoting
    cif_chem_comp = cif_chem_comp.map(lambda s: f'\n;{s};' if ('\n' in s) else (f'"{s}"' if '\'' in s else (f'\'{s}\'' if ((' ' in s) or ('"' in s) or ('[' in s) or (']' in s)) else s)))
    # justify left
    cols = chem_comp.columns
    cif_chem_comp_str = cif_chem_comp[cols].apply(ljust).to_string(index=False, header=False)
    cif_chem_comp_str = (
        'loop_\n' +
        '\n'.join([f'_chem_comp.{x}' for x in cols]) + '\n' +
        cif_chem_comp_str + '\n' +
        '#\n'
    )

    # write atom coordinates
    cif = atoms.copy()
    # convert numbers to string
    cif['Cartn_x'] = cif['Cartn_x'].apply(lambda a: '{:<.3f}'.format(a))
    cif['Cartn_y'] = cif['Cartn_y'].apply(lambda a: '{:<.3f}'.format(a))
    cif['Cartn_z'] = cif['Cartn_z'].apply(lambda a: '{:<.3f}'.format(a))
    cif['occupancy'] = cif['occupancy'].apply(lambda a: '{:<.2f}'.format(a))
    cif['B_iso_or_equiv'] = cif['B_iso_or_equiv'].apply(lambda a: '{:<.2f}'.format(a))
    cif = cif.astype({
        'pdbx_formal_charge': 'Int64',
        'label_seq_id': 'Int64',
    })
    # make sure everything is string
    cif = cif.astype('string')
    # fill missing values
    cif = cif.fillna({
        'label_alt_id': '.',
        'label_seq_id': '.',
    }).fillna('?')
    # rudimentary quoting
    cif = cif.map(lambda s: f'\n;{s};' if ('\n' in s) else (f'"{s}"' if '\'' in s else (f'\'{s}\'' if ((' ' in s) or ('"' in s) or ('[' in s) or (']' in s)) else s)))

    # justify left
    cols = atom_site.columns
    cif_str = cif[cols].apply(ljust).to_string(index=False, header=False)
    cif_str = (
        'loop_\n' +
        '\n'.join([f'_atom_site.{x}' for x in cols]) + '\n' +
        cif_str + '\n' +
        '#\n'
    )

    fn = f'complex-strucs/{ligand.pdb_id}_{ligand.mon_id}_{ligand.asym_id}.complex.cif'
    with open(fn, 'wt', encoding='utf-8') as f:
        f.write('data_' + str(ligand.pdb_id) + '\n')
        f.write('#\n')
        f.write(cif_chem_comp_str)
        f.write(cif_str)

    center = np.mean(atoms.loc[atoms['label_asym_id'] == ligand.asym_id, ['Cartn_x', 'Cartn_y', 'Cartn_z']].to_numpy(), axis=0)

    json_fn = f'complex-strucs/{ligand.pdb_id}_{ligand.mon_id}_{ligand.asym_id}.ligand_info.json'
    logging.info(f'writing ligand info {json_fn}')
    nonpoly_scheme = pd.DataFrame(j['pdbx_nonpoly_scheme'])
    pdb_seq_num = int(nonpoly_scheme.loc[nonpoly_scheme['asym_id'] == ligand.asym_id, 'pdb_seq_num'].unique()[0])
    pdb_strand_id = nonpoly_scheme.loc[nonpoly_scheme['asym_id'] == ligand.asym_id, 'pdb_strand_id'].unique()[0]

    d = {
        'pdb_id': ligand.pdb_id,
        'asym_id': ligand.asym_id,
        'pdb_strand_id': pdb_strand_id,
        'pdb_chain_id': chainid_map[ligand.asym_id],
        'pdb_seq_num': pdb_seq_num,
        'comp_id': ligand.mon_id,
        'center_x': center[0],
        'center_y': center[1],
        'center_z': center[2],
        'receptor_fn': fn,
        'chain_map': chainid_map,
        'comp_id_map': {v: k for k, v in legacy_comp_id_map.items()},
    }
    with open(json_fn, 'wt', encoding='utf-8') as f:
        json.dump(d, f)
    
    pymol_ligand_selector = f'{ligand.pdb_id}_{ligand.mon_id}_{ligand.asym_id}.complex and chain {chainid_map[ligand.asym_id]}'
    # some structures, like 7r2g, do not contain a nearby protein for all ligand instances
    if len(included_poly_entities) > 0:
        pymol_polymer_selector = f'{ligand.pdb_id}_{ligand.mon_id}_{ligand.asym_id}.complex and chain {chainid_map[list(included_poly_entities)[0]]}'
    else:
        pymol_polymer_selector = ''
    pymol_selectors.append((pymol_ligand_selector, pymol_polymer_selector))

for s in pymol_selectors:
    print('\t'.join(s))
