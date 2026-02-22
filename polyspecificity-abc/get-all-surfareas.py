#!/usr/bin/env python3

import logging
logging.basicConfig(level=logging.INFO)

import sys
import json
import glob
import os
import numpy as np
import pandas as pd

columns = [
    'C', 'A', 'N', 'NA', 'ND', 'OA', 'OB', 'OD', 'P', 'SA', 'SB', 'SD',
]

ligand_df = []
for fn in glob.glob('complex-strucs/surface-area/*.surfarea2.tsv'):
    ligand_key = os.path.basename(fn).split('.')[0]
    #fn = f'complex-strucs/surface-area/{ligand_key}.surfarea2.tsv'
    logging.info(fn)
    df = pd.read_csv(fn, dialect='excel-tab', na_filter=False)

    # normalize by ligand surface - use analytic solvent excluded surface
    fn = f'complex-strucs/surface-area/{ligand_key}_ligand.area'
    logging.info(fn)
    ligand_atomic_surf = pd.read_table(fn, sep='\s+')
    ligand_surface = ligand_atomic_surf['sas_0'].sum()

    df['ses'] = df['ses'] / ligand_surface
    df['sas'] = df['sas'] / ligand_surface
    ligand_df.append(df)

ligand_df = pd.concat(ligand_df)
ligand_df = ligand_df.pivot(index='protein_key', columns='atom_type', values='sas').fillna(0.0)
# rearrange columns
ligand_df = ligand_df[[x for x in columns if x in ligand_df.columns]]
"""
ligand_df = ligand_df.merge(long_table3_3.set_index('ligand_key')[['group', 'protein_name', 'mon_id']],
                            how='left', left_index=True, right_index=True)
"""

ligand_df.to_csv(sys.stdout, sep='\t')

