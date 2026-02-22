#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ligand_key')
parser.add_argument('ligand_pdb_fn')
parser.add_argument('--renumber', action='store_true', default=False, help='Renumber atom serials starting from 1')

args = parser.parse_args()

import logging
logging.basicConfig(level=logging.INFO)

import os
import sys

from PDB import PDBFile

import json

fn = f'../{args.ligand_key}.complex.pdb'
with open(fn, 'rt', encoding='utf-8') as f:
    pdbf = PDBFile()
    pdbf.parsepdb(f)

fn = f'../{args.ligand_key}.ligand_info.json'
with open(fn, 'rt', encoding='utf-8') as f:
    ligand_info = json.load(f)

# pdb chain ID in the complex.pdb file
ligand_chain_id = ligand_info['pdb_chain_id']
ligand_res_id = ligand_info['pdb_seq_num']
ligand_resname = ligand_info['comp_id']

ligand_atoms = [a for a in pdbf.struc.atoms if (
    (a.chainid == ligand_chain_id) and 
    (a.pdbresidnum == ligand_res_id) and 
    (a.resname == ligand_resname)
)]

res_names = list(sorted(set(a.resname for a in ligand_atoms)))

if (len(res_names) != 1) or (ligand_info['comp_id'] not in res_names):
    logging.error('ligand %s not found in %s, res_names = %s', ligand_info['comp_id'], args.ligand_key, res_names)
else:
    if args.renumber:
        s = 1
        for a in ligand_atoms:
            a.pdbserial = s
            s += 1
    with open(args.ligand_pdb_fn, 'wt', encoding='utf-8') as f:
        for a in ligand_atoms:
            f.write(a.formatpdb())

