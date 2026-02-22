#!/usr/bin/env python

import numpy as np
import pandas as pd
from sklearn.svm import SVC, LinearSVC
import pickle as pkl

import logging
logging.basicConfig(level=logging.INFO)

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('surfarea_tsv', help='Surface areas in TSV format')
parser.add_argument('svm_model', help='SVM model file name')

args = parser.parse_args()

with open(args.svm_model, 'rb') as f:
    svc = pkl.load(f)

df = pd.read_csv(args.surfarea_tsv, dialect='excel-tab')

fields = [
    'A', 'C', 'N', 'NA', 'ND', 'OA', 'OB', 'P', 'SA', 'SB',
]
data = df.reindex(columns=fields).fillna(0.0)
svc_pred = svc.predict(data)

df['svc_pred'] = svc_pred

print(df.to_csv(sep='\t', index=False))

