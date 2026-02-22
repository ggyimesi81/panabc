#!/bin/bash

echo "pdb_id" >proximity_in.tsv
find test-pdb/mmcif -name "*.cif.gz" -printf "%f\n" | sed -e 's/\.cif\.gz$//' >>proximity_in.tsv
python proximity.py --in-tsv proximity_in.tsv --threshold 5.0 --tsv proximity.5.tsv

