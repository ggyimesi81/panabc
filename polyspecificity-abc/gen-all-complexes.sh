#!/bin/bash

for mon_id in 4C8 AIN FMM GBM H9F HT1 I6H MIX MTX R1H R1Q RHQ RS4 RTO STI SUZ TA1 TPF TTC URF VX7 X2N YMZ ; do
    python gen-complex-strucs.py --mon-id ${mon_id} ligand-entries/${mon_id}_ligands.tsv ligand-entries/${mon_id}.proximity.5.tsv
done

