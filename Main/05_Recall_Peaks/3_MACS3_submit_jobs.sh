#!/bin/bash

# File containing celltype names (one per line)
CELLTYPES=/g/data/ei56/od8037/TenK10K/PeakCalling/celltype_names.txt

# Script directory
SCRIPT_DIR=/g/data/ei56/od8037/TenK10K/PeakCalling/Scripts/3_MACS3

# Loop through each line in the celltype list and submit jobs
while read -r CELLTYPE; do
    # Submit the job
    qsub "$SCRIPT_DIR/${CELLTYPE}.qsub.sh"
done < "$CELLTYPES"
