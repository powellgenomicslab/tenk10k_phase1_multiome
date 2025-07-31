#!/bin/bash

# File containing library names (one per line)
LIBRARIES=/g/data/ei56/od8037/TenK10K/PeakCalling/library_names.txt

# Script directory
SCRIPT_DIR=/g/data/ei56/od8037/TenK10K/PeakCalling/Scripts/1_celltype_fragments

# Loop through each line in the library list and submit jobs
while read -r LIBRARY; do
    # Submit the job
    qsub "$SCRIPT_DIR/${LIBRARY}.qsub.sh"
done < "$LIBRARIES"