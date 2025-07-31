#!/bin/bash

# File containing celltype names (one per line)
CELLTYPES=/g/data/ei56/od8037/TenK10K/PeakCalling/celltype_names.txt

# Template script with the 'cell_type' placeholder
TEMPLATE=/g/data/ei56/od8037/TenK10K/PeakCalling/3_MACS3.qsub.sh

# Output directory for the generated scripts
OUTPUT_DIR=/g/data/ei56/od8037/TenK10K/PeakCalling/Scripts/3_MACS3

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through each line in the celltype list
while read -r CELLTYPE; do
    # Replace 'cell_type' with the actual celltype name and write to new file
    sed "s/cell_type/$CELLTYPE/g" "$TEMPLATE" > "$OUTPUT_DIR/${CELLTYPE}.qsub.sh"
done < "$CELLTYPES"