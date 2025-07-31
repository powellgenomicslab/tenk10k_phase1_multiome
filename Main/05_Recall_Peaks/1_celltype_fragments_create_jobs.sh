#!/bin/bash

# File containing library names (one per line)
LIBRARIES=/g/data/ei56/od8037/TenK10K/PeakCalling/library_names.txt

# Template script with the 'lib' placeholder
TEMPLATE=/g/data/ei56/od8037/TenK10K/PeakCalling/1_celltype_fragments.qsub.sh

# Output directory for the generated scripts
OUTPUT_DIR=/g/data/ei56/od8037/TenK10K/PeakCalling/Scripts/1_celltype_fragments

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through each line in the library list
while read -r LIBRARY; do
    # Replace 'lib' with the actual library name and write to new file
    sed "s/lib/$LIBRARY/g" "$TEMPLATE" > "$OUTPUT_DIR/${LIBRARY}.qsub.sh"
done < "$LIBRARIES"