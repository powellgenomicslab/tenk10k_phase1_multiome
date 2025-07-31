#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N combine_fragments
#PBS -l walltime=2:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=8GB
#PBS -l jobfs=100GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/ei56/od8037/TenK10K/PeakCalling/OutputLogs/2_combine_fragments
#PBS -e /g/data/ei56/od8037/TenK10K/PeakCalling/OutputLogs/2_combine_fragments

# File containing celltype names (one per line)
CELLTYPES=/g/data/ei56/od8037/TenK10K/PeakCalling/celltype_names.txt

# Output directory for the combined fragments
OUTPUT_DIR=/g/data/ei56/od8037/TenK10K/PeakCalling/CelltypeFragments/AllCombinedFragments

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through each cell type
while read -r CELLTYPE; do
    INPUT_DIR="/g/data/ei56/od8037/TenK10K/PeakCalling/CelltypeFragments/${CELLTYPE}"
    OUTPUT_FILE="${OUTPUT_DIR}/${CELLTYPE}.tsv.gz"

    echo "Combining files for $CELLTYPE..."

    # Concatenate all .tsv.gz files and write to output
    if [ -d "$INPUT_DIR" ]; then
        cat "$INPUT_DIR"/*.tsv.gz > "$OUTPUT_FILE"
        echo "Done: $OUTPUT_FILE"
    else
        echo "Warning: Directory not found for $CELLTYPE"
    fi

done < "$CELLTYPES"