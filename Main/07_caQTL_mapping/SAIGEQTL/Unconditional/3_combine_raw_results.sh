#!/bin/bash
#PBS -P ke02
#PBS -q normalbw
#PBS -N 3_combine_results
#PBS -l walltime=4:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=120GB
#PBS -l jobfs=100GB
#PBS -l ncpus=26
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/OutputLogs/3_combine_results/3_combine_results.OU
#PBS -e /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/OutputLogs/3_combine_results/3_combine_results.ER

CELLTYPE_LIST=/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/celltype_names.txt

process_celltype() {
    local CELLTYPE=$1
    echo "Processing $CELLTYPE..."

    local OUT="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/FinalResults/Raw/${CELLTYPE}.csv"
    local TARGET_DIR="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Results/${CELLTYPE}"
    
    if [ ! -d "$TARGET_DIR" ]; then
        echo "[Warn] Directory not found for $CELLTYPE"
        return
    fi

    cd "$TARGET_DIR" || return
    rm -f "$OUT"

    local FIRST_FILE
    FIRST_FILE=$(find . -type f ! -name "*.*" ! -name "$OUT" -print -quit)
    
    if [ -z "$FIRST_FILE" ]; then
        echo "[Warn] No valid files found in $CELLTYPE. Skipping."
        return
    fi

    head -n 1 "$FIRST_FILE" > "$OUT"

    for i in {1..22}; do
        echo "Processing chromosome $i for $CELLTYPE..."
        if [ -d "chr$i" ]; then
            find "chr$i" -type f ! -name "*.*" \
                -exec awk 'FNR>1' {} + \
                >> "$OUT"
        fi
    done
}
export -f process_celltype

module load parallel
cat "$CELLTYPE_LIST" | parallel -j 26 process_celltype {}

echo "Finished"
