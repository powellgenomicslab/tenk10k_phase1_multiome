#!/bin/bash
#PBS -P ke02
#PBS -q normalbw
#PBS -N celltype
#PBS -l walltime=4:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=64GB
#PBS -l jobfs=50GB
#PBS -l ncpus=22
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/OutputLogs/5_combine_results/celltype.OU
#PBS -e /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/OutputLogs/5_combine_results/celltype.ER

CELLTYPE="celltype"
TARGET_DIR="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/Results/${CELLTYPE}"
FINAL_OUT="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/FinalResults/Raw/${CELLTYPE}.csv"

echo "Processing $CELLTYPE..."

if [ ! -d "$TARGET_DIR" ]; then
    echo "[Warn] Directory not found for $CELLTYPE"
    exit 1
fi

cd "$TARGET_DIR" || exit 1

# Ensure the final output directory exists and clear old runs
mkdir -p "$(dirname "$FINAL_OUT")"
rm -f "$FINAL_OUT"

# Extract header from the first available peak file
FIRST_FILE=$(find . -mindepth 2 -type f ! -name "*.*" | head -n 1)
if [ -z "$FIRST_FILE" ]; then
    echo "[Warn] No valid files found in $CELLTYPE. Exiting."
    exit 1
fi
head -n 1 "$FIRST_FILE" > "$FINAL_OUT"

# Define function to process a single chromosome using the fast local JOBFS
process_chromosome() {
    local CHR=$1
    local TMP_CHR_OUT="${TMPDIR}/chr${CHR}_tmp.csv"
    
    if [ -d "chr${CHR}" ]; then
        find "chr${CHR}" -type f ! -name "*.*" -exec awk 'FNR>1' {} + > "$TMP_CHR_OUT"
    fi
}
export -f process_chromosome

# Parallelize across chromosomes 1 to 22
module load parallel
seq 1 22 | parallel -j 22 process_chromosome {}

# Concatenate the temporary JOBFS chromosome files into the final Lustre file in order
echo "Merging chromosomes..."
for i in {1..22}; do
    TMP_CHR_OUT="${TMPDIR}/chr${i}_tmp.csv"
    if [ -f "$TMP_CHR_OUT" ]; then
        cat "$TMP_CHR_OUT" >> "$FINAL_OUT"
    fi
done

echo "Finished $CELLTYPE"
