#!/bin/bash
#PBS -P ke02
#PBS -q normalbw
#PBS -N celltype
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l ncpus=22
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/OutputLogs/3.5_combine_raw/celltype.OU
#PBS -e /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/OutputLogs/3.5_combine_raw/celltype.ER

CELLTYPE="celltype" 
BASE_DIR="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Results/${CELLTYPE}"
OUT_DIR="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/FinalResults/SingleAssoc/${CELLTYPE}"

mkdir -p "$OUT_DIR"

process_chromosome() {
    local CHR=$1
    local TARGET_DIR=$2
    local OUT_DIR=$3
    
    local CHR_DIR="${TARGET_DIR}/chr${CHR}"
    local OUT_FILE="${OUT_DIR}/chr${CHR}.tsv"

    if [ ! -d "$CHR_DIR" ]; then
        echo "[Warn] Directory $CHR_DIR not found. Skipping."
        return
    fi

    echo "Processing chromosome $CHR..."

    # 1. Get header and prepend "PEAK" followed by a literal tab (\t)
    local FIRST_FILE=$(find "$CHR_DIR" -maxdepth 1 -name "*.singleAssoc.txt" -print -quit)
    if [ -z "$FIRST_FILE" ]; then
        return
    fi
    printf "PEAK\t%s\n" "$(head -n 1 "$FIRST_FILE")" > "$OUT_FILE"

    # 2. Process all files, extract PEAK, and prepend using Tab as the separator
    find "$CHR_DIR" -maxdepth 1 -name "*.singleAssoc.txt" -exec awk '
        BEGIN { 
            FS="\t";    # Input is tab-separated
            OFS="\t";   # Output should be tab-separated
        }
        FNR > 1 {
            # Extract filename from path
            n = split(FILENAME, path_parts, "/")
            fname = path_parts[n]
            
            # Strip the suffix
            sub(/\.singleAssoc\.txt$/, "", fname)
            
            # Print Peak ID then the original line
            print fname, $0
        }
    ' {} + >> "$OUT_FILE"

    echo "Finished chromosome $CHR. Output: $OUT_FILE"
}
export -f process_chromosome

module load parallel
printf "%s\n" {1..22} | parallel -j 22 process_chromosome {} "$BASE_DIR" "$OUT_DIR"

mkdir -p "${BASE_DIR}/empty_tmp"
process_and_cleanup() {
    local CHR=$1
    local TARGET_DIR=$2
    local OUT_DIR=$3
    
    local CHR_DIR="${TARGET_DIR}/chr${CHR}"
    local OUT_FILE="${OUT_DIR}/chr${CHR}.tsv"
    local TAR_FILE="${TARGET_DIR}/chr${CHR}.tar"

    if [ ! -d "$CHR_DIR" ]; then
        return
    fi

    echo "[$(date +%T)] Starting Chromosome $CHR..."

    # 1. Merge Data (PEAK prepended TSV)
    local FIRST_FILE=$(find "$CHR_DIR" -maxdepth 1 -name "*.singleAssoc.txt" -print -quit)
    if [ -n "$FIRST_FILE" ]; then
        printf "PEAK\t%s\n" "$(head -n 1 "$FIRST_FILE")" > "$OUT_FILE"
        find "$CHR_DIR" -maxdepth 1 -name "*.singleAssoc.txt" -exec awk '
            BEGIN { FS="\t"; OFS="\t"; }
            FNR > 1 {
                n = split(FILENAME, p, "/"); fname = p[n];
                sub(/\.singleAssoc\.txt$/, "", fname);
                print fname, $0
            }
        ' {} + >> "$OUT_FILE"
    fi

    # 2. Tar the directory (no compression)
    echo "[$(date +%T)] Tarring chr${CHR}..."
    tar -cf "$TAR_FILE" -C "$TARGET_DIR" "chr${CHR}"

    # 3. Fast Delete using rsync
    echo "[$(date +%T)] Wiping chr${CHR}..."
    rsync -a --delete "${TARGET_DIR}/empty_tmp/" "$CHR_DIR/"
    rmdir "$CHR_DIR"

    echo "[$(date +%T)] Finished chr${CHR}."
}
export -f process_and_cleanup

printf "%s\n" {1..22} | parallel -j 22 process_and_cleanup {} "$BASE_DIR" "$OUT_DIR"
rmdir "${BASE_DIR}/empty_tmp"

echo "Finished"
