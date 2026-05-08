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
#PBS -o /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/OutputLogs/3.5_combine_raw/celltype.OU
#PBS -e /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/OutputLogs/3.5_combine_raw/celltype.ER

CELLTYPE="celltype" 
BASE_DIR="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/Results/${CELLTYPE}"
OUT_DIR="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/FinalResults/SingleAssoc/${CELLTYPE}"

mkdir -p "$OUT_DIR"
mkdir -p "${BASE_DIR}/empty_tmp"

process_and_cleanup() {
    local CHR=$1
    local TARGET_DIR=$2
    local OUT_DIR=$3
    
    local CHR_DIR="${TARGET_DIR}/chr${CHR}"
    local OUT_FILE="${OUT_DIR}/chr${CHR}.tsv"
    local TAR_FILE="${TARGET_DIR}/chr${CHR}.tar.gz" 

    if [ ! -d "$CHR_DIR" ]; then
        echo "[Warn] Directory $CHR_DIR not found. Skipping."
        return
    fi

    echo "[$(date +%T)] Starting Chromosome $CHR..."

    # 1. Merge Data (PEAK prepended TSV)
    local FIRST_FILE=$(find "$CHR_DIR" -maxdepth 1 -name "*.singleAssoc.txt" -print -quit)
    if [ -n "$FIRST_FILE" ]; then
        printf "PEAK\t%s\n" "$(head -n 1 "$FIRST_FILE")" > "$OUT_FILE"
        
        # Run awk and capture the exit status
        find "$CHR_DIR" -maxdepth 1 -name "*.singleAssoc.txt" -exec awk '
            BEGIN { FS="\t"; OFS="\t"; }
            FNR > 1 {
                n = split(FILENAME, p, "/"); fname = p[n];
                sub(/\.singleAssoc\.txt$/, "", fname);
                print fname, $0
            }
        ' {} + >> "$OUT_FILE"
        
        # SAFEGUARD 1: Did the merge succeed?
        if [ $? -ne 0 ]; then
            echo "[ERROR] Awk merge failed for chr${CHR}. Aborting." >&2
            return 1
        fi
    fi

    # 2. Tar and compress the directory 
    echo "[$(date +%T)] Tarring and compressing chr${CHR}..."
    
    # SAFEGUARD 2: Did tar succeed AND is the archive uncorrupted?
    if tar -czf "$TAR_FILE" -C "$TARGET_DIR" "chr${CHR}" && gzip -t "$TAR_FILE"; then
        echo "[$(date +%T)] Archive verified. Wiping chr${CHR}..."
        
        # 3. Fast Delete using rsync
        rsync -a --delete "${TARGET_DIR}/empty_tmp/" "$CHR_DIR/"
        rmdir "$CHR_DIR"
        
        echo "[$(date +%T)] Finished chr${CHR}. Output: $OUT_FILE"
    else
        echo "[ERROR] Tar creation or verification failed for chr${CHR}! Aborting deletion." >&2
        return 1
    fi
}

export -f process_and_cleanup

module load parallel
printf "%s\n" {1..22} | parallel -j 22 process_and_cleanup {} "$BASE_DIR" "$OUT_DIR"

rmdir "${BASE_DIR}/empty_tmp"

echo "Finished"