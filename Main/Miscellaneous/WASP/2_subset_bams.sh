#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N poolname
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=64GB
#PBS -l jobfs=50GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP/OutputLogs/2_subset_bams/poolname.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP/OutputLogs/2_subset_bams/poolname.ER

source /g/data/ei56/od8037/software/miniconda3/etc/profile.d/conda.sh
conda activate py

POOL=poolname
SEARCH_ROOT="/g/data/fy54/data/atac/atac_count_outs/force_cells_test"
BAM=""
for dir in "${SEARCH_ROOT}"/*; do
    if [ -d "${dir}/${POOL}" ]; then
        BAM="${dir}/${POOL}/outs/possorted_bam.bam"
        break
    fi
done
OUTDIR=/g/data/fy54/od8037/TenK10K/WASP/BAMs
MAP_FILE=/g/data/fy54/od8037/TenK10K/WASP/CellMaps/${POOL}.tsv

echo "Starting single-pass split ..."
sinto filterbarcodes \
    --bam ${BAM} \
    --cells ${MAP_FILE} \
    --outdir ${OUTDIR} \
    --nproc 1

echo "Done"
