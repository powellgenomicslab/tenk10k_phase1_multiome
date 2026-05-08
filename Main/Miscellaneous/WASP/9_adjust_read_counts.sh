#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N 9_adjust_read_counts
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=256GB
#PBS -l jobfs=100GB
#PBS -l ncpus=28
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/9_adjust_read_counts/9_adjust_read_counts.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/9_adjust_read_counts/9_adjust_read_counts.ER

# Adjust read counts
source /g/data/ei56/od8037/software/miniconda3/etc/profile.d/conda.sh
conda activate wasp

export OUTDIR=/g/data/fy54/od8037/TenK10K/WASP_16/Outputs
update_total_depth () {
    local CELLTYPE=${1}

    IN_FILE=${OUTDIR}/${CELLTYPE}/input_files.txt
    OUT_FILE=${OUTDIR}/${CELLTYPE}/output_files.txt
    ls ${OUTDIR}/${CELLTYPE}/haplotype_read_counts* | grep -v adjusted > $IN_FILE
    cat $IN_FILE | sed "s/.txt/.adjusted.txt/" >  $OUT_FILE

    python /g/data/fy54/od8037/software/WASP/CHT/update_total_depth.py \
        --seq /g/data/fy54/od8037/TenK10K/WASP/HDF5/seq.h5 $IN_FILE $OUT_FILE

}
export -f update_total_depth

module load parallel
parallel -j 28 --verbose \
    update_total_depth {1} \
    :::: /g/data/fy54/od8037/TenK10K/WASP_16/celltype_names.txt

echo "Finished"
