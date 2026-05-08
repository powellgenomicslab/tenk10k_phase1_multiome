#!/bin/bash
#PBS -P ke02
#PBS -q normalbw
#PBS -N 12_run_CHT
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l ncpus=26
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/12_run_CHT/12_run_CHT.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/12_run_CHT/12_run_CHT.ER

# Estimate overdispersions
source /g/data/ei56/od8037/software/miniconda3/etc/profile.d/conda.sh
conda activate wasp

export OUTDIR=/g/data/fy54/od8037/TenK10K/WASP_16/Outputs
run_CHT () {
    local CELLTYPE=${1}

    python /g/data/fy54/od8037/software/WASP/CHT/combined_test.py \
        --min_as_counts 10 \
        --bnb_disp ${OUTDIR}/${CELLTYPE}/cht_bnb_coef.txt \
        --as_disp ${OUTDIR}/${CELLTYPE}/cht_as_coef.txt \
        ${OUTDIR}/${CELLTYPE}/cht_input_files.txt \
        ${OUTDIR}/${CELLTYPE}/cht_results.txt

}
export -f run_CHT

module load parallel
parallel -j 26 --verbose \
    run_CHT {1} \
    :::: /g/data/fy54/od8037/TenK10K/WASP_16/celltype_names.txt

echo "Finished"
