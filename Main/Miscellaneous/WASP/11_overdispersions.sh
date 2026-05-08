#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N 11_overdispersions
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l ncpus=28
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/11_overdispersions/11_overdispersions.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/11_overdispersions/11_overdispersions.ER

# Estimate overdispersions
source /g/data/ei56/od8037/software/miniconda3/etc/profile.d/conda.sh
conda activate wasp

export OUTDIR=/g/data/fy54/od8037/TenK10K/WASP_16/Outputs
fit_coefficients () {
    local CELLTYPE=${1}

    ls ${OUTDIR}/${CELLTYPE}/*.adjusted.hetp.txt.gz > ${OUTDIR}/${CELLTYPE}/cht_input_files.txt

    python /g/data/fy54/od8037/software/WASP/CHT/fit_as_coefficients.py \
	    ${OUTDIR}/${CELLTYPE}/cht_input_files.txt \
	    ${OUTDIR}/${CELLTYPE}/cht_as_coef.txt
    python /g/data/fy54/od8037/software/WASP/CHT/fit_bnb_coefficients.py \
	    --min_counts 50 \
	    --min_as_counts 10 \
        ${OUTDIR}/${CELLTYPE}/cht_input_files.txt \
	    ${OUTDIR}/${CELLTYPE}/cht_bnb_coef.txt

}
export -f fit_coefficients

module load parallel
parallel -j 28 --verbose \
    fit_coefficients {1} \
    :::: /g/data/fy54/od8037/TenK10K/WASP_16/celltype_names.txt

echo "Finished"
