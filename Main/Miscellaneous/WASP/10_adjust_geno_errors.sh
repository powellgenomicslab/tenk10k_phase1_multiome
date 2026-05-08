#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N 10_Ppartnum
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l ncpus=28
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/10_adjust_geno_errors/Part_partnum.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/10_adjust_geno_errors/Part_partnum.ER

# Adjust to account for genotyping errors 
source /g/data/ei56/od8037/software/miniconda3/etc/profile.d/conda.sh
conda activate wasp

update_het_probs (){
    local POOL=${1}
    local CELLTYPE=${2}
    local DONOR=${3}

    local ATACDIR=/g/data/fy54/od8037/TenK10K/WASP_16/Outputs/${CELLTYPE} 

    python /g/data/fy54/od8037/software/WASP/CHT/update_het_probs.py \
        --ref_as_counts ${ATACDIR}/ref_as_counts.${DONOR}.h5  \
        --alt_as_counts ${ATACDIR}/alt_as_counts.${DONOR}.h5 \
        ${ATACDIR}/haplotype_read_counts.${DONOR}.adjusted.txt.gz \
        ${ATACDIR}/haplotype_read_counts.${DONOR}.adjusted.hetp.txt.gz
}
export -f update_het_probs

module load parallel
parallel -j 28 --colsep '\t' --verbose \
    update_het_probs {1} {2} {3} \
    :::: /g/data/fy54/od8037/TenK10K/WASP_16/Combinations/combinations_partnum.tsv

echo "Finished"
