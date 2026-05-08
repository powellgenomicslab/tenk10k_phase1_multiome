#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N 8_Ppartnum
#PBS -l walltime=9:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l ncpus=28
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/8_make_CHT_input/Part_partnum.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/8_make_CHT_input/Part_partnum.ER
#PBS -M o.dong@garvan.org.au

# Create a CHT input file for each individual x celltype 
source /g/data/ei56/od8037/software/miniconda3/etc/profile.d/conda.sh
conda activate wasp

export HDF5DIR=/g/data/fy54/od8037/TenK10K/WASP_16/HDF5

extract_haplotype_read_counts(){
    local POOL=${1}
    local CELLTYPE=${2}
    local DONOR=${3}

    local ATACDIR=/g/data/fy54/od8037/TenK10K/WASP_16/Outputs/${CELLTYPE} 

    python /g/data/fy54/od8037/software/WASP/CHT/extract_haplotype_read_counts.py \
        --chrom /g/data/fy54/od8037/TenK10K/WASP/autosome_length.txt \
        --snp_index ${HDF5DIR}/snp_index.h5 \
        --snp_tab ${HDF5DIR}/snp_tab.h5 \
        --geno_prob ${HDF5DIR}/geno_probs.h5 \
        --haplotype ${HDF5DIR}/haplotypes.h5 \
        --samples /g/data/fy54/od8037/TenK10K/WASP_16/donor_names.txt  \
        --individual ${DONOR} \
        --ref_as_counts ${ATACDIR}/ref_as_counts.${DONOR}.h5 \
        --alt_as_counts ${ATACDIR}/alt_as_counts.${DONOR}.h5 \
        --other_as_counts ${ATACDIR}/other_as_counts.${DONOR}.h5 \
        --read_counts ${ATACDIR}/read_counts.${DONOR}.h5 \
        ${ATACDIR}/target_regions.txt.gz \
        | gzip > ${ATACDIR}/haplotype_read_counts.${DONOR}.txt.gz
}
export -f extract_haplotype_read_counts

module load parallel
parallel -j 28 --colsep '\t' --verbose \
    extract_haplotype_read_counts {1} {2} {3} \
    :::: /g/data/fy54/od8037/TenK10K/WASP_16/Combinations/combinations_partnum.tsv

echo "Finished"
