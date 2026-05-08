#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N 6_Ppartnum
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l ncpus=28
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/6_convert_BAM/Part_partnum.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/6_convert_BAM/Part_partnum.ER
#PBS -M o.dong@garvan.org.au

# Convert BAM files to HDF5 format
source /g/data/ei56/od8037/software/miniconda3/etc/profile.d/conda.sh
conda activate wasp

export HDF5DIR=/g/data/fy54/od8037/TenK10K/WASP_16/HDF5
export OUTDIR=/g/data/fy54/od8037/TenK10K/WASP_16/Outputs

convert_bam() {
    local POOL=${1}
    local CELLTYPE=${2}
    local DONOR=${3}

    mkdir -p ${OUTDIR}/${CELLTYPE}
    python /g/data/fy54/od8037/software/WASP/CHT/bam2h5.py \
        --chrom /g/data/fy54/od8037/TenK10K/WASP/autosome_length.txt \
        --snp_index ${HDF5DIR}/snp_index.h5 \
        --snp_tab ${HDF5DIR}/snp_tab.h5 \
        --haplotype ${HDF5DIR}/haplotypes.h5 \
        --individual ${DONOR} \
        --ref_as_counts ${OUTDIR}/${CELLTYPE}/ref_as_counts.${DONOR}.h5 \
        --alt_as_counts ${OUTDIR}/${CELLTYPE}/alt_as_counts.${DONOR}.h5 \
        --other_as_counts ${OUTDIR}/${CELLTYPE}/other_as_counts.${DONOR}.h5 \
        --read_counts ${OUTDIR}/${CELLTYPE}/read_counts.${DONOR}.h5 \
        /g/data/fy54/od8037/TenK10K/WASP/Outputs/FilteredBAMs/${CELLTYPE}-${DONOR}.bam
}
export -f convert_bam

module load parallel
parallel -j 28 --colsep '\t' --verbose \
    convert_bam {1} {2} {3} \
    :::: /g/data/fy54/od8037/TenK10K/WASP_16/Combinations/combinations_partnum.tsv

echo "Finished"
