#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N 7_get_regions
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=256GB
#PBS -l jobfs=100GB
#PBS -l ncpus=28
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/7_get_regions/7_get_regions.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP_16/OutputLogs/7_get_regions/7_get_regions.ER
#PBS -M o.dong@garvan.org.au

# Identify target regions for CHT 
source /g/data/ei56/od8037/software/miniconda3/etc/profile.d/conda.sh
conda activate wasp

get_target_regions() {
    local CELLTYPE=${1}

    local ATACDIR=/g/data/fy54/od8037/TenK10K/WASP_16/Outputs/${CELLTYPE} 
    cd ${ATACDIR}
    local HDF5DIR=/g/data/fy54/od8037/TenK10K/WASP_16/HDF5

    ls ${ATACDIR}/read_counts.*.h5 | \
        xargs -n 1 basename | sed -E s/read_counts.// | \
        sed -E s/.h5// >  ${ATACDIR}/individuals.txt

    python /g/data/fy54/od8037/software/WASP/CHT/get_target_regions.py \
        --target_region_size 2000 \
        --min_as_count 10 \
        --min_read_count 100 \
        --min_het_count 1 \
        --min_minor_allele_count 1 \
        --chrom /g/data/fy54/od8037/TenK10K/WASP/autosome_length.txt \
        --read_count_dir ${ATACDIR} \
        --individuals ${ATACDIR}/individuals.txt \
        --samples /g/data/fy54/od8037/TenK10K/WASP_16/donor_names.txt \
        --snp_tab ${HDF5DIR}/snp_tab.h5 \
        --snp_index ${HDF5DIR}/snp_index.h5 \
        --haplotype ${HDF5DIR}/haplotypes.h5 \
        --output_file ${ATACDIR}/target_regions.txt.gz
}
export -f get_target_regions

module load parallel
parallel -j 28 --verbose \
    get_target_regions {1} \
    :::: /g/data/fy54/od8037/TenK10K/WASP_16/celltype_names.txt

echo "Finished"
