#!/bin/bash
#PBS -P ke02
#PBS -q normalbw
#PBS -N 1_copy_index_VCF
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=50GB
#PBS -l ncpus=22
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/OutputLogs/1_index_VCF.OU
#PBS -e /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/OutputLogs/1_index_VCF.ER

module load parallel
module load bcftools

export TARGET_DIR="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/Genotypes"
export SOURCE_DIR="/g/data/fy54/angxue/brenner_backup/ScratchGeneral/proj/caQTL/genotype"
mkdir -p ${TARGET_DIR}

process_chr() {
    local CHR_NUM=$1
    local FILENAME="TenK10K_TOB_ATAC_renamed_chr${CHR_NUM}_all_variants.vcf.gz"
    
    echo "Copying chromosome ${CHR_NUM}..."
    cp "${SOURCE_DIR}/${FILENAME}" "${TARGET_DIR}/${FILENAME}"
    
    echo "Indexing chromosome ${CHR_NUM}..."
    bcftools index -f -c "${TARGET_DIR}/${FILENAME}"
}

export -f process_chr

parallel -j 22 process_chr ::: {1..22}
echo "Finished"
