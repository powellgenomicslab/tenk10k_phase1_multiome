#!/bin/bash
#PBS -P ke02
#PBS -q normalbw
#PBS -N chrchrnum
#PBS -l walltime=1:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=80GB
#PBS -l jobfs=50GB
#PBS -l ncpus=16
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/OutputLogs/0_convert_VCF/chrchrnum.OU
#PBS -e /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/OutputLogs/0_convert_VCF/chrchrnum.ER

CHR_NUM=chrnum
PREFIX="TenK10K_TOB_ATAC_renamed_chr${CHR_NUM}_all_variants"

TARGET_DIR="/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/Genotypes"
SOURCE_DIR="/g/data/fy54/angxue/brenner_backup/ScratchGeneral/proj/caQTL/genotype"

echo "Converting chromosome ${CHR_NUM}..."
/g/data/ei56/od8037/software/PLINK2/plink2 \
    --vcf "${SOURCE_DIR}/${PREFIX}.vcf.gz" \
    --make-bed \
    --threads 16 \
    --out "${TARGET_DIR}/${PREFIX}"

echo "Finished"
