#!/bin/bash
# Convert SNP data to HDF5 format
module load bcftools
module load parallel

IN_DIR="/g/data/fy54/od8037/TenK10K/WASP/FromBrenner/VCFs"
OUT_DIR="/g/data/fy54/od8037/TenK10K/WASP_16/VCFs/ChrVCFs"
ls ${IN_DIR}/S0{220..221}.vcf.gz > ${OUT_DIR}/vcf_list.txt

# Merge VCF files
bcftools merge --file-list "/g/data/fy54/od8037/TenK10K/WASP_16/VCFs/ChrVCFs/vcf_list.txt" -Ou | \
    bcftools view -v snps -m2 -M2 -Oz -o "/g/data/fy54/od8037/TenK10K/WASP_16/VCFs/merged.snps.vcf.gz"
bcftools index -t "/g/data/fy54/od8037/TenK10K/WASP_16/VCFs/merged.snps.vcf.gz"

split_vcf() {
    local chr=${1}
    local input_vcf=${2}
    local outdir=${3}

    bcftools view -r chr${chr} ${input_vcf} -W -Oz -o "${outdir}/chr${chr}.vcf.gz"
    bcftools index -t "${outdir}/chr${chr}.vcf.gz"
}
export -f split_vcf

parallel -j 11 --verbose split_vcf {} "/g/data/fy54/od8037/TenK10K/WASP_16/VCFs/merged.snps.vcf.gz" "${OUT_DIR}" ::: {1..22}

OUTDIR=/g/data/fy54/od8037/TenK10K/WASP_16/HDF5
/g/data/fy54/od8037/software/WASP/snp2h5/snp2h5 \
    --chrom /g/data/fy54/od8037/TenK10K/WASP/autosome_length.txt \
    --format vcf \
    --geno_prob ${OUTDIR}/geno_probs.h5 \
    --haplotype ${OUTDIR}/haplotypes.h5 \
    --snp_index ${OUTDIR}/snp_index.h5 \
    --snp_tab ${OUTDIR}/snp_tab.h5 \
    /g/data/fy54/od8037/TenK10K/WASP_16/VCFs/ChrVCFs/chr*.vcf.gz

