#!/bin/bash
# Build STAR reference genome 
/g/data/fy54/od8037/software/STAR-2.7.11b/STAR/bin/Linux_x86_64_static/STAR \
    --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir /g/data/fy54/od8037/TenK10K/WASP/Genome \
    --genomeFastaFiles /g/data/fy54/od8037/TenK10K/WASP/FromBrenner/genome.fa \
    --sjdbGTFfile /g/data/fy54/od8037/TenK10K/WASP/FromBrenner/genes.gtf \
    --genomeSuffixLengthMax 300 

# Split VCFs (from '/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/data/vcf/')
module load bcftools
module load parallel

IN_DIR="/g/data/fy54/od8037/TenK10K/WASP/FromBrenner/VCFs"
OUT_DIR="/g/data/fy54/od8037/TenK10K/WASP/VCFs/DonorVCFs/"

seq 220 227 | parallel -j 8 \
    "bcftools view -v snps -m2 -M2 $IN_DIR/S0{}.vcf.gz -Ou | \
    bcftools +split -O v -o $OUT_DIR"
