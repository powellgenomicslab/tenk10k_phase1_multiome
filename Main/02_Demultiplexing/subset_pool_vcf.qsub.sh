#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N subset_TOB_vcf_by_ATAC_pools
#PBS -l walltime=1:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=10GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd

cd /g/data/ei56/projects/TenK10K/genotype/

VCF=/g/data/ei56/projects/TenK10K/genotype/Merged_MAF0.05_hg38_chr.vcf.gz

# i=${PBS_ARRAY_INDEX};
i=${var1}

pool="S0${i}"

touch ./logs/${pool}.log
echo ${pool} >> ./logs/${pool}.log
cat ${pool}_sample_list.txt >> ./logs/${pool}.log

/usr/bin/time -v bcftools view -S ${pool}_sample_list.txt -o ${pool}.vcf -Ov $VCF --threads 1 >> ./logs/${pool}.log 2>&1 

####
