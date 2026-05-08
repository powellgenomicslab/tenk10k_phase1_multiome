#!/bin/bash
#PBS -P ke02
#PBS -q normalbw
#PBS -N ct_name_chr_num
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=64GB
#PBS -l jobfs=50gb
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/MultiColoc/caQTL2GWAS/OutputLogs/condition/ct_name/chrchr_num.OU
#PBS -e /g/data/fy54/od8037/TenK10K/MultiColoc/caQTL2GWAS/OutputLogs/condition/ct_name/chrchr_num.ER

cd /g/data/fy54/od8037/TenK10K/MultiColoc/caQTL2GWAS

Rscript=/g/data/ei56/od8037/software/miniconda3/envs/r/bin/Rscript

touch Logs/condition/ct_name/chrchr_num.log
> Logs/condition/ct_name/chrchr_num.log

# Main command
/usr/bin/time -v $Rscript 3_coloc.R ct_name condition ${TMPDIR} chr_num >> Logs/condition/ct_name/chrchr_num.log 2>&1

####


