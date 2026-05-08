#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N cond_index
#PBS -l walltime=2:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=50GB
#PBS -l jobfs=32GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/ColocNew/DiseaseTraits/OutputLogs/1_process_GWAS/condition_index.OU
#PBS -e /g/data/fy54/od8037/TenK10K/ColocNew/DiseaseTraits/OutputLogs/1_process_GWAS/condition_index.ER

cd /g/data/fy54/od8037/TenK10K/ColocNew/DiseaseTraits

Rscript=/g/data/ei56/od8037/software/miniconda3/envs/r/bin/Rscript

touch Logs/1_process_GWAS/condition_index.log
> Logs/1_process_GWAS/condition_index.log

# Main command
/usr/bin/time -v $Rscript 1_process_GWAS.R index >> Logs/1_process_GWAS/condition_index.log 2>&1

####


