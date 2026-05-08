#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N ct_name_chr_num
#PBS -l walltime=14:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=80GB
#PBS -l jobfs=32GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/ei56/od8037/NewGenotypes/Coloc/caQTL2eQTL/Output_Logs/ct_name_chr_num.OU
#PBS -e /g/data/ei56/od8037/NewGenotypes/Coloc/caQTL2eQTL/Output_Logs/ct_name_chr_num.ER

cd /g/data/ei56/od8037/NewGenotypes/Coloc/caQTL2eQTL

Rscript=/g/data/ei56/od8037/software/miniconda3/envs/r/bin/Rscript

touch Coloc_Logs/ct_name_chr_num.log
> Coloc_Logs/ct_name_chr_num.log

# Main command
/usr/bin/time -v $Rscript 1_coloc.R ct_name chr_num >> Coloc_Logs/ct_name_chr_num.log 2>&1

####


