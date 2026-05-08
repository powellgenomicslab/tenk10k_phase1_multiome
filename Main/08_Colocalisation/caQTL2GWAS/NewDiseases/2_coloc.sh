#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N ct_name_condition
#PBS -l walltime=2:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=16GB
#PBS -l jobfs=32GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/ColocNew/NewDiseases/OutputLogs/condition/ct_name.OU
#PBS -e /g/data/fy54/od8037/TenK10K/ColocNew/NewDiseases/OutputLogs/condition/ct_name.ER

cd /g/data/fy54/od8037/TenK10K/ColocNew/NewDiseases

Rscript=/g/data/ei56/od8037/software/miniconda3/envs/r/bin/Rscript

touch Logs/condition/ct_name.log
> Logs/condition/ct_name.log

# Main command
/usr/bin/time -v $Rscript 2_coloc.R ct_name condition >> Logs/condition/ct_name.log 2>&1

####
