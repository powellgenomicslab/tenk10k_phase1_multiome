#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N celltype
#PBS -l walltime=2:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=64GB
#PBS -l jobfs=64GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/FineMapping/OutputLogs/1_process_caQTL/celltype.OU
#PBS -e /g/data/fy54/od8037/TenK10K/FineMapping/OutputLogs/1_process_caQTL/celltype.ER

cd /g/data/fy54/od8037/TenK10K/FineMapping

Rscript=/g/data/ei56/od8037/software/miniconda3/envs/r/bin/Rscript

touch Logs/1_process_caQTL/celltype.log
> Logs/1_process_caQTL/celltype.log

# Main command
/usr/bin/time -v $Rscript 1_process_caQTL.R celltype >> Logs/1_process_caQTL/celltype.log 2>&1

####


