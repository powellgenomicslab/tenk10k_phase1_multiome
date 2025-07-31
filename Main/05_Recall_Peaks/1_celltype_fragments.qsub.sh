#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N lib_fragments
#PBS -l walltime=2:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/ei56/od8037/TenK10K/PeakCalling/OutputLogs/1_celltype_fragments
#PBS -e /g/data/ei56/od8037/TenK10K/PeakCalling/OutputLogs/1_celltype_fragments

cd /g/data/ei56/od8037/TenK10K/PeakCalling/

Rscript=/g/data/ei56/od8037/software/miniconda3/envs/r/bin/Rscript

touch ./Logs/1_celltype_fragments/lib.log
> ./Logs/1_celltype_fragments/lib.log

# Main command
/usr/bin/time -v $Rscript 1_celltype_fragments.R lib >> ./Logs/1_celltype_fragments/lib.log 2>&1

####
