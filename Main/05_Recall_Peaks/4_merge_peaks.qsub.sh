#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N lib_fragments
#PBS -l walltime=1:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=32GB
#PBS -l jobfs=50GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/ei56/od8037/TenK10K/PeakCalling/OutputLogs/4_merge_peaks
#PBS -e /g/data/ei56/od8037/TenK10K/PeakCalling/OutputLogs/4_merge_peaks

cd /g/data/ei56/od8037/TenK10K/PeakCalling/

Rscript=/g/data/ei56/od8037/software/miniconda3/envs/r/bin/Rscript

touch ./Logs/4_merge_peaks/merge_peaks.log
> ./Logs/4_merge_peaks/merge_peaks.log

# Main command
/usr/bin/time -v $Rscript 4_merge_peaks.R >> ./Logs/4_merge_peaks/merge_peaks.log 2>&1

####
