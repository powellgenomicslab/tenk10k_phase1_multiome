#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N cell_type
#PBS -l walltime=4:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=55GB
#PBS -l jobfs=100GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -o /g/data/ei56/od8037/TenK10K/PeakCalling/OutputLogs/3_MACS3/cell_type.OU
#PBS -e /g/data/ei56/od8037/TenK10K/PeakCalling/OutputLogs/3_MACS3/cell_type.ER

cd /g/data/ei56/od8037/TenK10K/PeakCalling/

# Activate conda environment
source /g/data/ei56/od8037/software/miniconda3/etc/profile.d/conda.sh
conda activate py

touch ./Logs/3_MACS3/cell_type.log
> ./Logs/3_MACS3/cell_type.log

# Main command
/usr/bin/time -v python -u 3_MACS3.py cell_type >> ./Logs/3_MACS3/cell_type.log 2>&1

####