#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N atac_02_qc_pool_2
#PBS -l walltime=2:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=64GB
#PBS -l ncpus=1
#PBS -r y
#PBS -l wd
#PBS -M a.xue@garvan.org.au
#PBS -m ae

cd /g/data/ei56/ax3061/proj/tenk10k/caQTL/

Rscript=/g/data/ei56/ax3061/software/miniconda3/envs/multiome/bin/Rscript

touch ./output/logs/atac_02_qc_pool_2.log
> ./output/logs/atac_02_qc_pool_2.log

# Main command and record the time&mem usage
/usr/bin/time -v $Rscript atac_02_qc.R pool_2 >> ./output/logs/atac_02_qc_pool_2.log 2>&1


####
