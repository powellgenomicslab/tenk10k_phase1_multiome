#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N chrnum_celltype
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=256GB
#PBS -l jobfs=100GB
#PBS -l ncpus=16
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/FineMapping/OutputLogs/2_finemapping/celltype/chrchrnum.OU
#PBS -e /g/data/fy54/od8037/TenK10K/FineMapping/OutputLogs/2_finemapping/celltype/chrchrnum.ER

cd /g/data/fy54/od8037/TenK10K/FineMapping

Rscript=/g/data/ei56/od8037/software/miniconda3/envs/r/bin/Rscript

touch Logs/2_finemapping/celltype/chrchrnum.log
> Logs/2_finemapping/celltype/chrchrnum.log

# Main command
/usr/bin/time -v $Rscript 2_finemapping.R celltype chrnum ${TMPDIR} >> Logs/2_finemapping/celltype/chrchrnum.log 2>&1

####


