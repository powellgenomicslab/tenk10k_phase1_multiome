#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N computeEpiage-cd14Mono
#$ -o logs/epiAge_cd14Mono.out
#$ -e logs/epiAge_cd14Mono.err
#$ -l mem_requested=2000G
#$ -l h_rt=48:00:00

# Load necessary modules
source /home/petall/.bashrc

# Activate virtual environment
conda activate signac

# Change to project directory
cd /directflow/SCCGGroupShare/projects/petall/tenk10k_epiage_tensorQTL/

echo "Processing cell type: CD14_Mono"
touch logs/2_celltype_subsetEpiage/20250718-cd14Mono-epitrace.log

# Run R script
/usr/bin/time -v Rscript src/2_subset_fullTenk10k/2b_epitrace_classical_monocytes.R 2>&1 | tee "logs/2_celltype_subsetEpiage/20250718-cd14Mono-epitrace.log"