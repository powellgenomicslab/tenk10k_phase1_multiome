#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N computeEpiage
#$ -o logs/epiAge_$TASK_ID.out
#$ -e logs/epiAge_$TASK_ID.err
#$ -pe smp 4
#$ -l mem_requested=500G
#$ -l h_rt=48:00:00
#$ -t 1-28

# Load necessary modules
source /home/petall/.bashrc

# Activate virtual environment
conda activate signac

# Change to project directory
cd /directflow/SCCGGroupShare/projects/petall/tenk10k_epiage_tensorQTL/

# Get list of cell types (optionally exclude CD14_Mono)
celltypes=($(ls -d output/20250630_celltype_subsets/*/ | xargs -n 1 basename | grep -v '^CD14_Mono$'))

# Print for debugging
echo "Available cell types: ${celltypes[@]}"

# Select current cell type
celltype=${celltypes[$((SGE_TASK_ID - 1))]}

# Check for invalid task ID
if [ -z "$celltype" ]; then
  echo "Error: No cell type found for SGE_TASK_ID=$SGE_TASK_ID"
  exit 1
fi

echo "Processing cell type: $celltype"
touch logs/2_celltype_subsetEpiage/${celltype}.log

# Run R script
/usr/bin/time -v Rscript src/2_subset_fullTenk10k/2_runEpiAge.R "$celltype" 2>&1 | tee "logs/2_celltype_subsetEpiage/${celltype}.log"
