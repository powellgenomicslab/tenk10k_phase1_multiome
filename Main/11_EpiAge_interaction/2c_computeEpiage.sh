#!/bin/bash
#
# SGE Job Script for Epigenetic Age Computation
#
# This script submits array jobs to compute epigenetic age estimates for each
# cell type using the EpiTrace algorithm. Each task processes one cell type.
#
# Author: Peter C Allen
#
#$ -S /bin/bash
#$ -cwd
#$ -N computeEpiage
#$ -o logs/epiAge_$TASK_ID.out
#$ -e logs/epiAge_$TASK_ID.err
#$ -pe smp 4
#$ -l mem_requested=500G
#$ -l h_rt=48:00:00
#$ -t 1-28

echo "Starting epigenetic age computation task $SGE_TASK_ID"
echo "Job started at: $(date)"

# Load conda environment
source ~/.bashrc
conda activate signac

# Change to project directory (update path as needed)
cd /path/to/project/tenk10k_epiage_tensorQTL/

# Get list of available cell types
celltypes=($(ls -d output/20250630_celltype_subsets/*/ | xargs -n 1 basename))

echo "Available cell types: ${celltypes[@]}"
echo "Total cell types: ${#celltypes[@]}"

# Select cell type for this task
celltype=${celltypes[$((SGE_TASK_ID - 1))]}

# Validate task ID
if [ -z "$celltype" ]; then
    echo "Error: Invalid task ID $SGE_TASK_ID (no corresponding cell type)"
    exit 1
fi

echo "Processing cell type: $celltype"

# Create log directory for this cell type
mkdir -p logs/2_celltype_subsetEpiage
touch "logs/2_celltype_subsetEpiage/${celltype}.log"

# Run epigenetic age computation
echo "Starting EpiTrace analysis..."
/usr/bin/time -v Rscript src/2_subset_fullTenk10k/2_runEpiAge.R "$celltype" 2>&1 | tee "logs/2_celltype_subsetEpiage/${celltype}.log"

echo "Task completed at: $(date)"
echo "Cell type $celltype processing finished"
