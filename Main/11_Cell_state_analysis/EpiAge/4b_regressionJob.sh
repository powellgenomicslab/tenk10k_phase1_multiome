#!/bin/bash
#
# Genotype × Epigenetic Age Interaction Analysis Job Submission
#
# This script submits SGE jobs to perform genotype-epigenetic age interaction
# analysis for each cell type and chromosome combination.
#
# Author: Peter C Allen

echo "Submitting genotype × epigenetic age interaction analysis jobs..."

# Configuration
CELLTYPE_DIR="output/celltype_subsets"  # Update path as needed
LOG_DIR="logs/interaction_analysis"

# Create log directory
mkdir -p "$LOG_DIR"

# Check if cell type directory exists
if [ ! -d "$CELLTYPE_DIR" ]; then
    echo "Error: Cell type directory not found: $CELLTYPE_DIR"
    exit 1
fi

# Loop through all cell type subdirectories
for cell_path in "$CELLTYPE_DIR"/*; do
    if [ -d "$cell_path" ]; then
        cell_type=$(basename "$cell_path")
        echo "Submitting jobs for cell type: $cell_type"

        # Submit jobs for chromosomes 1-22
        for chr in {1..22}; do
            echo "  Submitting chr$chr..."
            
            # Submit SGE job
            qsub -N "interaction_${cell_type}_chr${chr}" <<EOF
#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -o $LOG_DIR/regression_${cell_type}_chr${chr}.log
#$ -N interaction_${cell_type}_chr${chr}
#$ -l mem_requested=100G
#$ -l h_rt=2:00:00

echo "Starting genotype × epigenetic age interaction analysis"
echo "Cell type: $cell_type"
echo "Chromosome: chr$chr"
echo "Job started at: \$(date)"

# Load conda environment
source ~/.bashrc
conda activate signac

# Run the regression analysis
Rscript src/2_subset_fullTenk10k/4_interactionRegression.R "$cell_type" "chr$chr"

echo "Job completed at: \$(date)"
EOF
        done
    fi
done

echo "All jobs submitted successfully!"
echo "Monitor job status with: qstat"

        done
    fi
done
