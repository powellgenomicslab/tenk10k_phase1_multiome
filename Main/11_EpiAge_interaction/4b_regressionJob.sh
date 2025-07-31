#!/bin/bash

# Define the path to your cell type directories
CELLTYPE_DIR="output/2_celltype_subsets"

# Loop through all subdirectories (each representing a cell type)
for cell_path in "$CELLTYPE_DIR"/*; do
    if [ -d "$cell_path" ]; then
        cell_type=$(basename "$cell_path")

        # Loop through chromosomes 1 to 22
        for chr in {1..22}; do
            # Submit an SGE job
            qsub -N "job_${cell_type}_chr${chr}" <<EOF
#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -o logs/2_celltype_subsetEpiage/20250723-regression-${cell_type}-chr${chr}.log
#$ -N genotypeSNP_${cell_type}_chr${chr}
#$ -l mem_requested=100G
#$ -l h_rt=1:00:00


# Activate conda environment
source ~/.bashrc
conda activate signac

# Run the R script
Rscript src/2_subset_fullTenk10k/4_interactionRegression.R "$cell_type" chr"$chr"
EOF

        done
    fi
done
