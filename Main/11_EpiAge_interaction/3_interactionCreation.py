#!/usr/bin/env python3
"""
Epigenetic Age Interaction Data Creation

This script processes epigenetic age estimates computed for different cell types
and creates summary statistics (median and standard deviation) per donor for 
downstream interaction analysis with genotype data.

Author: Peter C Allen
"""

import anndata as ad
import os
import pandas as pd
import glob

# Load the merged dataset containing all cell types and metadata
adata = ad.read_h5ad('output/20250630_celltype_subsets/repeat1_merged_dataset.h5ad')

# Extract cell metadata
obs = adata.obs

# Save metadata as reference
print(f"Loaded metadata for {len(obs)} cells across {len(obs['predicted.id'].unique())} cell types")
obs.to_csv('output/20250630_celltype_subsets/repeat1_merged_dataset_meta.csv')

# Process each cell type's epigenetic age estimates
epitrace_files = glob.glob('output/20250630_celltype_subsets/*epitrace_age.csv')
print(f"Found {len(epitrace_files)} epigenetic age files to process")

for epitrace_file in epitrace_files:
    # Load epigenetic age data for this cell type
    epitrace_df = pd.read_csv(epitrace_file, index_col=0)

    # Extract cell type name from filename
    cell_type = '_'.join(os.path.basename(epitrace_file).split('_')[:-2])
    print(f"Processing {cell_type}...")

    # Merge epigenetic age data with cell metadata
    obs_celltype = epitrace_df.merge(obs, left_index=True, right_index=True, how='left')

    # Compute median epigenetic age per donor
    median_epitrace = obs_celltype.groupby('donor_id')['EpiTraceAge_iterative'].median().reset_index()
    median_epitrace.columns = ['individual', 'median_epiage']

    # Save median epigenetic age per donor
    output_path = f"data/20250630_epiages/{cell_type}_epiage.txt"
    median_epitrace.to_csv(output_path, sep='\t', index=False)
    print(f"  Saved median epigenetic age for {len(median_epitrace)} donors")

    # Compute standard deviation of epigenetic age per donor  
    sd_epitrace = obs_celltype.groupby('donor_id')['EpiTraceAge_iterative'].std().reset_index()
    sd_epitrace.columns = ['individual', 'sd_epiage']

    # Save standard deviation of epigenetic age per donor
    output_path = f"data/20250630_epiages/{cell_type}_epiageSD.txt"
    sd_epitrace.to_csv(output_path, sep='\t', index=False)
    print(f"  Saved epigenetic age variability for {len(sd_epitrace)} donors")

print("Epigenetic age interaction data creation completed.")
