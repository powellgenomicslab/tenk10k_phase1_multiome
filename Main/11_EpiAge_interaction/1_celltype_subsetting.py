#!/usr/bin/env python3
"""
Cell Type Subsetting for Single-Cell ATAC-seq Analysis

This script processes merged single-cell ATAC-seq data and creates cell type-specific
subsets for downstream epigenetic age analysis. It performs donor-stratified sampling
to ensure balanced representation across individuals and selects highly variable peaks.

Author: Peter C Allen
"""

import scanpy as sc
import numpy as np
import pandas as pd
import os
import scipy.io

# Path to your directory of .h5ad files
data_dir = "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/ATAC_Final/output/New_Peak_scanpy/"
h5ad_files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if f.endswith('_1.h5ad')]

# Load all files
adatas = [sc.read_h5ad(f) for f in h5ad_files]

# Merge 
merged_adata = adatas[0].concatenate(adatas[1:], join='inner', batch_key='sample_pool')

merged_adata.write("output/20250630_celltype_subsets/repeat1_merged_dataset.h5ad") # AnnData object with n_obs × n_vars = 3546117 × 440996

np.random.seed(602)

# Configuration parameters
total_cells_target = 250000  # Maximum cells per cell type
top_features = 200000        # Number of most variable peaks to retain
output_root = "output/20250630_celltype_subsets"

# Create output directory
os.makedirs(output_root, exist_ok=True)

# Load the merged dataset (already processed and saved)
print(f"Dataset: {merged_adata.n_obs:,} cells × {merged_adata.n_vars:,} peaks")

# Process each cell type separately
cell_types = merged_adata.obs['predicted.id'].unique()
print(f"Processing {len(cell_types)} cell types...")

for cell_type in cell_types:
    print(f"\nProcessing {cell_type}...")

    # Extract cells for this cell type
    ct_data = merged_adata[merged_adata.obs['predicted.id'] == cell_type].copy()
    n_cells = ct_data.n_obs
    
    print(f"  Found {n_cells:,} cells")

    if n_cells <= total_cells_target:
        print(f"  Using all {n_cells:,} cells (below target)")
        subset_data = ct_data
    else:
        # Perform donor-stratified sampling to maintain representation balance
        donor_ids = ct_data.obs['donor_id'].unique()
        cells_per_donor = total_cells_target // len(donor_ids)
        sampled_indices = []

        print(f"  Sampling {cells_per_donor:,} cells per donor from {len(donor_ids)} donors")
        
        for donor in donor_ids:
            donor_data = ct_data[ct_data.obs['donor_id'] == donor]
            n = min(cells_per_donor, donor_data.n_obs)
            sampled = np.random.choice(donor_data.obs_names, n, replace=False)
            sampled_indices.extend(sampled)

        subset_data = ct_data[sampled_indices].copy()
        print(f"  Sampled {len(sampled_indices):,} cells total")

    # Select most variable peaks to reduce dimensionality
    sc.pp.highly_variable_genes(subset_data, n_top_genes=top_features, subset=True)
    print(f"  Selected {top_features:,} most variable peaks")

    # Save data in MTX format for downstream analysis
    celltype_dir = os.path.join(output_root, cell_type.replace(" ", "_"))
    os.makedirs(celltype_dir, exist_ok=True)

    base_name = f"{cell_type}_subset"

    # Save count matrix
    scipy.io.mmwrite(os.path.join(celltype_dir, f"{base_name}_matrix.mtx"), 
                     subset_data.X.astype("float32"))

    # Save cell barcodes
    pd.DataFrame(subset_data.obs.index).to_csv(
        os.path.join(celltype_dir, f"{base_name}_rep1_barcodes.tsv"), 
        index=False, header=False)

    # Save peak names
    pd.DataFrame(subset_data.var.index).to_csv(
        os.path.join(celltype_dir, f"{base_name}_rep1_features.tsv"), 
        index=False, header=False)

    # Save cell metadata
    subset_data.obs.to_csv(
        os.path.join(celltype_dir, f"{base_name}_rep1_metadata.csv"))

    print(f"  Saved to {celltype_dir}")

print("\nCell type subsetting completed successfully!")
