#!/usr/bin/env python3
"""
UMAP Visualization of Epigenetic Age by Cell Type

This script creates UMAP visualizations colored by epigenetic age for each 
cell type to explore the relationship between cellular embedding and 
epigenetic aging patterns.

Author: Peter C Allen
"""

import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# Configuration
base_dir = "output/20250630_celltype_subsets"
fig_dir = "figures/4-epigeneticAge-Celltype"
os.makedirs(fig_dir, exist_ok=True)

print("Creating UMAP visualizations colored by epigenetic age...")

# Find available cell types
cell_types = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
print(f"Found {len(cell_types)} cell types to process")

for cell_type in cell_types:
    print(f"Processing {cell_type}...")
    
    ct_dir = os.path.join(base_dir, cell_type)
    
    # Define file paths
    file_paths = {
        'mtx': os.path.join(ct_dir, f"{cell_type}_subset_matrix.mtx"),
        'barcodes': os.path.join(ct_dir, f"{cell_type}_subset_rep1_barcodes.tsv"),
        'features': os.path.join(ct_dir, f"{cell_type}_subset_rep1_features.tsv"),
        'metadata': os.path.join(ct_dir, f"{cell_type}_subset_rep1_metadata.csv"),
        'epiage': os.path.join(base_dir, f"{cell_type}_epitrace_age.csv")
    }
    
    # Check that all required files exist
    missing_files = [name for name, path in file_paths.items() if not os.path.exists(path)]
    if missing_files:
        print(f"  Skipping {cell_type} - missing files: {', '.join(missing_files)}")
        continue

    # Load single-cell data
    adata = sc.read_mtx(file_paths['mtx'])
    # Add cell and feature names
    adata.obs_names = pd.read_csv(file_paths['barcodes'], header=None)[0].values
    adata.var_names = pd.read_csv(file_paths['features'], header=None)[0].values

    # Load and merge metadata
    metadata = pd.read_csv(file_paths['metadata'], index_col=0)
    adata.obs = adata.obs.join(metadata, how="left")

    # Load and merge epigenetic age estimates
    epiage_data = pd.read_csv(file_paths['epiage'], index_col=0)
    if "EpiTraceAge_iterative" not in epiage_data.columns:
        print(f"  Skipping {cell_type} - EpiTraceAge_iterative column not found")
        continue
    
    adata.obs = adata.obs.join(epiage_data[["EpiTraceAge_iterative"]], how="left")
    
    # Check if epigenetic age data was successfully merged
    n_with_epiage = adata.obs["EpiTraceAge_iterative"].notna().sum()
    print(f"  Epigenetic age available for {n_with_epiage}/{adata.n_obs} cells")

    # Perform standard preprocessing for UMAP
    print(f"  Computing UMAP for {adata.n_obs:,} cells...")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Create UMAP plot colored by epigenetic age
    plt.figure(figsize=(8, 6))
    sc.pl.umap(
        adata,
        color="EpiTraceAge_iterative",
        show=False,
        color_map="viridis",
        title=f"{cell_type} - Epigenetic Age",
        frameon=False,
        size=15
    )
    
    # Save the plot
    output_path = os.path.join(fig_dir, f"{cell_type}_epiage_umap.png")
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    
    print(f"  Saved UMAP to: {output_path}")

print("UMAP visualization completed!")