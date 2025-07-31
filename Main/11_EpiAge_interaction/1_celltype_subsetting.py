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

# Merge them – inner join to keep only common genes, or outer to keep all
merged_adata = adatas[0].concatenate(adatas[1:], join='inner', batch_key='sample_pool')

# Optional: save the merged file
merged_adata.write("output/20250630_celltype_subsets/repeat1_merged_dataset.h5ad") # AnnData object with n_obs × n_vars = 3546117 × 440996

# merged_adata = load_and_merge_h5ad_files(directory) # AnnData object with n_obs × n_vars = 3546117 × 440996

merged_adata = sc.read_h5ad("output/20250630_celltype_subsets/repeat1_merged_dataset.h5ad")

np.random.seed(602)

total_cells_target = 250000
top_features = 200000
output_root = "output/20250630_celltype_subsets"

os.makedirs(output_root, exist_ok=True)

# Iterate through unique predicted.id values
for cell_type in merged_adata.obs['predicted.id'].unique():
    print(f"Processing {cell_type}...")

    # Filter for the current cell type
    ct_data = merged_adata[merged_adata.obs['predicted.id'] == cell_type].copy()
    n_cells = ct_data.n_obs

    if n_cells <= total_cells_target:
        print(f"  Using all {n_cells} cells (less than or equal to target).")
        subset_data = ct_data
    else:
        donor_ids = ct_data.obs['donor_id'].unique()
        cells_per_donor = total_cells_target // len(donor_ids)
        sampled_indices = []

        for donor in donor_ids:
            donor_data = ct_data[ct_data.obs['donor_id'] == donor]
            n = min(cells_per_donor, donor_data.n_obs)
            sampled = np.random.choice(donor_data.obs_names, n, replace=False)
            sampled_indices.extend(sampled)

        subset_data = ct_data[sampled_indices].copy()
        print(f"  Sampled {len(sampled_indices)} cells from {len(donor_ids)} donors.")

    # Subset to top highly variable genes/features
    sc.pp.highly_variable_genes(subset_data, n_top_genes=top_features, subset=True)
    print(f"  Selected top {top_features} variable features.")

    # Prepare output directory
    celltype_dir = os.path.join(output_root, cell_type.replace(" ", "_"))
    os.makedirs(celltype_dir, exist_ok=True)

    # Base filename prefix
    base_name = f"{cell_type}_subset"

    # Save MTX matrix
    scipy.io.mmwrite(os.path.join(celltype_dir, f"{base_name}_matrix.mtx"), subset_data.X.astype("float32"))

    # Save barcodes (cell names)
    pd.DataFrame(subset_data.obs.index).to_csv(os.path.join(celltype_dir, f"{base_name}_rep1_barcodes.tsv"), index=False, header=False)

    # Save features (gene or peak names)
    pd.DataFrame(subset_data.var.index).to_csv(os.path.join(celltype_dir, f"{base_name}_rep1_features.tsv"), index=False, header=False)

    # Save metadata
    subset_data.obs.to_csv(os.path.join(celltype_dir, f"{base_name}_rep1_metadata.csv"))

    print(f"  Saved MTX + metadata to {celltype_dir}\n")
