# Load packages
import scanpy as sc
import palantir
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load data and normalise using Palantir
print("Reading combined data...")
ad = sc.read_h5ad("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CD4_combined.h5ad")
sc.pp.normalize_total(ad)
palantir.preprocess.log_transform(ad)

# Add PCA data
print("Loading PCA data...")
pca_df = pd.read_csv("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CSVs/PCA_75PCs.csv", index_col=0)
pca_df.columns = [f"PC{i+1}" for i in range(pca_df.shape[1])]
ad.obsm["X_pca"] = pca_df.loc[ad.obs_names].to_numpy()

# Add UMAP data
print("Loading UMAP data...")
umap_df = pd.read_csv("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/UMAP/UMAP_PC75_KNN30.csv", index_col=0)
ad.obsm["X_umap"] = umap_df.loc[ad.obs_names].to_numpy()

# Diffusion maps
print("Calculating diffusion maps...")
dm_res = palantir.utils.run_diffusion_maps(ad, n_components=100)

# Plot eigenvalues
print("Plotting eigenvalues...")
eigen_vals = ad.uns["DM_EigenValues"]
plt.plot(range(1, len(eigen_vals) + 1), eigen_vals, marker="o")
plt.xlabel("Diffusion Component")
plt.ylabel("Eigenvalue")
plt.title("Eigenvalue Decay")
plt.savefig("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/Plots/Palantir/eigenvalues.png", dpi=300, bbox_inches="tight")
plt.close()

# Plot derivative of eigenvalues
print("Plotting derivative of eigenvalues...")
ev = ad.uns["DM_EigenValues"]
plt.plot(np.diff(ev))
plt.title("First Derivative of Eigenvalues")
plt.savefig("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/Plots/Palantir/eigenvalues_derivative.png", dpi=300, bbox_inches="tight")
plt.close()

# Determine multiscale space
print("Determining multiscale space...")
ms_data = palantir.utils.determine_multiscale_space(ad)

# Visualise diffusion maps
print("Plotting diffusion components...")
palantir.plot.plot_diffusion_components(ad)
plt.savefig("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/Plots/Palantir/diffusion_components.png", dpi=300, bbox_inches="tight")
plt.close()

# Find early cell
early_cell = palantir.utils.early_cell(
    ad, 
    celltype = "CD4 Naive", 
    celltype_column = "cell_type", 
    eigvec_key = "DM_EigenVectors_multiscaled"
)
print(f"Early cell identified: {early_cell}")

# Find terminal states
terminal_states = palantir.utils.find_terminal_states(
    ad,
    celltypes=["Treg", "CD4 CTL"],
    celltype_column="cell_type",
    eigvec_key="DM_EigenVectors_multiscaled"
)
print(f"Terminal states identified: {terminal_states}")

# Save object
print("Saving AnnData object with diffusion maps...")
ad.write("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CD4_combined_diffusion.h5ad")
print("Finished!")
