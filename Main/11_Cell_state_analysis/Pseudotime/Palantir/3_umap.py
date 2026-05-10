import sys
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import pandas as pd
import os 

# Capture command line arguments
inputArgs = sys.argv[1:]
num_pcs = int(inputArgs[0])
print(f"Number of PCs: {num_pcs}")

# Load data and normalise
print("Reading combined data...")
ad = sc.read_h5ad("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CD4_combined.h5ad")
print("Normalising...")
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)

# Identify highly variable genes, regress out effects and scale data
print("Finding highly variable genes, regressing out effects and scaling...")
sc.pp.highly_variable_genes(ad)
ad = ad[:, ad.var.highly_variable].copy()
sc.pp.regress_out(ad, ["nCount_RNA"])
sc.pp.scale(ad, max_value=10)

# Perform PCA
print("Performing PCA ...")
sc.tl.pca(ad, svd_solver="arpack", n_comps=num_pcs)
print("Performing Harmony batch correction ...")
sce.pp.harmony_integrate(ad, "library")
ad.obsm["X_pca_raw"] = ad.obsm["X_pca"]
ad.obsm["X_pca"] = ad.obsm["X_pca_harmony"]

# Generate PCA plots
plot_dir = f"/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/Plots/{num_pcs}PCs"
os.makedirs(plot_dir, exist_ok=True)
sc.settings.set_figure_params(dpi=1500, frameon=False, figsize=(5, 4))
with plt.rc_context():
    sc.pl.pca_variance_ratio(ad, log=True, n_pcs=num_pcs, show=False)
    plt.savefig(f"{plot_dir}/pca_variance_ratio_log.png", bbox_inches="tight")

with plt.rc_context():
    sc.pl.pca_variance_ratio(ad, log=False, n_pcs=num_pcs, show=False)
    plt.savefig(f"{plot_dir}/pca_variance_ratio.png", bbox_inches="tight")

with plt.rc_context():
    sc.pl.embedding(
        ad,
        basis="X_pca", 
        color="cell_type",
        frameon=False,
        show=False
    )
    plt.savefig(f"{plot_dir}/pca_plot.png", bbox_inches="tight")

# Compute UMAP
print("Generating KNN graph ...")
sc.pp.neighbors(ad, n_pcs=num_pcs, n_neighbors=30)
print("Computing UMAP ...")
sc.tl.umap(ad)

with plt.rc_context():
    sc.pl.embedding(
        ad,
        basis="X_umap",
        color="cell_type",
        frameon=False,
        show=False
    )
    plt.savefig(f"{plot_dir}/umap_plot.png", bbox_inches="tight")

# Save UMAP coordinates
print("Saving UMAP coordinates ...")
umap_df = pd.DataFrame(
    ad.obsm["X_umap"],
    index=ad.obs_names,
    columns=["UMAP1", "UMAP2"]
)

umap_dir = "/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/UMAP"
os.makedirs(umap_dir, exist_ok=True)
umap_df.to_csv(f"{umap_dir}/UMAP_PC{num_pcs}_KNN30.csv")

print("Finished!")
