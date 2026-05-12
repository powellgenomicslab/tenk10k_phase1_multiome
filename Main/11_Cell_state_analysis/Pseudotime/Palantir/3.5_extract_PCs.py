import sys
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import pandas as pd

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
sc.tl.pca(ad, svd_solver="arpack", n_comps=75)
print("Performing Harmony batch correction ...")
sce.pp.harmony_integrate(ad, "library")
ad.obsm["X_pca_raw"] = ad.obsm["X_pca"]
ad.obsm["X_pca"] = ad.obsm["X_pca_harmony"]

# Save PCA coordinates
print("Saving PCA coordinates ...")
pca_df = pd.DataFrame(
    ad.obsm["X_pca"],
    index=ad.obs_names
)
pca_df.to_csv("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CSVs/PCA_75PCs.csv")

print("Finished!")