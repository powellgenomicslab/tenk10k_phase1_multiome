# Load packages
import scanpy as sc
import pandas as pd

# Load object with diffusion maps
print("Reading combined data...")
ad = sc.read_h5ad("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CD4_combined_diffusion.h5ad")

# Extract eigenvectors
print("Extracting eigenvectors...")
eigen_vecs = ad.obsm["DM_EigenVectors"]

# Save eigenvectors
print("Saving eigenvectors...")
eigen_df = pd.DataFrame(eigen_vecs, index=ad.obs_names)
eigen_df.to_csv("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CSVs/DM_eigenvectors.csv")

# Extract multiscale space
print("Extracting multiscale space...")
ms_data = ad.obsm["DM_EigenVectors_multiscaled"]

# Save multiscale space
print("Saving multiscale space...")
ms_df = pd.DataFrame(ms_data, index=ad.obs_names)
ms_df.to_csv("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CSVs/DM_multiscale_space.csv")
print("Finished!")
