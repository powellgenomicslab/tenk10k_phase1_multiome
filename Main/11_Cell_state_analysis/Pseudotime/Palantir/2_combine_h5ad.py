import anndata as ad
import os
import glob
import scanpy as sc

folder = "/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/Objects"
files = sorted(glob.glob(os.path.join(folder, "*.h5ad")))

def read_file(f):
    print(f"Reading {f}")
    return ad.read_h5ad(f)

print(f"Reading {len(files)} files")
ad_list = [read_file(f) for f in files]

print("Concatenating files")
ad_master = ad.concat(ad_list)

print("Filtering genes")
sc.pp.filter_genes(ad_master, min_counts=3)

print("Saving combined object")
ad_master.write_h5ad("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CD4_combined.h5ad", compression="gzip")

print("Finished!")
