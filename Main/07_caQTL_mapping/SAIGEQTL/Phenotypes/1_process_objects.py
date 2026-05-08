import snapatac2 as snap
import sys
import os
from tqdm.auto import tqdm

inputArgs = sys.argv[1:]
pool = inputArgs[0]
print(f"Processing pool: {pool}")

with open("/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/OldRun/celltype_names.txt") as f:
    celltype_list = [line.strip() for line in f if line.strip()]

print("Loading object")
obj_dir = "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/ATAC_Final/output/New_Peak_scanpy"
adata = snap.read(f"{obj_dir}/{pool}.h5ad", backed = None)

for celltype in tqdm(celltype_list):
    print(f"Processing celltype: {celltype}")
    with open(f"/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/ValidPeaks/{celltype}.txt") as f:
        valid_peaks = [line.strip() for line in f if line.strip()]
        
    adata_sub = adata[
        adata.obs["predicted.id"] == celltype, 
        adata.var_names.isin(valid_peaks)
    ]

    if adata_sub.n_obs == 0:
        print(f"No cells for {celltype}, skipping.")
        continue

    df = adata_sub.to_df()
    df.insert(0, "donor_id", adata_sub.obs["donor_id"].values)
    df.insert(1, "barcode", adata_sub.obs["library"].astype(str) + "_" + adata_sub.obs_names)
    df.columns = df.columns.str.replace("[:\\-]", "_", regex=True)

    print(f"Saving celltype: {celltype}")
    out_dir = f"/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Parquets/{celltype}"
    os.makedirs(out_dir, exist_ok = True)
    df.to_parquet(f"{out_dir}/{pool}.parquet", engine = "pyarrow", compression = "snappy")

print("Finished!")
