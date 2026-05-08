import scanpy as sc
import pandas as pd
import numpy as np
import os
from joblib import Parallel, delayed
from tqdm.auto import tqdm

def process_pool_batch(pools_subset, donor_to_idx, n_peaks, n_total_donors, obj_dir, target_celltypes):

    batch_results = {}
    for ct in target_celltypes:
        batch_results[ct] = {
            "counts": np.zeros(n_peaks, dtype=np.int64),
            "mask": np.zeros((n_peaks, n_total_donors), dtype=bool),
            "total_cells": 0,
            "seen_donors": set() 
        }

    for pool in pools_subset:
        try:
            adata = sc.read(f"{obj_dir}/{pool}.h5ad", backed='r')
            obs = adata.obs
            present_celltypes = obs['predicted.id'].unique()
            
            for ct in present_celltypes:     
                # Get subset for this celltype
                ct_mask = obs['predicted.id'] == ct
                ct_indices = np.where(ct_mask)[0]
                
                # Update total cells of this celltype
                batch_results[ct]["total_cells"] += len(ct_indices)
                
                # Identify which donors are involved in this celltype
                ct_donors = obs.loc[ct_mask, 'donor_id'].unique()
                
                # Iterate donors within this celltype
                for donor in ct_donors:
                    d_idx = donor_to_idx[donor]
                    batch_results[ct]["seen_donors"].add(d_idx)
                    
                    donor_ct_indices = np.where((obs["donor_id"] == donor) & (obs["predicted.id"] == ct))[0]
                    donor_matrix = adata.X[donor_ct_indices, :]

                    peak_counts = donor_matrix.sum(axis=0).A1.astype(np.int64)
                    batch_results[ct]["counts"] += peak_counts
                    
                    present_mask = peak_counts > 0
                    batch_results[ct]["mask"][present_mask, d_idx] = True
            print(pool)

        except Exception as e:
            print(f"Error processing {pool}: {e}")
            
    return batch_results

print("Loading metadata...")
with open("/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/donor_list.txt") as f:
    donor_list = [line.strip() for line in f if line.strip()]
donor_to_idx = {d: i for i, d in enumerate(donor_list)}
n_total_donors = len(donor_list)

obj_dir = "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/ATAC_Final/output/New_Peak_scanpy"
pools = sorted([f[:-5] for f in os.listdir(obj_dir) if f.endswith(".h5ad")])

ref_adata = sc.read(f"{obj_dir}/S0220_1.h5ad", backed='r')
all_peaks = ref_adata.var_names.to_numpy()
n_total_peaks = len(all_peaks)

with open("/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/OldRun/celltype_names.txt") as f:
    target_celltypes = [line.strip() for line in f if line.strip()]
print(f"Processing {len(target_celltypes)} cell types across {len(pools)} pools.")

# Run in parallel
pool_batches = np.array_split(pools, 16)
results = Parallel(n_jobs=16, verbose=50)(
    delayed(process_pool_batch)(
        batch, donor_to_idx, n_total_peaks, n_total_donors, obj_dir, target_celltypes
    ) for batch in pool_batches
)

print("Aggregating results by cell type...")
final_stats = {}
for ct in target_celltypes:
    final_stats[ct] = {
        "counts": np.zeros(n_total_peaks, dtype=np.int64),
        "mask": np.zeros((n_total_peaks, n_total_donors), dtype=bool),
        "total_cells": 0,
        "total_donors": set() 
    }
    
for batch_res in tqdm(results, desc="Merging Batches"):
    for ct, stats in batch_res.items():
        final_stats[ct]["counts"] += stats["counts"]
        final_stats[ct]["mask"] |= stats["mask"]
        final_stats[ct]["total_cells"] += stats["total_cells"]
        final_stats[ct]["total_donors"].update(stats["seen_donors"])
        
print("Filtering and saving...")
peak_series = pd.Series(all_peaks)
peak_chroms = peak_series.str.split(":").str[0]

for ct in tqdm(target_celltypes, desc="Cell Types"):
    ct_data = final_stats[ct]
    n_cells_ct = ct_data["total_cells"]
    n_donors_ct = len(ct_data["total_donors"])
    
    if n_cells_ct == 0:
        continue

    unique_donors_per_peak = ct_data["mask"].sum(axis=1)
    
    df = pd.DataFrame({
        "peak": all_peaks,
        "chrom": peak_chroms,
        "pct_cells": ct_data["counts"] / n_cells_ct,
        "pct_donors": unique_donors_per_peak / n_donors_ct
    })
    

    valid_df = df[
        (df["pct_cells"] > 0.01) & (df["pct_donors"] > 0.01)
    ]
    valid_df["peak"].to_csv(
        f"/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/ValidPeaks/{ct}.txt", 
        index=False, header=False
    )
