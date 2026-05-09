import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad
import snapatac2 as snap
import os

# Get list of file paths
obj_dir = "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/ATAC_Final/output/New_Peak_scanpy"
all_files = sorted(os.listdir(obj_dir))
lib_list = [f for f in all_files if f.endswith("_1.h5ad")]
adata_dirs = [f"{obj_dir}/{x}" for x in lib_list]

adatas = []
for path in adata_dirs:
    print(f"Reading: {path}")
    try:
        adata = snap.read(path, backed=None)
        adata = adata[~adata.obs["predicted.id"].isin(["CD8_Proliferating", "ILC"])].copy()
        adatas.append(adata)
    except Exception as e:
        print(f"Failed to read {path}: {e}")

print(f"Successfully read {len(adatas)} files. Concatenating...")
adata_merge = ad.concat(adatas)
adata_merge.obs_names_make_unique()

# Find marker regions
print("Finding marker regions")
marker_peaks = snap.tl.marker_regions(adata_merge, groupby="predicted.id", pvalue=0.01)

# Compute motif enrichment
print("Computing motif enrichment")
motifs = snap.tl.motif_enrichment(
    motifs=snap.datasets.cis_bp(unique=True),
    regions=marker_peaks,
    genome_fasta=snap.genome.hg38
)

# Stack FC and FDR values across all cell types
fc = np.vstack([df['log2(fold change)'] for df in motifs.values()])
fdr = np.vstack([df['adjusted p-value'] for df in motifs.values()])

filter1 = np.apply_along_axis(lambda x: np.any(np.abs(x) >= 0), 0, fc)
filter2 = np.apply_along_axis(lambda x: np.any(x >= 0), 0, fdr)
passed = np.logical_and(filter1, filter2)

# Signed -ln(-log10(p))
sign = np.sign(fc[:, passed])
pvals = np.vstack([df['p-value'].to_numpy()[passed] for df in motifs.values()])
minval = np.min(pvals[np.nonzero(pvals)])
pvals = np.clip(pvals, minval, None)
pvals = sign * np.log(-np.log10(pvals))

# Create DataFrame of transformed scores
df = pd.DataFrame(
    pvals.T,
    columns=list(motifs.keys()),
    index=next(iter(motifs.values()))['id'].to_numpy()[passed],
)

# Define function to select top motifs for each cell type
def get_top_motif_ids(motifs_dict, cell_types, top_n=5):
    chosen_sets = []
    
    for cell_type in cell_types:
        if cell_type not in motifs_dict:
            print(f"Warning: {cell_type} not found in motifs.")
            continue
            
        # Keep it in Polars for speed
        df = motifs_dict[cell_type]
        
        # Filter p-value, sort by absolute log2FC, and take top N
        top_ids = (
            df.with_columns(abs_lfc = pl.col("log2(fold change)").abs())
            .sort(
                by=["p-value", "abs_lfc"], 
                descending=[False, True]  
            )
            .head(top_n)
            .get_column("id")
            .to_list()
        )
        chosen_sets.append(set(top_ids))
    
    # Combine all sets using union
    return set().union(*chosen_sets)

# target_cells = ["CD14_Mono", "CD16_Mono", "HSPC"]
target_cells = list(motifs.keys())
chosen_ids = get_top_motif_ids(motifs, target_cells, 5)

# Function to filter rows
def keep_row(row, target_cells=target_cells):
    max_celltype = row.idxmax()
    motif_id = row.name
    if max_celltype not in target_cells:
        return True
    return motif_id in chosen_ids

df_filtered = df[df.index.isin(chosen_ids)]
df_filtered.index = df_filtered.index.str.split("+").str[0]
df_filtered.columns = df_filtered.columns.str.replace("_", " ")

celltype_titles = [
    'ASDC', 'B<sub>intermediate</sub>', 'B<sub>memory</sub>', 'B<sub>naive</sub>',
    'CD4<sub>CTL</sub>', 'CD4<sub>Naive</sub>', 'CD4<sub>Proliferating</sub>',
    'CD4<sub>TCM</sub>', 'CD4<sub>TEM</sub>', 'CD8<sub>Naive</sub>',
    'CD8<sub>TCM</sub>', 'CD8<sub>TEM</sub>', 'CD14<sub>Mono</sub>',
    'CD16<sub>Mono</sub>', 'HSPC', 'MAIT', 'NK', 'NK<sub>CD56bright</sub>',
    'NK<sub>Proliferating</sub>', 'Plasmablast', 'T<sub>reg</sub>',
    'cDC1', 'cDC2', 'dnT', 'gdT', 'pDC'
]

celltype_order = [
    'CD4<sub>Naive</sub>', 'CD4<sub>TEM</sub>', 'CD4<sub>TCM</sub>', 'CD4<sub>CTL</sub>', 'CD4<sub>Proliferating</sub>', 'T<sub>reg</sub>',
    'gdT', 'dnT', 'MAIT',
    'CD8<sub>Naive</sub>', 'CD8<sub>TEM</sub>', 'CD8<sub>TCM</sub>', 
    'NK', 'NK<sub>CD56bright</sub>', 'NK<sub>Proliferating</sub>',
    'B<sub>naive</sub>', 'B<sub>intermediate</sub>', 'B<sub>memory</sub>', 'Plasmablast',
    'CD14<sub>Mono</sub>', 'CD16<sub>Mono</sub>',
    'cDC1', 'cDC2', 'pDC', 'ASDC',
    'HSPC'
]

df = df_filtered.T
df.index = celltype_titles
df = df.reindex(celltype_order)
df.index = [name.replace("<sub>", "$_{").replace("</sub>", "}$") for name in df.index]

g = sns.clustermap(
    df,
    row_cluster=False, 
    col_cluster=True,
    method="ward",
    cmap="RdBu_r",
    figsize=(19, 9),
    cbar_pos=(0.22, 1.05, 0.7, 0.03), # (left, bottom, width, height)
    cbar_kws={
        'label': r'$\ln(-\log_{10}(p))$',
        'orientation': 'horizontal' 
    } 
)

g.cax.tick_params(labelsize=14) 
g.cax.xaxis.label.set_size(16)

plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=12)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=12)

out_path = "/g/data/ei56/od8037/Plotting/Motif_Plot/motif_plot.png"
g.savefig(out_path, dpi=600, bbox_inches="tight")
