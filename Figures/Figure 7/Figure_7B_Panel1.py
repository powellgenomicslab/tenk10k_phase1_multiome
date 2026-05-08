
import os
import numpy as np
import pandas as pd
import scglue

import matplotlib.pyplot as plt
import seaborn as sns

PATH = "/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/10Xdata_diff_eQTL/combined/s08_violin_plot"
os.makedirs(PATH, exist_ok=True)

cell_type = "CD14_Mono"
conn_A_path = '/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/10Xdata_diff_eQTL/' + cell_type + '/save/s04_infer_gene_tf/gene_peak_conn.pkl.gz'
conn_B_path = '/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/10Xdata_diff_eQTL/GTEx_v8/save/s04_infer_gene_tf/gene_peak_conn.pkl.gz'
conn_A = pd.read_pickle(conn_A_path)
conn_B = pd.read_pickle(conn_B_path)
eqtl_source_pool = ['CD4_Naive', 'CD14_Mono', 'Treg']
for eqtl_source in eqtl_source_pool:
    eqtl_source = pd.read_pickle('/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/10Xdata_diff_eQTL/' + eqtl_source + '/save/s04_infer_gene_tf/gene_peak_conn.pkl.gz')['eqtl']
    conn_A['eqtl'] += eqtl_source
conn_A['GLUE_with_GTEx_eQTL'] = conn_B['glue']
conn_A['GLUE_with_TenK10K_eQTL'] = conn_A['glue']
conn_A['Peak-to-Gene pairs with support'] = conn_A['eqtl']

import matplotlib.pyplot as plt
import seaborn as sns

plot_df = conn_A.melt(
    id_vars="Peak-to-Gene pairs with support",
    value_vars=["GLUE_with_GTEx_eQTL", "GLUE_with_TenK10K_eQTL"],
    var_name="Method",
    value_name="GLUE Score"
)
plot_df_true = plot_df[plot_df["Peak-to-Gene pairs with support"] == True]
fig, ax = plt.subplots(figsize=(7, 5.2))

custom_colors = {
    "GLUE_with_GTEx_eQTL": "#FF6F61",       # a coral/red
    "GLUE_with_TenK10K_eQTL": "#4DB6AC"  # a teal/green
}

sns.violinplot(data=plot_df_true, x="Method", y="GLUE Score", scale="width", palette=custom_colors)

ax.set_ylabel("Inferred Peak-to-Gene Regulatory Score", fontsize=16)
ax.set_xlabel("")
ax.tick_params(axis='both', labelsize=14)
ax.set_title("10X Multiome -- \nComparing Different eQTLs as Prior Guidance", fontsize=16)
plt.tight_layout()
fig.savefig(f"{PATH}/violin_plot.png", dpi=300, bbox_inches="tight")