
import os
import numpy as np
import pandas as pd
import scglue

import matplotlib.pyplot as plt
import seaborn as sns

PATH = "save/s08_violin_plot"
os.makedirs(PATH, exist_ok=True)

# We can load in the saved 'gene_peak_conn.pkl.gz' file to pump in the save glue scores, and then make comparison.
conn_A_path = '/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_celltype_eQTL/CD14_Mono/save/s04_infer_gene_tf/gene_peak_conn_sig.pkl.gz'
conn_B_path = '/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_celltype_noeQTL/CD14_Mono/save/s04_infer_gene_tf/gene_peak_conn.pkl.gz'

conn_A = pd.read_pickle(conn_A_path)
conn_B = pd.read_pickle(conn_B_path)

conn_df = conn_A
conn_df['GLUE_with_eQTL_guidance'] = conn_A['glue']

min_GLUE, max_GLUE = conn_df['glue'].min(), conn_df['glue'].max()
GLUE_noeQTL_score = conn_B['glue']
GLUE_noeQTL_score_normed = (GLUE_noeQTL_score - GLUE_noeQTL_score.min()) / (GLUE_noeQTL_score.max() - GLUE_noeQTL_score.min()) * (max_GLUE - min_GLUE) + min_GLUE

conn_df['GLUE_base'] = GLUE_noeQTL_score_normed
conn_df['Peak-to-Gene pairs with support'] = conn_df['eqtl_pval'] < 5e-4

conn_df = conn_df[['GLUE_base', 'GLUE_with_eQTL_guidance', 'Peak-to-Gene pairs with support']]

import matplotlib.pyplot as plt
import seaborn as sns

plot_df = conn_df.melt(
    id_vars="Peak-to-Gene pairs with support",
    value_vars=[ "GLUE_base", "GLUE_with_eQTL_guidance"],
    var_name="Method",
    value_name="GLUE Score"
)
fig, ax = plt.subplots(figsize=(7, 5.2))

sns.violinplot(data=plot_df, x="Method", y="GLUE Score",
    hue="Peak-to-Gene pairs with support", scale="width")

ax.set_ylabel("Inferred Peak-to-Gene Regulatory Score", fontsize=16)
ax.set_xlabel("")
ax.tick_params(axis='both', labelsize=14)
ax.legend(title="Peak-to-Gene pairs with supporting evidence", fontsize=13, title_fontsize=14)
ax.set_title("TenK10K CD14_Mono", fontsize=16)
plt.tight_layout()
fig.savefig(f"{PATH}/violin_plot.png", dpi=300, bbox_inches="tight")

