
import os
import numpy as np
import pandas as pd
import scglue

import matplotlib.pyplot as plt
import seaborn as sns

PATH = "/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_smalltest/Treg_data_diff_eQTL/combined/s08_violin_plot"

conn_A_path = '/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_smalltest/Treg_data_diff_eQTL/combined/s04_infer_gene_tf/save/Treg/gene_peak_conn.pkl.gz'
conn_B_path = '/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_smalltest/Treg_data_diff_eQTL/combined/s04_infer_gene_tf/save/CD14_Mono/gene_peak_conn.pkl.gz'

conn_A = pd.read_pickle(conn_A_path)
conn_B = pd.read_pickle(conn_B_path)

matched_meanglue_eqtl = (conn_A['glue'][conn_A['eqtl']]).mean()
matched_meanglue_noeqtl = (conn_A['glue'][conn_A['eqtl'] == False]).mean()
unmatched_meanglue_eqtl = (conn_B['glue'][conn_A['eqtl']]).mean()
unmatched_meanglue_noeqtl = (conn_B['glue'][conn_A['eqtl'] == False]).mean()

conn_A['GLUE_with_CD14_Mono_eQTL'] = conn_B['glue']
conn_A['GLUE_with_Treg_eQTL'] = conn_A['glue']
conn_A['Peak-to-Gene pairs with support'] = conn_A['eqtl']

import matplotlib.pyplot as plt
import seaborn as sns

plot_df = conn_A.melt(
    id_vars="Peak-to-Gene pairs with support",
    value_vars=["GLUE_with_CD14_Mono_eQTL", "GLUE_with_Treg_eQTL"],
    var_name="Method",
    value_name="GLUE Score"
)
plot_df_true = plot_df[plot_df["Peak-to-Gene pairs with support"] == True]
fig, ax = plt.subplots(figsize=(7, 5.2))

custom_colors = {
    "GLUE_with_CD14_Mono_eQTL": "#17becf",  # dark cyan
    "GLUE_with_Treg_eQTL": "#e377c2"       # raspberry/magenta
}

sns.violinplot(data=plot_df_true, x="Method", y="GLUE Score", scale="width", palette=custom_colors)

ax.set_ylabel("Inferred Peak-to-Gene Regulatory Score", fontsize=16)
ax.set_xlabel("")
ax.tick_params(axis='both', labelsize=14)
ax.set_title("TenK10K Treg -- \nComparing Different eQTLs as Prior Guidance", fontsize=16)

plt.tight_layout()
fig.savefig(f"{PATH}/violin_plot.png", dpi=300, bbox_inches="tight")
