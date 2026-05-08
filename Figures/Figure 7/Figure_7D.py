
import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
import utils

PATH = "save/s12_plot_genepeakscore_eqtl_pvalue"
os.makedirs(PATH, exist_ok=True)

# We can load in the saved 'gene_peak_conn.pkl.gz' file to pump in the save glue scores, and then make comparison.
conn_A_path = '/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_celltype_eQTL/CD14_Mono/save/s04_infer_gene_tf/gene_peak_conn_sig.pkl.gz'
conn_B_path = '/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_celltype_noeQTL/CD14_Mono/save/s04_infer_gene_tf/gene_peak_conn.pkl.gz'

conn_A = pd.read_pickle(conn_A_path)
conn_B = pd.read_pickle(conn_B_path)

conn_A['glue_noeQTL'] = conn_B['glue']
conn_A['logp_eqtl'] = -np.log10(conn_A['eqtl_pval'])

# %%
conn_A_eqtl = conn_A[conn_A['eqtl'] == True]

DIST_BINS = [0, 25, 50, 75, 100, 125, 150]
conn_A_eqtl["dist_bin"] = utils.make_dist_bins(conn_A_eqtl["dist"], bins=DIST_BINS)

conn_A_eqtl_glue = conn_A_eqtl.copy()
conn_A_eqtl_glue = conn_A_eqtl_glue.drop('glue_noeQTL', axis=1)
conn_A_eqtl_glue['if_eqtl'] = "with"

conn_A_eqtl_woglue = conn_A_eqtl.copy()
conn_A_eqtl_woglue['glue'] = conn_A_eqtl_woglue['glue_noeQTL']
conn_A_eqtl_woglue = conn_A_eqtl_woglue.drop('glue_noeQTL', axis=1)
conn_A_eqtl_woglue['if_eqtl'] = "w/o"

conn_A_eqtl_gain = conn_A_eqtl.copy()
conn_A_eqtl_gain['glue'] = conn_A_eqtl_gain['glue'] - conn_A_eqtl_gain['glue_noeQTL']
conn_A_eqtl_gain = conn_A_eqtl_gain.drop('glue_noeQTL', axis=1)
conn_A_eqtl_gain['if_eqtl'] = "gain"

# conn_A_eqtl_combined = pd.concat([conn_A_eqtl_glue, conn_A_eqtl_woglue, conn_A_eqtl_gain])
conn_A_eqtl_combined = pd.concat([conn_A_eqtl_glue, conn_A_eqtl_woglue])

plt.rcParams.update({
    'font.size': 14,           # General font size
    'axes.labelsize': 15,      # Axis labels
    'axes.titlesize': 15,      # Axis titles
    'legend.fontsize': 14,     # Legend
    'xtick.labelsize': 12,     # X-axis tick labels
    'ytick.labelsize': 12      # Y-axis tick labels
})
# g = utils.boxplot(x="dist_bin", y="glue", hue="if_eqtl", data=conn_A_eqtl_combined, hue_order=['w/o', 'with', 'gain'])
g = utils.boxplot(x="dist_bin", y="glue", hue="if_eqtl", data=conn_A_eqtl_combined, hue_order=['w/o', 'with'])
g.fig.set_size_inches(8, 6)
handles, labels = g.ax_joint.get_legend_handles_labels()
g.ax_joint.get_legend().remove()
g.ax_joint.set_xlabel("Genomic Distance")
g.ax_joint.set_ylabel("GLUE Regulatory Score")
for item in g.ax_joint.get_xticklabels():
    item.set_rotation(67.5)
legend = g.fig.legend(
    handles, labels,
    loc="center left",
    bbox_to_anchor=(0.72, 0.5),
    frameon=False,
    title="Using eQTL guidance\nfor improving\npeak-to-gene\nregulatory inference"
)
legend.get_title().set_multialignment('center')
plt.tight_layout()
g.fig.suptitle("TenK10K CD14_Mono", fontsize=15, y=0.92, x=0.4)
g.fig.subplots_adjust(right=0.77)  # leave room for legend on the right
g.fig.savefig(f"{PATH}/eqtl_dist_binned_glue_nogain.png", bbox_inches='tight', dpi=300)
