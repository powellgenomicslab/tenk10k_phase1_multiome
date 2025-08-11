import os
import anndata
import pandas as pd
import numpy as np
import networkx as nx
from matplotlib import rcParams
import matplotlib.pyplot as plt
from networkx.algorithms.bipartite import biadjacency_matrix
from networkx.drawing.nx_agraph import graphviz_layout

import scglue

PATH = "save/s04_compare"
os.makedirs(PATH, exist_ok=True)

# %%
genes = anndata.read_h5ad(
    "save/s01_preprocessing/omics_data/rna.h5ad", backed="r"
).var.query("dcq_highly_variable").index.to_numpy()

# %%
tfs = np.loadtxt("save/s03_infer_gene_tf/tfs.txt", dtype=str)

# %%
allATAC_glue_merged = nx.read_graphml("save/s03_infer_gene_tf/glue_merged.graphml.gz")
allATAC_glue_mat = biadjacency_matrix(allATAC_glue_merged, tfs, genes)
allATAC_glue_flag = allATAC_glue_mat.toarray().ravel()

# %%
partATAC_glue_merged = nx.read_graphml("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/scGLUE/src_celltype_eQTL_OneK1K/MonoC/save/s03_infer_gene_tf_corrected/glue_merged.graphml.gz")
partATAC_glue_mat = biadjacency_matrix(partATAC_glue_merged, tfs, genes)
partATAC_glue_flag = partATAC_glue_mat.toarray().ravel()

# %%
glue_crosstab = pd.crosstab(
    pd.Series(allATAC_glue_flag, name="allATAC_glue"),
    pd.Series(partATAC_glue_flag, name="partATAC_glue")
).iloc[::-1, ::-1]

print(glue_crosstab)


#glue_merged_csv = pd.read_csv('s03_infer_gene_tf/glue_merged.csv')
#glue_merged_GTEx_csv = pd.read_csv('../../Original_eQTL/s03_infer_gene_tf/glue_merged.csv')
#
#glue_merged_csv_all = glue_merged_csv.merge(glue_merged_GTEx_csv.drop_duplicates(), on=['TF', 'Target gene'], how='left', indicator=True)
#residual_glue = glue_merged_csv.loc[glue_merged_csv_all['_merge'] == 'left_only']
#
## %%
#g = nx.DiGraph()
#for ind, row in residual_glue.iterrows():
#    g.add_edge(row['TF'], row['Target gene'])
#nx.set_node_attributes(g, "target", name="type")
## nx.write_graphml(g, f"{PATH}/glue_merged.graphml.gz")
#
## %%
#scglue.plot.set_publication_params()
#rcParams['figure.figsize'] = (10, 10)
#plt.figure(figsize=(10, 10))
#nx.draw(g, graphviz_layout(g), with_labels=True)
#plt.savefig(f"{PATH}/residual_eQTL_graph.png", dpi=300)
#plt.close()
