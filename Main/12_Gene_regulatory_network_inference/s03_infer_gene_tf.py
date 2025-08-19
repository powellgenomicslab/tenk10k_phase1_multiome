
import functools
import gc
import itertools
import operator
import os
from math import ceil
import subprocess

import anndata
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import scipy.stats
import seaborn as sns
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
from networkx.drawing.nx_agraph import graphviz_layout

import scglue
import utils

# %%
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (25, 25)

PATH = "save/s03_infer_gene_tf"
os.makedirs(PATH, exist_ok=True)

np.random.seed(0)

# %% [markdown]
# # Read data

# %%
rna = anndata.read_h5ad("save/s01_preprocessing/omics_data/rna.h5ad")
atac = anndata.read_h5ad("save/s01_preprocessing/omics_data/atac.h5ad")

rna.var["name"] = rna.var_names
atac.var["name"] = atac.var_names

# %%
genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names).query("dcq_highly_variable"))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names).query("dcq_highly_variable"))
tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)
flanks = tss.expand(500, 500)

# %%
dist_graph = nx.read_graphml("save/s01_preprocessing/network/dist.graphml.gz")  # Serves as genomic windows
eqtl_graph = nx.read_graphml("save/s01_preprocessing/network/eqtl.graphml.gz")

# %%
chip = scglue.genomics.read_bed("/directflow/SCCGGroupShare/projects/jayfan/Projects/scGLUE/DEMO_eQTL/metafile/ENCODE-TF-ChIP-hg38.bed.gz")
tfs = scglue.genomics.Bed(rna.var.loc[np.intersect1d(np.unique(chip["name"]), rna.var_names), :])
tfs.index.name = "tfs"

# %% [markdown]
# # SCENIC: coexpression network

# %%
# The following blocks only need to be executed at their first run.

### Already been run!
rna.obs['cells'] = rna.obs_names
rna[:, np.union1d(genes.index, tfs.index)].write_loom(f"{PATH}/rna.loom")
np.savetxt(f"{PATH}/tfs.txt", tfs.index, fmt="%s")

import subprocess

subprocess.run(f'singularity run --bind {PATH}:/data /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/scGLUE/aertslab-pyscenic-0.12.0.sif pyscenic grn /data/rna.loom /data/tfs.txt \
  -o /data/scenic_grn.csv --seed 0 --num_workers 20 --cell_id_attribute cells --gene_attribute features', shell = True, executable="/bin/bash")

# %%
### Already been run!

scenic_grn = pd.read_csv(f"{PATH}/scenic_grn.csv")
orphan_tfs = set(tfs.index).difference(genes.index)  # When treated as target genes cannot be included in cis-regulatory rankings
scenic_grn = scenic_grn.loc[[item not in orphan_tfs for item in scenic_grn["target"]], :]
scenic_grn.to_csv(f"{PATH}/scenic_grn.csv", index=False)

# %% [markdown]
# # Gene-peak connection

# %% [markdown]
# ## Distance

# %%
dist = biadjacency_matrix(dist_graph, genes.index, peaks.index, weight="dist", dtype=np.float32)

# %% [markdown]
# ## eQTL

# %%
eqtl = biadjacency_matrix(eqtl_graph, genes.index, peaks.index, weight="weight", dtype=np.float32)

# %% [markdown]
# ## GLUE

# %%
feature_embedding = pd.read_csv(f"save/s02_glue/ckpt/prior:dcq/feature_embeddings.csv", header=None, index_col=0)

## %%
glue_graph = scglue.genomics.regulatory_inference(
   feature_embedding.index, feature_embedding.to_numpy(), dist_graph.subgraph([*genes.index, *peaks.index]),
   alternative="greater", random_state=0
)
glue = biadjacency_matrix(glue_graph, genes.index, peaks.index, weight="score", dtype=np.float32)
pval = biadjacency_matrix(glue_graph, genes.index, peaks.index, weight="pval", dtype=np.float32)
qval = biadjacency_matrix(glue_graph, genes.index, peaks.index, weight="qval", dtype=np.float32)


## For plotting the gene tracks
## %%
#scglue.genomics.Bed(atac.var).write_bed(f"{PATH}/peaks.bed", ncols=3)
#scglue.genomics.write_links(
#   glue_graph,
#   scglue.genomics.Bed(rna.var).strand_specific_start_site(),
#   scglue.genomics.Bed(atac.var),
#   f"{PATH}/gene2peak.links", keep_attrs=["score"]
#)
#keep_attrs=["score"]
#gene2peak_link = nx.to_pandas_edgelist(
#        glue_graph
#    ).merge(
#        scglue.genomics.Bed(rna.var).strand_specific_start_site().df.iloc[:, :4], how="left", left_on="source", right_on="name"
#    ).merge(
#        scglue.genomics.Bed(atac.var).df.iloc[:, :4], how="left", left_on="target", right_on="name"
#    ).loc[:, [
#        "chrom_x", "chromStart_x", "chromEnd_x",
#        "chrom_y", "chromStart_y", "chromEnd_y",
#        *(keep_attrs or [])
#    ]].dropna()
#gene2peak_link["chromStart_x"] = gene2peak_link["chromStart_x"].astype('int')
#gene2peak_link["chromEnd_x"] = gene2peak_link["chromEnd_x"].astype('int')
#gene2peak_link["chromStart_y"] = gene2peak_link["chromStart_y"].astype('int')
#gene2peak_link["chromEnd_y"] = gene2peak_link["chromEnd_y"].astype('int')
#gene2peak_link.to_csv(f"{PATH}/gene2peak.links", sep="\t", index=False, header=False)
# gene2peak_link = pd.read_csv('/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/scGLUE/src_celltype_eQTL/CD14_Mono/save/s03_infer_gene_tf/gene2peak.links')
# subprocess.run('pyGenomeTracks --tracks tracks.ini --region chr1:629000-929000 --outFileName tracks.png', shell = True, executable="/bin/bash")


# Add eQTL info in the pygenomics plot
#gene_id_mapping = {ens: name for ens, name in zip(rna.var["gene_id_trimmed"], rna.var_names)}
#
## %%
#eqtl_tsv = pd.read_csv("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/scGLUE/src_celltype_eQTL_OneK1K/MonoC/save/s001_eQTL_summary/MonoC_cis_eqtls_210211.tsv", sep='\t', comment='#')
#eqtl_tsv = eqtl_tsv.loc[eqtl_tsv['chromStart'].str.isdigit()]
#eqtl_tsv = scglue.genomics.Bed(eqtl_tsv)
#
#eqtl_genes = eqtl_tsv['name']
#eqtl_genes = set(eqtl_genes)
#rna.var["in_eqtl"] = [item in eqtl_genes for item in rna.var_names]
## rna.var["in_eqtl"].sum()
#
## %%
#eqtl_links = eqtl_tsv.df.iloc[:, :4].merge(tss.df.iloc[:, :4], how="left", on="name").dropna().assign(score=1)
#eqtl_links["chromStart_y"] = eqtl_links["chromStart_y"].astype('int')
#eqtl_links["chromEnd_y"] = eqtl_links["chromEnd_y"].astype('int')
#eqtl_links = eqtl_links.query("chrom_x == chrom_y")
#eqtl_links["name"] = eqtl_links.pop("name")
#eqtl_links.to_csv(f"{PATH}/eqtl.annotated_links", sep="\t", index=False, header=False)
# subprocess.run("grep -w LINC01409 save/s03_infer_gene_tf/eqtl.annotated_links | cut -f1-7 > save/s03_infer_gene_tf/eqtl.links", shell = True, executable="/bin/bash")
# subprocess.run('pyGenomeTracks --tracks tracks_eqtl.ini --region chr1:629000-929000 --outFileName tracks_eqtl.png', shell = True, executable="/bin/bash")


# loc = rna.var.loc['RYR1']
# chrom = loc["chrom"]
# chromLen = loc["chromEnd"] - loc["chromStart"]
# chromStart = loc["chromStart"] - chromLen
# chromEnd = loc["chromEnd"] + chromLen
# REMEMBER!!! We need to ground the chromStart and End locations from float to integer.
# subprocess.run('pyGenomeTracks --tracks tracks.ini --region chr8:128245846-128724112 --outFileName tracks.png', shell = True, executable="/bin/bash")

# %% [markdown]
# ## Windowing

# %%
window = biadjacency_matrix(
   dist_graph, genes.index, peaks.index, weight=None, dtype=np.float32
).tocoo()

dist = window.multiply(dist.toarray())
eqtl = window.multiply(eqtl.toarray())
glue = window.multiply(glue.toarray())
pval = window.multiply(pval.toarray())
qval = window.multiply(qval.toarray())

for mat in (dist, eqtl, glue, qval):
   assert np.all(window.row == mat.row)
   assert np.all(window.col == mat.col)

# %%
gene_peak_conn = pd.DataFrame({
   "gene": genes.index[window.row],
   "peak": peaks.index[window.col],
   "dist": dist.data.astype(int),
   "eqtl": eqtl.data.astype(bool),
   "glue": glue.data,
   "pval": pval.data,
   "qval": qval.data
})
# gene_peak_conn["pchic"] = pd.Categorical(gene_peak_conn["pchic"], categories=[False, True])
# gene_peak_conn["eqtl"] = pd.Categorical(gene_peak_conn["eqtl"], categories=[False, True])

#### Already been run!
gene_peak_conn.to_pickle(f"{PATH}/gene_peak_conn.pkl.gz")
# gene_peak_conn = pd.read_pickle(f"{PATH}/gene_peak_conn.pkl.gz")


# %% [markdown]
# # Filtering gene-peak connection

# %% [markdown]
# ## GLUE

# %%
qval_cutoff = 0.05
_ = sns.scatterplot(
    x="glue", y="qval", data=gene_peak_conn.sample(n=2000),
    edgecolor=None, s=5, alpha=0.1
).axhline(y=qval_cutoff, c="darkred", ls="--")

# %%
# if gene_peak_conn.query(f'qval < {qval_cutoff}').shape[0] > 0:
#     gene_peak_conn_glue = gene_peak_conn.query(f"qval < {qval_cutoff}")
# elif gene_peak_conn.query('pval < 0.05').shape[0] > 5:
#     gene_peak_conn_glue = gene_peak_conn.query('pval < 0.05')
# else:
#     gene_peak_conn_glue = gene_peak_conn.query('pval < 0.1')
gene_peak_conn_glue = gene_peak_conn.query(f"qval < {qval_cutoff}")
# gene_peak_conn_glue.shape[0]
print('---------------')
print('Number of gene_peak pairs with GLUE qval < 0.05: ' + str(gene_peak_conn_glue.shape[0]))
print('Number of gene_peak pairs with GLUE qval < 0.1: ' + str(gene_peak_conn.query("qval < 0.1").shape[0]))
print('---------------')

# %%
frac_pos = gene_peak_conn_glue.shape[0] / gene_peak_conn.shape[0]
# frac_pos

# %%
glue_cutoff = gene_peak_conn.query(f"qval < {qval_cutoff}")["glue"].min()
# glue_cutoff

# %%
### Already been run!
#glue_all_links = gene_peak_conn.loc[:, ["gene", "peak", "glue"]].merge(
#    tss.df.iloc[:, :4], how="left", left_on="gene", right_index=True
#).merge(
#    peaks.df.iloc[:, :4], how="left", left_on="peak", right_index=True
#).loc[:, [
#    "chrom_x", "chromStart_x", "chromEnd_x",
#    "chrom_y", "chromStart_y", "chromEnd_y",
#    "glue", "gene"
#]]
#glue_all_links.to_csv(f"{PATH}/glue_all.annotated_links", sep="\t", index=False, header=False)
## del glue_all_links
#
## %%
#### Already been run!
#glue_links = gene_peak_conn_glue.loc[:, ["gene", "peak", "glue"]].merge(
#    tss.df.iloc[:, :4], how="left", left_on="gene", right_index=True
#).merge(
#    peaks.df.iloc[:, :4], how="left", left_on="peak", right_index=True
#).loc[:, [
#    "chrom_x", "chromStart_x", "chromEnd_x",
#    "chrom_y", "chromStart_y", "chromEnd_y",
#    "glue", "gene"
#]]
#glue_links.to_csv(f"{PATH}/glue.annotated_links", sep="\t", index=False, header=False)
# del glue_links

# %% [markdown]
# ## Distance

# %%
dist_cutoff = np.quantile(gene_peak_conn["dist"], frac_pos)
# dist_cutoff

# %%
gene_peak_conn_dist = gene_peak_conn.query(f"dist <= {dist_cutoff}")
# gene_peak_conn_dist.shape[0]

# %% [markdown]
# ## eQTL

# %%
gene_peak_conn_eqtl = gene_peak_conn.query("eqtl")
# gene_peak_conn_eqtl.shape[0]

# %% [markdown]
# # TF binding

# %% [markdown]
# ## Flanks

# %%
### Already been run!
flank_tf_binding = scglue.genomics.window_graph(flanks, chip, 0, right_sorted=True)
flank_tf_binding = nx.to_pandas_edgelist(flank_tf_binding, source="flank", target="tf")
# flank_tf_binding.shape

# %%
s = set(tfs.index)
flank_tf_binding = flank_tf_binding.loc[[item in s for item in flank_tf_binding["tf"]], :]
# flank_tf_binding.shape

# %%
flank_tf_binding.to_pickle(f"{PATH}/flank_tf_binding.pkl.gz")

# flank_tf_binding = pd.read_pickle(f"{PATH}/flank_tf_binding.pkl.gz")

# %% [markdown]
# ## Peaks

# %%
### Already been run!
peak_tf_binding = scglue.genomics.window_graph(peaks, chip, 0, right_sorted=True)
peak_tf_binding = nx.to_pandas_edgelist(peak_tf_binding, source="peak", target="tf")
# peak_tf_binding.shape

# %%
s = set(tfs.index)
peak_tf_binding = peak_tf_binding.loc[[item in s for item in peak_tf_binding["tf"]], :]
# peak_tf_binding.shape

# %%
peak_tf_binding.to_pickle(f"{PATH}/peak_tf_binding.pkl.gz")

# peak_tf_binding = pd.read_pickle(f"{PATH}/peak_tf_binding.pkl.gz")

# %% [markdown]
# # Cis-regulatory ranking

# %% [markdown]
# ## Flank

# %%
observed_flank_tf = scipy.sparse.coo_matrix((
    np.ones(flank_tf_binding.shape[0], dtype=np.int16), (
        flanks.index.get_indexer(flank_tf_binding["flank"]),
        tfs.index.get_indexer(flank_tf_binding["tf"]),
    )
), shape=(flanks.index.size, tfs.index.size)).toarray()

# %%
rank_flank_tf = pd.DataFrame(
    scipy.stats.rankdata(-observed_flank_tf, axis=0),
    index=flanks.index, columns=tfs.index
)
# rank_flank_tf.iloc[:5, :5]

# %% [markdown]
# ## Distance

# %%
### Already been run!
enrichment_gene_tf_dist, rank_gene_tf_dist = utils.cis_regulatory_ranking(
   gene_peak_conn_dist, peak_tf_binding,
   genes, peaks, tfs, n_samples=1000, random_seed=0
)

# %%
enrichment_gene_tf_dist.to_pickle(f"{PATH}/enrichment_gene_tf_dist.pkl.gz")
rank_gene_tf_dist.to_pickle(f"{PATH}/rank_gene_tf_dist.pkl.gz")

# enrichment_gene_tf_dist = pd.read_pickle(f"{PATH}/enrichment_gene_tf_dist.pkl.gz")
# rank_gene_tf_dist = pd.read_pickle(f"{PATH}/rank_gene_tf_dist.pkl.gz")

# %% [markdown]
# ## eQTL

# %%
### Already been run!
enrichment_gene_tf_eqtl, rank_gene_tf_eqtl = utils.cis_regulatory_ranking(
   gene_peak_conn_eqtl, peak_tf_binding,
   genes, peaks, tfs, n_samples=1000, random_seed=0
)

# %%
enrichment_gene_tf_eqtl.to_pickle(f"{PATH}/enrichment_gene_tf_eqtl.pkl.gz")
rank_gene_tf_eqtl.to_pickle(f"{PATH}/rank_gene_tf_eqtl.pkl.gz")

# enrichment_gene_tf_eqtl = pd.read_pickle(f"{PATH}/enrichment_gene_tf_eqtl.pkl.gz")
# rank_gene_tf_eqtl = pd.read_pickle(f"{PATH}/rank_gene_tf_eqtl.pkl.gz")

# %% [markdown]
# ## GLUE

# %%
### Already been run!
enrichment_gene_tf_glue, rank_gene_tf_glue = utils.cis_regulatory_ranking(
   gene_peak_conn_glue, peak_tf_binding,
   genes, peaks, tfs, n_samples=1000, random_seed=0
)

# %%
enrichment_gene_tf_glue.to_pickle(f"{PATH}/enrichment_gene_tf_glue.pkl.gz")
rank_gene_tf_glue.to_pickle(f"{PATH}/rank_gene_tf_glue.pkl.gz")

# enrichment_gene_tf_glue = pd.read_pickle(f"{PATH}/enrichment_gene_tf_glue.pkl.gz")
# rank_gene_tf_glue = pd.read_pickle(f"{PATH}/rank_gene_tf_glue.pkl.gz")

# %% [markdown]
# # SCENIC: cisTarget pruning

# %%
ctx_annotation = pd.concat([
    pd.DataFrame({
        "#motif_id": tfs.index + "_atac",
        "gene_name": tfs.index
    }),
    pd.DataFrame({
        "#motif_id": tfs.index + "_flank",
        "gene_name": tfs.index
    })
]).assign(
    motif_similarity_qvalue=0.0,
    orthologous_identity=1.0,
    description="placeholder"
)
ctx_annotation.to_csv(f"{PATH}/ctx_annotation.tsv", sep="\t", index=False)

# %%
flank_feather = rank_flank_tf.T
flank_feather = flank_feather.loc[np.unique(flank_feather.index), np.unique(flank_feather.columns)].astype(np.int16)
flank_feather.index += "_flank"
flank_feather.index.name = "tracks"
flank_feather.columns.name = None
columns = flank_feather.columns.tolist()
flank_feather = flank_feather.reset_index()
flank_feather = flank_feather.loc[:, [*columns, "tracks"]]
flank_feather.to_feather(f"{PATH}/flank_ctx.genes_vs_tracks.rankings.feather")

subprocess.run(f'singularity run --bind {PATH}:/data /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/scGLUE/aertslab-pyscenic-0.12.0.sif pyscenic ctx /data/scenic_grn.csv /data/flank_ctx.genes_vs_tracks.rankings.feather \
  --annotations_fname /data/ctx_annotation.tsv --expression_mtx_fname /data/rna.loom --output /data/scenic_flank_reg.csv --rank_threshold 1500 \
  --min_genes 6 --num_workers 20 --cell_id_attribute cells --gene_attribute features', shell = True, executable="/bin/bash")

# %% tags=[]
flank_merged = pd.read_csv(f"{PATH}/scenic_flank_reg.csv", header=None, skiprows=3, usecols=[0, 8], names=["tf", "targets"])
flank_merged["targets"] = flank_merged["targets"].map(lambda x: set(i[0] for i in eval(x)))
# The following sentence merge repeated tfs and their corresponding targets
flank_merged = flank_merged.groupby("tf").aggregate({"targets": lambda x: functools.reduce(set.union, x)})
flank_merged["n_targets"] = flank_merged["targets"].map(len)
flank_merged = flank_merged.sort_values("n_targets", ascending=False)
flank_merged

# %%
g = nx.DiGraph()
for tf, row in flank_merged.iterrows():
    for target in row["targets"]:
        g.add_edge(tf, target)
nx.set_node_attributes(g, "target", name="type")
for tf in flank_merged.index:
    g.nodes[tf]["type"] = "TF"
nx.write_graphml(g, f"{PATH}/flank_merged.graphml.gz")

# %%
nx.to_pandas_edgelist(
    g, source="TF", target="Target gene"
).to_csv(f"{PATH}/flank_merged.csv", index=False)

plt.figure(figsize=(25, 25))
nx.draw(g, graphviz_layout(g), with_labels=True)
plt.savefig(f"{PATH}/flank_merged_graph.png", dpi=300)
plt.close()

# %% [markdown]
# ## Distance

# %%
dist_feather = rank_gene_tf_dist.T
dist_feather = dist_feather.loc[np.unique(dist_feather.index), np.unique(dist_feather.columns)].astype(np.int16)
dist_feather.index += "_atac"
dist_feather.index.name = "tracks"
dist_feather.columns.name = None
columns = dist_feather.columns.tolist()
dist_feather = dist_feather.reset_index()
dist_feather = dist_feather.loc[:, [*columns, "tracks"]]
dist_feather.to_feather(f"{PATH}/dist_ctx.genes_vs_tracks.rankings.feather")

subprocess.run(f'singularity run --bind {PATH}:/data /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/scGLUE/aertslab-pyscenic-0.12.0.sif pyscenic ctx /data/scenic_grn.csv /data/dist_ctx.genes_vs_tracks.rankings.feather /data/flank_ctx.genes_vs_tracks.rankings.feather \
  --annotations_fname /data/ctx_annotation.tsv --expression_mtx_fname /data/rna.loom --output /data/scenic_dist_reg.csv --rank_threshold 1500 \
  --min_genes 6 --num_workers 20 --cell_id_attribute cells --gene_attribute features', shell = True, executable="/bin/bash")

# %% tags=[]
dist_merged = pd.read_csv(f"{PATH}/scenic_dist_reg.csv", header=None, skiprows=3, usecols=[0, 8], names=["tf", "targets"])
dist_merged["targets"] = dist_merged["targets"].map(lambda x: set(i[0] for i in eval(x)))
dist_merged = dist_merged.groupby("tf").aggregate({"targets": lambda x: functools.reduce(set.union, x)})
dist_merged["n_targets"] = dist_merged["targets"].map(len)
dist_merged = dist_merged.sort_values("n_targets", ascending=False)
dist_merged

# %%
g = nx.DiGraph()
for tf, row in dist_merged.iterrows():
    for target in row["targets"]:
        g.add_edge(tf, target)
nx.set_node_attributes(g, "target", name="type")
for tf in dist_merged.index:
    g.nodes[tf]["type"] = "TF"
nx.write_graphml(g, f"{PATH}/dist_merged.graphml.gz")

# %%
nx.to_pandas_edgelist(
    g, source="TF", target="Target gene"
).to_csv(f"{PATH}/dist_merged.csv", index=False)

plt.figure(figsize=(25, 25))
nx.draw(g, graphviz_layout(g), with_labels=True)
plt.savefig(f"{PATH}/dist_merged_graph.png", dpi=300)
plt.close()

# %% [markdown]
# ## eQTL

# %%
eqtl_feather = rank_gene_tf_eqtl.T
eqtl_feather = eqtl_feather.loc[np.unique(eqtl_feather.index), np.unique(eqtl_feather.columns)].astype(np.int16)
eqtl_feather.index += "_atac"
eqtl_feather.index.name = "tracks"
eqtl_feather.columns.name = None
columns = eqtl_feather.columns.tolist()
eqtl_feather = eqtl_feather.reset_index()
eqtl_feather = eqtl_feather.loc[:, [*columns, "tracks"]]
eqtl_feather.to_feather(f"{PATH}/eqtl_ctx.genes_vs_tracks.rankings.feather")

subprocess.run(f'singularity run --bind {PATH}:/data /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/scGLUE/aertslab-pyscenic-0.12.0.sif pyscenic ctx /data/scenic_grn.csv /data/eqtl_ctx.genes_vs_tracks.rankings.feather /data/flank_ctx.genes_vs_tracks.rankings.feather \
  --annotations_fname /data/ctx_annotation.tsv --expression_mtx_fname /data/rna.loom --output /data/scenic_eqtl_reg.csv --rank_threshold 1500 \
  --min_genes 6 --num_workers 20 --cell_id_attribute cells --gene_attribute features', shell = True, executable="/bin/bash")

# %% tags=[]
eqtl_merged = pd.read_csv(f"{PATH}/scenic_eqtl_reg.csv", header=None, skiprows=3, usecols=[0, 8], names=["tf", "targets"])
eqtl_merged["targets"] = eqtl_merged["targets"].map(lambda x: set(i[0] for i in eval(x)))
eqtl_merged = eqtl_merged.groupby("tf").aggregate({"targets": lambda x: functools.reduce(set.union, x)})
eqtl_merged["n_targets"] = eqtl_merged["targets"].map(len)
eqtl_merged = eqtl_merged.sort_values("n_targets", ascending=False)
eqtl_merged

# %%
g = nx.DiGraph()
for tf, row in eqtl_merged.iterrows():
    for target in row["targets"]:
        g.add_edge(tf, target)
nx.set_node_attributes(g, "target", name="type")
for tf in eqtl_merged.index:
    g.nodes[tf]["type"] = "TF"
nx.write_graphml(g, f"{PATH}/eqtl_merged.graphml.gz")

# %%
nx.to_pandas_edgelist(
    g, source="TF", target="Target gene"
).to_csv(f"{PATH}/eqtl_merged.csv", index=False)

plt.figure(figsize=(25, 25))
nx.draw(g, graphviz_layout(g), with_labels=True)
plt.savefig(f"{PATH}/eqtl_merged_graph.png", dpi=300)
plt.close()

# %% [markdown]
# ## GLUE

# %%
glue_feather = rank_gene_tf_glue.T
glue_feather = glue_feather.loc[np.unique(glue_feather.index), np.unique(glue_feather.columns)].astype(np.int16)
glue_feather.index += "_atac"
glue_feather.index.name = "tracks"
glue_feather.columns.name = None
columns = glue_feather.columns.tolist()
glue_feather = glue_feather.reset_index()
glue_feather = glue_feather.loc[:, [*columns, "tracks"]]
glue_feather.to_feather(f"{PATH}/glue_ctx.genes_vs_tracks.rankings.feather")

subprocess.run(f'singularity run --bind {PATH}:/data /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/scGLUE/aertslab-pyscenic-0.12.0.sif pyscenic ctx /data/scenic_grn.csv /data/glue_ctx.genes_vs_tracks.rankings.feather /data/flank_ctx.genes_vs_tracks.rankings.feather \
  --annotations_fname /data/ctx_annotation.tsv --expression_mtx_fname /data/rna.loom --output /data/scenic_glue_reg.csv --rank_threshold 1500 \
  --min_genes 6 --num_workers 20 --cell_id_attribute cells --gene_attribute features', shell = True, executable="/bin/bash")

# %% tags=[]
glue_merged = pd.read_csv(f"{PATH}/scenic_glue_reg.csv", header=None, skiprows=3, usecols=[0, 8], names=["tf", "targets"])
glue_merged["targets"] = glue_merged["targets"].map(lambda x: set(i[0] for i in eval(x)))
glue_merged = glue_merged.groupby("tf").aggregate({"targets": lambda x: functools.reduce(set.union, x)})
glue_merged["n_targets"] = glue_merged["targets"].map(len)
glue_merged = glue_merged.sort_values("n_targets", ascending=False)
glue_merged

# %%
g = nx.DiGraph()
for tf, row in glue_merged.iterrows():
    for target in row["targets"]:
        g.add_edge(tf, target)
nx.set_node_attributes(g, "target", name="type")
for tf in glue_merged.index:
    g.nodes[tf]["type"] = "TF"
nx.write_graphml(g, f"{PATH}/glue_merged.graphml.gz")

# %%
nx.to_pandas_edgelist(
    g, source="TF", target="Target gene"
).to_csv(f"{PATH}/glue_merged.csv", index=False)

plt.figure(figsize=(25, 25))
nx.draw(g, graphviz_layout(g), with_labels=True)
plt.savefig(f"{PATH}/glue_merged_graph.png", dpi=300)
plt.close()




