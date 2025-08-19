import functools
import itertools
import os

import anndata
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix

import scglue
from muon import atac as ac

# %%
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

os.makedirs('save/', exist_ok = True)
os.makedirs('save/s01_preprocessing/', exist_ok = True)

# %%
rna = anndata.read_h5ad("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/data_processing/cell_type_spec/RNA_classification/output_celltype/B_memory/RNA_B_memory_ALL.h5ad")
# %%
atac = anndata.read_h5ad("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/data_processing/cell_type_spec_TenK10K_phase1/ATAC_classification/output_celltype_238lib/B_memory/ATAC_B_memory.h5ad")
# %%

# Annotate rna with chrom, chromStart, and chromEnd
scglue.data.get_gene_annotation(
    # rna, gtf="/directflow/SCCGGroupShare/projects/jayfan/Projects/scGLUE/DEMO/gencode/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz",
    rna, gtf="/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/scGLUE/metadata/gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz",
    gtf_by="gene_name"
)

# Some of the gene expression reads included in TOB measurements cannot find the corresponding gene name in GENCODE.
# We therefore need to at first filter it.

nan_genes = pd.isna(rna.var['chrom'])
exist_gene_list = [name for name in rna.var_names if not nan_genes[name]]
rna = rna[:, exist_gene_list]

# Add chrom indexes for ATAC-seq data.

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)

# %%
sc.pp.filter_genes(rna, min_counts=1)
rna.obs_names += "-RNA"

# %%
sc.pp.filter_genes(atac, min_counts=1)
atac.obs_names += "-ATAC"

# %%
genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names))
tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)

# %% [markdown]
# # RNA

# %%
rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=6000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, n_comps=100, use_highly_variable=True, svd_solver="auto")
# Using harmony for batch correction
sce.pp.harmony_integrate(adata = rna, key = 'sequencing_library', basis = 'X_pca')
sc.pp.neighbors(rna, n_pcs=100, use_rep="X_pca_harmony", metric="cosine")
sc.tl.umap(rna)

# %%
rna.X = rna.layers["counts"]
# del rna.layers["counts"]

# %%
os.makedirs('save/s01_preprocessing/figure/', exist_ok = True)
# fig = sc.pl.umap(rna, color="wg2_scpred_prediction", title="scRNA-seq cell type", return_fig=True)
# fig.savefig("save/s01_preprocessing/figure/rna_ct.pdf")
# fig = sc.pl.umap(rna, color="sequencing_library", title="scRNA-seq sequencing library", return_fig=True)
# fig.savefig("save/s01_preprocessing/figure/rna_lib.pdf")

merge_mapping = {}
rna_libraries = np.unique(rna.obs['sequencing_library'])
for rna_library in rna_libraries:
    merge_mapping[rna_library] = rna_library[:-1]
rna.obs['merged_library'] = rna.obs['sequencing_library'].map(merge_mapping)
fig = sc.pl.umap(rna, color="merged_library", title="scRNA-seq merged library", return_fig=True)
fig.savefig("save/s01_preprocessing/figure/rna_merge_lib.pdf")

# %% [markdown]
# # ATAC

atac.layers["counts"] = atac.X
# TF-IDF
# ac.pp.tfidf(atac, scale_factor=1e4)
# sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=.5)
# np.sum(atac.var.highly_variable)
# atac = atac[:, atac.var["highly_variable"]].copy()

# %%
scglue.data.lsi(atac, n_components=100, use_highly_variable=False, n_iter=15)
atac.obsm['X_lsi'] = atac.obsm['X_lsi'][:, 1:]
# Using harmony for batch correction
sce.pp.harmony_integrate(adata = atac, key = 'library', basis = 'X_lsi')
sc.pp.neighbors(atac, n_pcs=99, use_rep="X_pca_harmony", metric="cosine")
sc.tl.umap(atac)

# %%
fig = sc.pl.umap(atac, color="predicted.id", title="scATAC-seq cell type", return_fig=True)
fig.savefig("save/s01_preprocessing/figure/atac_ct.pdf")
fig = sc.pl.umap(atac, color="library", title="scATAC-seq sequencing library", return_fig=True)
fig.savefig("save/s01_preprocessing/figure/atac_lib.pdf")

merge_mapping = {}
atac_libraries = np.unique(atac.obs['library'])
for atac_library in atac_libraries:
    merge_mapping[atac_library] = atac_library[:-2]
atac.obs['merged_library'] = atac.obs['library'].map(merge_mapping)
fig = sc.pl.umap(atac, color="merged_library", title="scATAC-seq merged library", return_fig=True)
fig.savefig("save/s01_preprocessing/figure/atac_merge_lib.pdf")

# %% [markdown]
# # Build graph

# %% [markdown]
# ## Overlap

# %%
overlap_graph = scglue.genomics.window_graph(
    genes.expand(2000, 0), peaks, 0,
    attr_fn=lambda l, r, d: {
        "weight": 1.0,
        "type": "overlap"
    }
)
overlap_graph = nx.DiGraph(overlap_graph)
overlap_graph.number_of_edges()

# %% [markdown]
# ## Genomic distance

# %%
dist_graph = scglue.genomics.window_graph(
    promoters, peaks, 150000,
    attr_fn=lambda l, r, d: {
        "dist": abs(d),
        "weight": scglue.genomics.dist_power_decay(abs(d)),
        "type": "dist"
    }
)
dist_graph = nx.DiGraph(dist_graph)
dist_graph.number_of_edges()

# %% [markdown]
# ## eQTL

# %%
rna.var["gene_id_trimmed"] = rna.var["gene_id"].map(scglue.genomics.ens_trim_version)
gene_id_mapping = {ens: name for ens, name in zip(rna.var["gene_id_trimmed"], rna.var_names)}

# %%
eqtl = pd.read_csv(eqtl_path, sep='\t', comment='#')
eqtl = eqtl.loc[eqtl['chromStart'].str.isdigit()]
eqtl = scglue.genomics.Bed(eqtl)

# %%
eqtl_graph = scglue.genomics.window_graph(
    eqtl, peaks, 0, left_sorted=False,
    attr_fn=lambda l, r, d: {
        "weight": 1.0,
        "type": "eqtl"
    }
)
eqtl_graph = nx.DiGraph(eqtl_graph)
eqtl_graph.number_of_edges()

# %%
eqtl_genes = eqtl['name']
eqtl_genes = set(eqtl_genes)
rna.var["in_eqtl"] = [item in eqtl_genes for item in rna.var_names]
rna.var["in_eqtl"].sum()

# %%
eqtl_links = eqtl.df.iloc[:, :4].merge(tss.df.iloc[:, :4], how="left", on="name").assign(score=1)
eqtl_links = eqtl_links.query("chrom_x == chrom_y")
eqtl_links["name"] = eqtl_links.pop("name")
# eqtl_links.to_csv(f"{PATH}/eqtl.annotated_links", sep="\t", index=False, header=False)

# %% [markdown]
# # Update highly variable genes

# %%
rna.var["o_highly_variable"] = rna.var["highly_variable"]
print("RNA: Number of o_highly_variable -- ", rna.var["o_highly_variable"].sum())

# %%
# The following block aims to filter out genes/promoters with at least one peaks associated.
# 'A1' actually flattens the matrix into a vector.
rna.var["in_cicero"] = biadjacency_matrix(
    scglue.genomics.window_graph(promoters, peaks, 0),
    genes.index
).sum(axis=1).A1 > 0
print("RNA: Counts within cicero -- ", rna.var["in_cicero"].sum())

# %%
rna.var["d_highly_variable"] = functools.reduce(np.logical_and, [
    rna.var["highly_variable"],
    rna.var["in_eqtl"],
    rna.var["in_cicero"]
])
print("RNA: Number of d_highly_variable -- ", rna.var["d_highly_variable"].sum())

# %%
rna.var["dcq_highly_variable"] = rna.var["highly_variable"]
print("RNA: Number of dcq_highly_variable -- ", rna.var["dcq_highly_variable"].sum())


# %% [markdown]
# # Combine graphs into priors

# %% [markdown]
# ## Overlap

# %%
o_prior = overlap_graph.copy()

# %%
hvg_reachable = scglue.graph.reachable_vertices(o_prior, rna.var.query("o_highly_variable").index)

# %%
atac.var["o_highly_variable"] = [item in hvg_reachable for item in atac.var_names]
print("ATAC: Number of o_highly_variable -- ", atac.var["o_highly_variable"].sum())

# %%
o_prior = scglue.graph.compose_multigraph(o_prior, o_prior.reverse())
for item in itertools.chain(atac.var_names, rna.var_names):
    o_prior.add_edge(item, item, weight=1.0, type="self-loop")
nx.set_edge_attributes(o_prior, 1, "sign")

# %%
o_prior = o_prior.subgraph(hvg_reachable)

# %% [markdown]
# ## Genomic distance

# %%
d_prior = dist_graph.copy()

# %%
hvg_reachable = scglue.graph.reachable_vertices(d_prior, rna.var.query("d_highly_variable").index)

# %%
atac.var["d_highly_variable"] = [item in hvg_reachable for item in atac.var_names]
print("ATAC: Number of d_highly_variable -- ", atac.var["d_highly_variable"].sum())

# %%
d_prior = scglue.graph.compose_multigraph(d_prior, d_prior.reverse())
for item in itertools.chain(atac.var_names, rna.var_names):
    d_prior.add_edge(item, item, weight=1.0, type="self-loop")
nx.set_edge_attributes(d_prior, 1, "sign")

# %%
d_prior = d_prior.subgraph(hvg_reachable)

# %% [markdown]
# ## Genomic distance + pcHi-C + eQTL

# %%
dcq_prior = scglue.graph.compose_multigraph(dist_graph, eqtl_graph)

# %%
hvg_reachable = scglue.graph.reachable_vertices(dcq_prior, rna.var.query("dcq_highly_variable").index)

# %%
atac.var["dcq_highly_variable"] = [item in hvg_reachable for item in atac.var_names]
print("ATAC: Number of dcq_highly_variable -- ", atac.var["dcq_highly_variable"].sum())

# %%
dcq_prior = scglue.graph.compose_multigraph(dcq_prior, dcq_prior.reverse())
for item in itertools.chain(atac.var_names, rna.var_names):
    dcq_prior.add_edge(item, item, weight=1.0, type="self-loop")
nx.set_edge_attributes(dcq_prior, 1, "sign")

# %%
dcq_prior = dcq_prior.subgraph(hvg_reachable)


# %% [markdown]
# # Write data

# %%
save_dir = 'save/s01_preprocessing/'
omics_dir = save_dir + 'omics_data/'
network_dir = save_dir + 'network/'
os.makedirs(omics_dir, exist_ok = True)
os.makedirs(network_dir, exist_ok = True)

rna.__dict__['_raw'].__dict__['_var'] = rna.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
# atac.__dict__['_raw'].__dict__['_var'] = atac.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
to_del_list = ['highly_variable_rank', 'hgnc_id', 'tag', 'havana_gene', 'artif_dupl']
for to_del_var in to_del_list:
    if (to_del_var in rna.var.keys()):
        del rna.var[to_del_var]

rna.write(omics_dir + "rna.h5ad", compression="gzip")
atac.write(omics_dir + "atac.h5ad", compression="gzip")

# %%
nx.write_graphml(overlap_graph, network_dir + "overlap.graphml.gz")
nx.write_graphml(dist_graph, network_dir + "dist.graphml.gz")
nx.write_graphml(eqtl_graph, network_dir + "eqtl.graphml.gz")

# %%
nx.write_graphml(o_prior, network_dir + "o_prior.graphml.gz")
nx.write_graphml(d_prior, network_dir + "d_prior.graphml.gz")
nx.write_graphml(dcq_prior, network_dir + "dcq_prior.graphml.gz")




