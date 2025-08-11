import os
import re

import anndata
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams

import scglue

# %%
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

# %% [markdown]
# # Read data

# %%
rna = anndata.read_h5ad("save/s01_preprocessing/omics_data/rna.h5ad")
atac = anndata.read_h5ad("save/s01_preprocessing/omics_data/atac.h5ad")

PRIOR = "dcq"
prior = nx.read_graphml(f"save/s01_preprocessing/network/{PRIOR}_prior.graphml.gz")

# %%
rna.var["highly_variable"] = rna.var[f"{PRIOR}_highly_variable"]
atac.var["highly_variable"] = atac.var[f"{PRIOR}_highly_variable"]
# rna.var["highly_variable"].sum(), atac.var["highly_variable"].sum()

# %% [markdown]

os.makedirs("save/s02_glue/", exist_ok=True)
PATH = "save/s02_glue/ckpt"
os.makedirs(PATH, exist_ok=True)

# # Train model

# %%
scglue.models.configure_dataset(rna, "NB", use_highly_variable=True, use_layer="counts", use_rep="X_pca_harmony")
scglue.models.configure_dataset(atac, "NB", use_highly_variable=True, use_layer="counts", use_rep="X_pca_harmony")

# %% tags=[]
glue = scglue.models.SCGLUEModel(
    {"rna": rna, "atac": atac}, sorted(prior.nodes)
)
glue.compile()
glue.fit(
    {"rna": rna, "atac": atac}, prior,
    align_burnin=np.inf, safe_burnin=False,
    directory=f"{PATH}/prior:{PRIOR}/pretrain"
)
glue.save(f"{PATH}/prior:{PRIOR}/pretrain/final.dill")

# %%
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
scglue.data.estimate_balancing_weight(
    rna, atac, use_rep="X_glue"
)

# %%
scglue.models.configure_dataset(rna, "NB", use_highly_variable=True, use_layer="counts", use_rep="X_pca_harmony", use_dsc_weight="balancing_weight")
scglue.models.configure_dataset(atac, "NB", use_highly_variable=True, use_layer="counts", use_rep="X_pca_harmony", use_dsc_weight="balancing_weight")

# %%
glue = scglue.models.SCGLUEModel(
    {"rna": rna, "atac": atac}, sorted(prior.nodes)
)
glue.adopt_pretrained_model(scglue.models.load_model(
    f"{PATH}/prior:{PRIOR}/pretrain/final.dill"
))
glue.compile()
glue.fit(
    {"rna": rna, "atac": atac}, prior,
    directory=f"{PATH}/prior:{PRIOR}/fine-tune"
)
# glue.save(f"{PATH}/prior:{PRIOR}/fine-tune/final.dill")

# Check whether the integration is reliable, which quantifies the consistency between the integration result and the guidance graph
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, prior
)
print("---Below we print dx, the integration consistency score---")
print(dx)

# %% [markdown]
# # Embeddings

# %% [markdown]
# ## Cell embeddings

# %%
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

# %%
combined = anndata.AnnData(
    obs=pd.concat([rna.obs, atac.obs], join="inner"),
    obsm={"X_glue": np.concatenate([rna.obsm["X_glue"], atac.obsm["X_glue"]])}
)

# %%
sc.pp.neighbors(combined, n_pcs=50, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)

# %%
rna.write(f"{PATH}/prior:{PRIOR}/rna_glue.h5ad", compression="gzip")
atac.write(f"{PATH}/prior:{PRIOR}/atac_glue.h5ad", compression="gzip")
combined.write(f"{PATH}/prior:{PRIOR}/combined_glue.h5ad", compression="gzip")

# %% [markdown]
# ## Feature embeddings

# %%
feature_embeddings = pd.DataFrame(
    glue.encode_graph(prior),
    index=glue.vertices
)
feature_embeddings.iloc[:5, :5]

# %%
feature_embeddings.to_csv(f"{PATH}/prior:{PRIOR}/feature_embeddings.csv", index=True, header=False)



