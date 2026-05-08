
import anndata
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import scglue
import os

PATH = "save/s04_infer_gene_tf"

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

# %%
feature_embedding = pd.read_csv(f"save/s02_glue/ckpt/prior:dcq/feature_embeddings.csv", header=None, index_col=0)

## %%
glue_graph = scglue.genomics.regulatory_inference(
   feature_embedding.index, feature_embedding.to_numpy(), dist_graph.subgraph([*genes.index, *peaks.index]),
   alternative="greater", random_state=0
)

glue_edges = nx.to_pandas_edgelist(glue_graph)

# %%
# Formatting the significant SMR used as guidance to the eQTL-like format.

def aggregate_SMR(dir):
    for chr in range(1, 23):
        if os.path.exists(dir + "Chr" + str(chr) + "_results.smr") == False:
            continue
        SMR_chr = pd.read_csv(dir + "Chr" + str(chr) + "_results.smr", sep = "\t")
        if chr == 1:
            SMR_agg = SMR_chr
        else:
            SMR_agg = pd.concat([SMR_agg, SMR_chr])
    return SMR_agg
    
smr_peak_gene = aggregate_SMR("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/SMR_run_rawID/CD14_Mono/")
smr_sig_peak_gene = smr_peak_gene[smr_peak_gene['p_SMR'] < 5e-8]

valid_gene_list = feature_embedding.index[:6000].tolist()
valid_peak_list = feature_embedding.index[6000:].tolist()

smr_sig_peak_gene_valid = smr_sig_peak_gene[(smr_sig_peak_gene['Outco_Gene'].isin(valid_gene_list)) & (smr_sig_peak_gene['Expo_Gene'].isin(valid_peak_list))]
smr_sig_peak_gene_valid = smr_sig_peak_gene_valid[['Outco_Gene', 'Expo_Gene']].rename(columns={'Outco_Gene': 'source', 'Expo_Gene': 'target'})

# %%
# Refine the GLUE score, only focus on the SMR highlighted gene-peak pairs
glue_edges_APOBEC3A = glue_edges[glue_edges['source'] == 'APOBEC3A']
keep_attrs=["score"]
gene2peak_link_glue = glue_edges_APOBEC3A.merge(
       scglue.genomics.Bed(rna.var).strand_specific_start_site().df.iloc[:, :4], how="left", left_on="source", right_on="name"
   ).merge(
       scglue.genomics.Bed(atac.var).df.iloc[:, :4], how="left", left_on="target", right_on="name"
   ).loc[:, [
       "chrom_x", "chromStart_x", "chromEnd_x",
       "chrom_y", "chromStart_y", "chromEnd_y",
       *(keep_attrs or [])
   ]].dropna()
gene2peak_link_glue["chromStart_x"] = gene2peak_link_glue["chromStart_x"].astype('int')
gene2peak_link_glue["chromEnd_x"] = gene2peak_link_glue["chromEnd_x"].astype('int')
gene2peak_link_glue["chromStart_y"] = gene2peak_link_glue["chromStart_y"].astype('int')
gene2peak_link_glue["chromEnd_y"] = gene2peak_link_glue["chromEnd_y"].astype('int')
gene2peak_link_glue.to_csv(f"{PATH}/gene2peak_APOBEC3A.links", sep="\t", index=False, header=False)

# Modify tracks_SMR.ini file here, change the glue score links file.
subprocess.run('pyGenomeTracks --tracks tracks_SMR.ini --region chr22:38940000-39000000 --outFileName tracks_SMR_APOBEC3A.png', shell = True, executable="/bin/bash")
