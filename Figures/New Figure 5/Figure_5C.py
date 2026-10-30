import pandas as pd
import os
import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout

os.chdir("../")
caQTL_glue_csv = pd.read_csv("SMR/save/s040_infer_gene_tf_SMR_pSMR_5e-8/glue_merged.csv")
raw_glue_csv = pd.read_csv("save/s04_infer_gene_tf/glue_merged.csv")

caQTL_glue_set = set(map(tuple, caQTL_glue_csv.to_records(index=False)))
raw_glue_set = set(map(tuple, raw_glue_csv.to_records(index=False)))

# diff_pairs returns all the TF-target gene pairs that only exist in the caQTL set, not the raw set
diff_pairs = [pair for pair in caQTL_glue_set if pair not in raw_glue_set]

# diff_pair_genes = list({gene for pair in diff_pairs for gene in pair})
diff_pair_target_genes = [pair[1] for pair in diff_pairs]

def aggregate_SMR(dir, trait_name):
    for chr in range(1, 23):
        if os.path.exists(dir + "Chr" + str(chr) + "_results.smr") == False:
            continue
        SMR_chr = pd.read_csv(dir + "Chr" + str(chr) + "_results.smr", sep = "\t")
        SMR_chr['trait'] = trait_name
        if chr == 1:
            SMR_agg = SMR_chr
        else:
            SMR_agg = pd.concat([SMR_agg, SMR_chr])
    return SMR_agg

# %% [SMR outputs]
eQTL_GWAS_output = "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_eQTL/CD14_Mono/disease_traits/"
smr_eQTL_GWAS = pd.concat([aggregate_SMR(eQTL_GWAS_output + trait + "/", trait) for trait in sorted(os.listdir(eQTL_GWAS_output))])
SMR_genes = smr_eQTL_GWAS['Gene']

common_genes = set(diff_pair_target_genes).intersection(set(SMR_genes))

# Returns the TF-target gene pairs, where the target gene is in the peak-to-gene lists given by SMR
diff_pairs_valid = [pair for pair in diff_pairs if pair[1] in common_genes]
diff_pairs_valid = [pair for pair in diff_pairs_valid if pair[0] != pair[1]]

g = nx.DiGraph(diff_pairs_valid)

plt.figure(figsize=(15, 15))
nx.draw(g, graphviz_layout(g), with_labels=True)
plt.savefig("SMR/save/s040_infer_gene_tf_SMR_pSMR_5e-8/glue_caQTL_specific.png", dpi=300)
plt.close()

# Only draw the subgraph where the gene of interest exists
target_gene = "APOBEC3A"
component = next(c for c in nx.weakly_connected_components(g) if target_gene in c)
g_sub = g.subgraph(component)

plt.figure(figsize=(10, 10))
pos = graphviz_layout(g)
sub_pos = {n: pos[n] for n in g_sub.nodes}
nx.draw(g_sub, pos=sub_pos, with_labels=True, font_size=25, arrowsize=25)
plt.savefig("SMR/save/s040_infer_gene_tf_SMR_pSMR_5e-8/glue_caQTL_specific_" + target_gene + ".png", dpi=300)
plt.close()

