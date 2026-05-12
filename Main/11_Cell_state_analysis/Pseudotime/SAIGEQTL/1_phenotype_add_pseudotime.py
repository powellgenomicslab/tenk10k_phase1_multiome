
import pandas as pd
import os
import scanpy as sc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--cell_type', required=True)
parser.add_argument('--chr', required=True)
args = parser.parse_args()

cell_type = args.cell_type
chr = args.chr

palantir_sc = sc.read_h5ad("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CD4_combined.h5ad")

palantir_df = pd.read_csv("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Phenotypes/Palantir_output_all.csv")
palantir_df['barcode'] = palantir_df['barcode'].str.replace(
    r'^(.+?)_(.+)$', r'\2_\1', regex=True
)

finaltsv_dir = "/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/FinalTSVs/"
save_dir = "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SAIGEQTL/Dynamic/Phenotypes/FinalTSVs/"
cell_type_list = sorted(os.listdir(finaltsv_dir))

os.makedirs(os.path.join(save_dir, cell_type), exist_ok = True)

phenotype_df = pd.read_csv(os.path.join(finaltsv_dir, cell_type, f"chr{chr}.tsv"), sep = '\t')
phenotype_df = phenotype_df.merge(palantir_df, on='barcode', how='inner')

cols = phenotype_df.columns.tolist()
pc5_idx = cols.index('PC5')
cols.insert(pc5_idx + 1, cols.pop(cols.index('pseudotime')))
phenotype_df = phenotype_df[cols]

phenotype_df.to_csv(os.path.join(save_dir, cell_type, f"chr{chr}.tsv"), sep='\t', index=False)