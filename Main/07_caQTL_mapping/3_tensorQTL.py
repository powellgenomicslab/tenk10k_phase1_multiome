# Import libraries
import numpy as np
import rpy2
from rpy2.robjects.packages import importr
qvalue = importr('qvalue')
print(np.array(qvalue.qvalue(rpy2.robjects.vectors.FloatVector([0.001,0.01,0.1,0.4,0.4,0.5]), **{'lambda':rpy2.robjects.vectors.FloatVector([0.5])}).rx2('qvalues')))

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")
import os
import sys
from pathlib import Path

# Convert the argument to integer
inputArgs = sys.argv[1:]
chr_num = int(inputArgs[0])
ct_name = inputArgs[1]
print(f"Chromosome: {chr_num}, Cell type: {ct_name}")

# Set directory to save results
result_dir = f"/g/data/ei56/od8037/NewGenotypes/caQTL/Runs/{ct_name}/Results"
Path(result_dir).mkdir(parents=True, exist_ok=True)
os.chdir(result_dir)

# Define paths to data
geno_path = "/g/data/ei56/ax3061/proj/tenk10k/caQTL/data/genotype"
plink_prefix_path = f"{geno_path}/TenK10K_TOB_ATAC_renamed_chr{chr_num}_common_variants_qced"
expression_bed = f"/g/data/ei56/od8037/Final_caQTL/Runs/{ct_name}/ExpressionBeds/chr{chr_num}.bed.gz"
covariates_file = f"/g/data/ei56/od8037/Final_caQTL/Runs/{ct_name}/covariates.txt" 
prefix = "TenK10K"

# Load phenotypes and covariates
print("Load phenotypes and covariates")
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
phenotype_df.columns = phenotype_df.columns.str.replace("X", "")

# Find NAs in covariates matrix
covariates_df = pd.read_csv(covariates_file, index_col=0).T
na_pos = covariates_df[covariates_df.isna().any(axis=1)].index.tolist()
# Remove corresponding columns from phenotype matrix
phenotype_df = phenotype_df.drop(columns = na_pos)

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
# Chromosome name changed format
variant_df.chrom = "chr" + variant_df.chrom


## cis-QTL: nominal p-values for all variant-phenotype pairs
# Map all cis-associations (results for each chromosome are written to file)

print("Match all the input files")
# Some donors are missing from genotype information, so do not analyse them
shared_donors = list(set(phenotype_df.columns) & set(genotype_df.columns))
covariates_df = covariates_df.loc[shared_donors, :]
genotype_df = genotype_df.loc[:, shared_donors]
phenotype_df = phenotype_df.loc[:, shared_donors]

# Genes on a specific chromosome
cis.map_nominal(
    genotype_df, 
    variant_df,
    phenotype_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
    phenotype_pos_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
    prefix, 
    covariates_df = covariates_df,
    write_top = True, 
    write_stats = True,
    maf_threshold = 0.05,
    window = 1000000
)

# Load results
print("Loading results")
pairs_df = pd.read_parquet(f"{prefix}.cis_qtl_pairs.chr{chr_num}.parquet")
print(pairs_df.head())

print("Saving all association pairs")
pairs_df.to_csv(f"{prefix}.cis_qtl_pairs.chr{chr_num}.csv", sep = "\t")


## cis-QTL: empirical p-values for phenotypes

print("Start permutation test")
# Genes on a specific chromosome
# This only output the top signal per gene
cis_df = cis.map_cis(
    genotype_df, 
    variant_df, 
    phenotype_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
    phenotype_pos_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
    covariates_df = covariates_df, 
    seed = 123456,
    maf_threshold = 0.05, 
    nperm = 10000,
    window = 1000000
)

# Calculate chromosome-wide FDR
# Need to understand why the author sets lamdba = 0.85
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda = 0.85)
print(cis_df.head())

# Save the results
print("Saving the significant association pairs")
cis_df.to_csv(f"{prefix}.sig_cis_qtl_pairs.chr{chr_num}.csv", sep = "\t")


# Conditionally independent QTLs
if any(cis_df.qval < 0.05):
	print("Identify conditionally independent caQTLs")
	indep_df = cis.map_independent(
		genotype_df, 
		variant_df, 
		cis_df,
        phenotype_df, 
		phenotype_pos_df, 
        covariates_df = covariates_df,
        maf_threshold = 0.05, 
		nperm = 10000,
        fdr = 0.05, 
		fdr_col = "qval",
        window = 1000000
	)

	print("Saving the conditionally independent caQTLs")
	indep_df.to_csv(f"{prefix}.independent_cis_qtl_pairs.chr{chr_num}.csv", sep = "\t")

else :
	print("No significant caQTLs with qval < 0.05. Do not identify conditionally independent QTLs.")

print("Analysis finished!")