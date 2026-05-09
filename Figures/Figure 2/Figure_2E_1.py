import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis

# Set parameters
peak = "chr13:24670806-24672096"
peak_name = peak.replace(":", "_").replace("-", "_")
chr_num = peak.split(":")[0].replace("chr", "")

cd14_mono_snp = "13:24570579:C:A"
cd4_tcm_snp = "13:24671328:T:C"

print(f"Processing peak: {peak_name}")
print(f"Chromosome: {chr_num}")
print(f"CD14_Mono SNP: {cd14_mono_snp}")
print(f"CD4_TCM SNP: {cd4_tcm_snp}")

# Genotype data
geno_path = "/g/data/ei56/ax3061/proj/tenk10k/caQTL/data/genotype"
plink_prefix_path = f"{geno_path}/TenK10K_TOB_ATAC_renamed_chr{chr_num}_common_variants_qced"

pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()

# Extract only the two SNPs needed for this example
genotype_df_subset = genotype_df.loc[[cd14_mono_snp, cd4_tcm_snp], :]

genotype_df_subset.to_csv(
    f"/g/data/ei56/od8037/Plotting/QTL_Boxplots/Example_2/genotype_{peak_name}.csv"
)

# Phenotype data
for ct_name in ["CD14_Mono", "CD4_TCM"]:
    print(f"Processing {ct_name} for peak {peak}")
    
    expression_bed = f"/g/data/ei56/od8037/Final_caQTL/Runs/{ct_name}/ExpressionBeds/chr{chr_num}.bed.gz"
    covariates_file = f"/g/data/ei56/od8037/Final_caQTL/Runs/{ct_name}/covariates.txt"

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    phenotype_df.columns = phenotype_df.columns.str.replace("X", "", regex=False)

    covariates_df = pd.read_csv(covariates_file, index_col=0).T

    # Remove samples with NA covariates
    na_pos = covariates_df[covariates_df.isna().any(axis=1)].index.tolist()
    phenotype_df = phenotype_df.drop(columns=na_pos)

    # Extract only the chosen peak
    phenotype_df_subset = phenotype_df.loc[[peak], :]

    phenotype_df_subset.to_csv(
        f"/g/data/ei56/od8037/Plotting/QTL_Boxplots/Example_2/{ct_name}_{peak_name}.csv"
    )
