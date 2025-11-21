# Load libraries
library(data.table)
library(glue)
library(fst)

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
condition <- args[1]

gwas_result_dir <- "/g/data/fy54/data/external_data/gwas/gymrek-ukbb-snp-gwas-catalogs_v6"

print(glue("Loading GWAS data for {condition} ..."))
gwas_df <- fread(glue(
    "{gwas_result_dir}/white_british_{condition}_snp_gwas_results_hg38.tab.gz"
))[grepl("_[ACGT]+_[ACGT]+$", snp)]

print("Processing GWAS data ...")
gwas_df <- gwas_df[
    (beta != 0 | varbeta != 0) & !is.na(beta) & !is.na(varbeta) & !is.na(position),
    .(chromosome, beta, varbeta, position, snp)
][!duplicated(snp)]

print("Saving GWAS data ...")
for (chr_num in 1:22) {
    print(glue("Chromosome {chr_num}"))
    dir_path <- glue("/g/data/ei56/od8037/NewGenotypes/Coloc/BloodTraits/GWAS_Results/{condition}")
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
    }
    gwas_df_subset <- gwas_df[chromosome == glue("chr{chr_num}"), .(beta, varbeta, position, snp)]
    write_fst(gwas_df_subset, glue("{dir_path}/chr{chr_num}.fst"))
}

print("Finished!")