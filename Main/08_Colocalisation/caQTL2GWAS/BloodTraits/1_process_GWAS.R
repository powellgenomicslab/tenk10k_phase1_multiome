# Load libraries
library(data.table)
library(glue)
library(fst)
library(stringi)

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
condition <- args[1]

# Load SNP mapping file
print("Loading SNP mapping file ...")
snp_map <- read_fst("/g/data/fy54/od8037/TenK10K/ColocNew/snp_id_map.fst", as.data.table = TRUE)
setkey(snp_map, swapped_style)

# Load GWAS data
print(glue("Loading GWAS data for {condition} ..."))
gwas_result_dir <- "/g/data/fy54/data/external_data/gwas/gymrek-ukbb-snp-gwas-catalogs_v6"
gwas_df <- fread(glue(
    "{gwas_result_dir}/white_british_{condition}_snp_gwas_results_hg38.tab.gz"
))[grepl("_[ACGT]_[ACGT]$", snp) & varbeta != 0]

# Re-map SNP IDs in GWAS data to PLINK format
print("Re-mapping SNP IDs in GWAS data to PLINK format ...")
setnames(gwas_df, "snp", "old_snp")
gwas_df[, snp := gsub("_", ":", sub("^chr", "", old_snp))]
gwas_df[snp_map, bim_style := i.bim_style, on = .(snp = swapped_style)]
gwas_df[!is.na(bim_style), `:=`(
    beta = -beta,
    snp = bim_style
)]
gwas_df <- na.omit(gwas_df[, .(chromosome, beta, varbeta, position, snp, p_value)])[!duplicated(snp)]

# Save the processed GWAS data for each chromosome
print("Saving GWAS data ...")
for (chr_num in 1:22) {
    print(glue("Chromosome {chr_num}"))
    dir_path <- glue("/g/data/fy54/od8037/TenK10K/ColocNew/BloodTraits/GWAS/{condition}")
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
    }
    gwas_df_subset <- gwas_df[chromosome == glue("chr{chr_num}"), .(beta, varbeta, position, snp)]
    write_fst(gwas_df_subset, glue("{dir_path}/chr{chr_num}.fst"))
}

print("Finished!")
