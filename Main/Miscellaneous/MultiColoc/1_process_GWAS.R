# Load libraries
library(data.table)
library(glue)
library(fst)

# Define named list of GWAS file paths
new_gwas_dir <- "/g/data/ei56/od8037/NewGenotypes/Coloc/caQTL2GWAS/NewDiseases/Raw_GWAS"
gwas_paths <- list(
    "ms" = "/g/data/fy54/analysis/tenk10k-causal/resources/pipeline_ma/ms.ma",
    "t1dm" = "/g/data/fy54/analysis/tenk10k-causal/resources/pipeline_ma/t1dm.ma",
    "asthma" = "/g/data1a/ei56/as8574/analysis/TenK10K_SMR/inputs/ma/asthma.ma",
    "cd" = glue("{new_gwas_dir}/CD_EAS_EUR_meta_Liu_et_al_QCed.txt.gz"),
    "ibd" = glue("{new_gwas_dir}/IBD_EAS_EUR_meta_Liu_et_al_QCed.txt.gz"),
    "uc" = glue("{new_gwas_dir}/UC_EAS_EUR_meta_Liu_et_al_QCed.txt.gz")
)
# condition_list <- names(gwas_paths)
condition_list <- c("ibd")

# Load SNP mapping file
print("Loading SNP mapping file ...")
snp_map <- read_fst("/g/data/fy54/od8037/TenK10K/ColocNew/snp_id_map.fst", as.data.table = TRUE)
setkey(snp_map, swapped_style)

# Loop through and process each GWAS file
for (condition_name in condition_list) {
    print(glue("Processing GWAS data for {condition_name}..."))
    out_dir <- glue("/g/data/fy54/od8037/TenK10K/MultiColoc/caQTL2GWAS/GWAS/{condition_name}")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    gwas_file <- gwas_paths[[condition_name]]
    df <- fread(gwas_file)[nchar(A1) == 1 & nchar(A2) == 1]

    if (condition_name %in% c("cd", "ibd", "uc")) {
        df <- df[CHR %in% as.character(1:22)]
        df[, SNP := gsub("^chr", "", SNP)]
    }

    # Re-map SNP IDs in GWAS data to PLINK format
    df[snp_map, bim_style := i.bim_style, on = .(SNP = swapped_style)]
    df[!is.na(bim_style), `:=`(
        b = -b,
        SNP = bim_style
    )]

    # Update columns for final output
    df[, c("snp", "beta", "varbeta", "chr", "position") := {
        splits <- tstrsplit(SNP, ":", fixed = TRUE, keep = 1:2, type.convert = TRUE)
        list(SNP, b, se^2, splits[[1]], splits[[2]])
    }]
    df <- na.omit(df[, .(chr, beta, varbeta, position, snp, P)])

    # Save each chromosome separately
    for (chr_num in 1:22) {
        print(glue("Processing chromosome {chr_num} ..."))
        df_subset <- df[
            chr == chr_num, .(beta, varbeta, position, snp, P)
        ][!duplicated(snp)]
        write_fst(df_subset, glue("{out_dir}/chr{chr_num}.fst"))
    }
}
