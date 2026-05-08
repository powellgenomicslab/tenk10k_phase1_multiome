# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(arrow)

# Generate cell-level covariates
pool_list <- readLines("/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/pool_list.txt")
meta_list <- list()
for (pool in pool_list) {
    print(glue("Processing {pool}"))
    meta_df <- fread(glue(
        "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/ATAC_Final/output/QCed/tob_atac_{pool}_annotated_meta_new_ref_nonPBMC_excluded.csv"
    )) %>%
        dplyr::mutate(barcode = glue("{library}_{V1}")) %>%
        dplyr::select(barcode, library, nCount_peaks)
    meta_list[[pool]] <- meta_df
}

meta_df_full <- rbindlist(meta_list) %>%
    dplyr::mutate(repeat_num = as.integer(sub(".*_(\\d+)$", "\\1", library))) %>%
    dplyr::select(barcode, repeat_num, nCount_peaks)

write_parquet(
    meta_df_full,
    "/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Covariates/cell_covariates.parquet"
)

# Generate donor-level covariates
covariates_df <- fread("/directflow/SCCGGroupShare/projects/oscdon/FirstPipeline/TensorQTL/unzipped_covariates/CD4_NC_covar_peer_factors_PF0.txt")

covariates_df <- t(covariates_df[1:8, ]) 
colnames(covariates_df) <- covariates_df["id", ]
covariates_df <- covariates_df[-1, ] %>% as.data.frame()
covariates_df[] <- lapply(covariates_df, as.numeric)
covariates_df <- rownames_to_column(covariates_df, "donor_id")
colnames(covariates_df) <- c(
    "donor_id", "sex", "geno_pc1", "geno_pc2", "geno_pc3", 
    "geno_pc4", "geno_pc5", "geno_pc6", "age"
)

write_parquet(
    covariates_df,
    "/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Covariates/donor_covariates.parquet"
)


