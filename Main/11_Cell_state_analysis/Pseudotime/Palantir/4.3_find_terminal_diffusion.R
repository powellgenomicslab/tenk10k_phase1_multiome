# Load libraries
library(data.table)
library(glue)
library(tidyverse)

# Load cell types
male_list <- fread("/g/data/fy54/proj_CH/MoChA/TenK10K_TOB_ATAC_231libraries_ID.tsv")[sex == 1]$individual %>% unique()
meta_df_list <- list()
for (cur_library in readLines("/g/data/ei56/od8037/TenK10K/PeakCalling/library_names.txt")) {
    print(paste("Processing library:", cur_library))
    meta_df <- fread(glue(
        "/g/data/ei56/jf1058/TenK10K/Multiome/data/annotated/annotated/tob_atac_{cur_library}_annotated_meta_new_ref_nonPBMC_excluded.csv"
    ))[predicted.id %in% c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Proliferating", "CD4 CTL", "Treg") & donor_id %in% male_list] %>%
        dplyr::mutate(barcode = paste0(V1, "_", cur_library)) %>%
        dplyr::select(barcode, predicted.id)
    if (nrow(meta_df) != 0) {
        meta_df_list[[cur_library]] <- meta_df
    }
}
full_meta_df <- rbindlist(meta_df_list)

# Load DM eigenvectors
df <- fread("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4/CSVs/DM_eigenvectors.csv", header = TRUE) %>% 
    dplyr::select(V1, 18) %>%
    dplyr::left_join(full_meta_df, by = c("V1" = "barcode")) %>%
    dplyr::filter(predicted.id == "Treg") 
colnames(df) <- c("barcode", "diffusion", "celltype")

top_barcode <- df %>%
    dplyr::arrange(diffusion) %>%
    dplyr::slice(1) %>%
    dplyr::pull(barcode)
print(top_barcode)
