# Load libraries
library(tidyverse)
library(glue)
library(data.table)

# Set path to metadata
annotated_dir <- "/g/data/ei56/jf1058/TenK10K/Multiome/data/annotated/annotated"

# Get list of libraries
library_list <- readLines("/g/data/fy54/od8037/TenK10K/WASP_16/library_names.txt")

# Get metadata for chosen libraries
full_meta <- data.frame()
for (cur_library in library_list) {
    print(cur_library)
    meta_df <- fread(glue(
        "{annotated_dir}/tob_atac_{cur_library}_annotated_meta_new_ref_nonPBMC_excluded.csv"
    )) %>% 
        dplyr::mutate(
            predicted.id = gsub(" ", "_", predicted.id),
            sinto_group = paste0(library, "-", predicted.id, "-", donor_id)
        ) %>%
        dplyr::rename(barcode = V1, celltype = predicted.id) 

    full_meta <- rbind(full_meta, meta_df)
}

# Generate 'donor_names.txt' file
donor_names <- sort(unique(full_meta$donor_id))
writeLines(donor_names, "/g/data/fy54/od8037/TenK10K/WASP_16/donor_names.txt")

# Generate 'celltype_names.txt' file
celltype_list <- sort(unique(full_meta$celltype))
writeLines(celltype_list, "/g/data/fy54/od8037/TenK10K/WASP_16/celltype_names.txt")

# Generate combinations file
full_meta %>%
    dplyr::select(library, celltype, donor_id) %>%
    dplyr::arrange(library, celltype, donor_id) %>%
    dplyr::distinct() %>%
    write.table(
        file = glue("/g/data/fy54/od8037/TenK10K/WASP_16/Combinations/combinations.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

