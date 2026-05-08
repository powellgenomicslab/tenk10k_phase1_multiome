library(tidyverse)
library(glue)
library(data.table)

celltype_list <- readLines("/g/data/ei56/od8037/Final_caQTL/celltype_names.txt")
for (celltype in celltype_list) {
    tensor_df <- fread(glue(
        "/g/data/fy54/angxue/brenner_backup/ScratchGeneral/proj/caQTL/genotype/{celltype}_top_caQTL_100kb.csv"
    )) %>%
        dplyr::mutate(peak = gsub(":", "_", phenotype_id)) %>%
        dplyr::mutate(
            peak = gsub("-", "_", peak),
            top_Tensor = variant_id
        ) %>%
        dplyr::select(peak, top_Tensor)
    dir.create(glue("/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/InputTSVs/{celltype}"), showWarnings = FALSE)

    for (chr_num in 1:22) {
        print(glue("Processing {celltype} chromosome {chr_num}..."))
        region_df <- fread(glue(
            "/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Regions/{celltype}/chr{chr_num}.csv"
        ))
        colnames(region_df) <- c("peak", "chr", "start", "end")

        comb_df <- dplyr::left_join(region_df, tensor_df, by = "peak") %>% 
            drop_na() %>%
            distinct() 
        
        if (nrow(comb_df) > nrow(region_df)) {
            print(glue("Warning: Found more rows in combined dataframe than original regions for {celltype} chromosome {chr_num}."))
        }

        fwrite(
            comb_df,
            file = glue("/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/InputTSVs/{celltype}/chr{chr_num}.csv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
    }
}

print("Finished!")
