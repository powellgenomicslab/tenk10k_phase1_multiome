source("/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/qval_functions.R")
library(data.table)
library(qvalue)
library(tidyverse)
library(glue)

celltype_list <- readLines("/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/celltype_names.txt")

for (celltype in celltype_list) {
    print(celltype)

    df <- as.data.frame(fread(glue(
        "/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/FinalResults/Raw/{celltype}.csv"
    ))) %>% dplyr::filter(Group == "Cauchy")
    df <- df[rowSums(is.na(df)) != ncol(df), ]

    df$Pvalue_noACAT <- mapply(
        function(p1, p2) get_CCT_pvalue(c(p1, p2)),
        df$Pvalue_Burden, df$Pvalue_SKAT
    )
    df$qvalue <- qvalue(df$Pvalue_noACAT)$qvalues

    regions_list <- list()
    for (chr_num in 1:22) {
        regions_list[[chr_num]] <- fread(glue(
            "/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Regions/{celltype}/chr{chr_num}.csv"
        )) %>% dplyr::mutate(Region = paste0(V2, "_", V3, "_", V4))
    }
    regions_df <- rbindlist(regions_list) %>% 
        dplyr::rename(Peak = V1) %>%
        dplyr::select(Region, Peak)

    other_cols <- setdiff(colnames(df), c("Region")) 
    final_df <- left_join(df, regions_df, by = "Region") %>%
        dplyr::select(all_of(c("Peak", other_cols))) %>%
        sort_peaks()
    fwrite(final_df, glue("/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/FinalResults/Processed/{celltype}.csv"))
}
