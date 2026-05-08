library(tidyverse)
library(data.table)
library(glue)

celltype_list <- readLines("/g/data/fy54/od8037/Brenner/SAIGEQTL/Phenotypes/celltype_names.txt")
for (celltype in celltype_list) {
    if (celltype == "CD8_Proliferating" || celltype == "ILC") { next }
    tensor_peak_list <- readLines(glue(
        "/g/data/ei56/od8037/TenK10K/caQTLNew/PeakLists/{celltype}.txt"
    ))
    saige_peak_list <- readLines(glue(
        "/g/data/fy54/od8037/Brenner/SAIGEQTL/Phenotypes/ValidPeaks/{celltype}.txt"
    ))
    final_list <- intersect(tensor_peak_list, saige_peak_list)
    
    df_list[[celltype]] <- data.frame(
        celltype = celltype,
        tensor = length(tensor_peak_list),
        saige = length(saige_peak_list),
        intersect = length(final_list)
    )
    writeLines(final_list, glue("/g/data/fy54/od8037/Brenner/SAIGEQTL/Phenotypes/FinalPeaks/{celltype}.txt"))
}
