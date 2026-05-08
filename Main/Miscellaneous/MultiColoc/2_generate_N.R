library(tidyverse)
library(data.table)
library(glue)

eqtl_path <- "/g/data/fy54/results/eqtl/december24_freeze"
celltype_list <- readLines("/g/data/fy54/od8037/TenK10K/MultiColoc/eQTL2GWAS/celltypes.txt")

df_list <- list()
for (celltype in celltype_list) {
    N <- fread(glue(
        "{eqtl_path}/{celltype}/{celltype}_common_all_cis_raw_pvalues.tsv"
    ), nrows = 1)$N
    df_list[[celltype]] <- data.table(
        celltype_name = celltype,
        num_donors = N
    )
}
df <- rbindlist(df_list)
fwrite(df, "/g/data/fy54/od8037/TenK10K/MultiColoc/eQTL2GWAS/num_donors.csv")
