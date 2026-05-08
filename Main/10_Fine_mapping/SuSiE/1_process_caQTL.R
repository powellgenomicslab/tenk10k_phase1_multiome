library(data.table)
library(glue)
library(arrow)
library(fst)

args <- commandArgs(trailingOnly = TRUE)
celltype_name <- args[1]
print(paste("Processing celltype:", celltype_name))

caqtl_dir <- glue("/g/data/ei56/od8037/TenK10K/caQTLNew/Results/{celltype_name}")
save_dir <- glue("/g/data/fy54/od8037/TenK10K/FineMapping/caQTLResults/1Mb/{celltype_name}")
if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
}

for (chr_num in 1:22) {
    print(glue("Processing chromosome {chr_num}..."))

    # Load caQTL results
    print("Loading caQTL results...")
    caqtl_df_path <- glue("{caqtl_dir}/TenK10K.cis_qtl_pairs.chr{chr_num}.parquet")
    if (!file.exists(caqtl_df_path)) {
        print(glue("File {caqtl_df_path} does not exist. Skipping chromosome {chr_num}."))
        next
    }
    caqtl_df <- as.data.table(read_parquet(caqtl_df_path))[
        , .(phenotype_id, variant_id, slope, slope_se)
    ]
    caqtl_df[, Z := slope / slope_se]

    # Save processed caQTL results
    print("Saving processed caQTL results...")
    write_fst(caqtl_df, glue("{save_dir}/chr{chr_num}.fst"))
}
print("Finished!")
