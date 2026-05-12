# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Load libraries
library(Seurat)
library(SeuratDisk)
library(Signac)
library(glue)
library(tidyverse)
library(data.table)

# Get library for current run
cur_library <- as.character(args[1])
print(paste("Processing library:", cur_library))
save_dir <- "/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/Objects"

# Load data and metadata
print("Loading data...")
data <- readRDS(glue(
    "/g/data/ei56/ax3061/proj/tenk10k/caQTL/data/QCed/tob_atac_{cur_library}_QCed_with_imputed_RNA.rds"
))
print("Loading metadata...")
meta_df <- fread(glue(
    "/g/data/ei56/jf1058/TenK10K/Multiome/data/annotated/annotated/tob_atac_{cur_library}_annotated_meta_new_ref_nonPBMC_excluded.csv"
))[predicted.id %in% c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Proliferating", "CD4 CTL", "Treg")]

if (nrow(meta_df) == 0) {
    stop("No valid cells found.")
}

print("Subsetting data...")
data_filtered <- subset(data, cells = meta_df$V1)

# Keep only RNA assay
print("Subsetting to RNA assay...")
data_rna <- DietSeurat(
    data_filtered,
    assays = "RNA",
    layers = c("counts"),
    features = NULL,
    dimreducs = NULL,
    graphs = NULL
)

# Update metadata
print("Updating metadata...")
meta <- data_rna@meta.data %>%
    dplyr::mutate(barcode = colnames(data_rna)) %>%
    dplyr::left_join(
        meta_df %>%
            dplyr::select(V1, predicted.id) %>%
            dplyr::rename(barcode = V1, cell_type = predicted.id),
        by = "barcode"
    ) %>%
    dplyr::mutate(barcode = paste0(barcode, "_", cur_library)) %>%
    dplyr::select(barcode, library, pool, donor_id, nCount_RNA, nFeature_RNA, cell_type)

data_rna@meta.data <- meta %>% tibble::column_to_rownames("barcode")

# Save as h5Seurat and convert to h5ad
print("Saving as h5ad...")
SaveH5Seurat(
    data_rna,
    filename = glue("{save_dir}/{cur_library}.h5Seurat"),
    overwrite = TRUE
)
print("Converting to h5ad...")
Convert(
    glue("{save_dir}/{cur_library}.h5Seurat"),
    dest = "h5ad",
    overwrite = TRUE
)
file.remove(glue("{save_dir}/{cur_library}.h5Seurat"))
