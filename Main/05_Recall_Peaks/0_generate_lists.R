# Load required libraries
library(data.table)
library(tidyverse)
library(glue)

# Set base directory for annotated RDS and metadata files
annotated_dir <- "/g/data/ei56/jf1058/TenK10K/Multiome/data/annotated/annotated"

# Get list of RDS files and extract library IDs 
rds_list <- list.files(annotated_dir, pattern = "_annotated.rds")
library_list <- sub(".*(S\\d{4}_\\d).*", "\\1", rds_list)

# Initialise vector to collect all unique cell types
all_celltypes <- c()

# Loop over each library to extract cell type labels from metadata CSVs
for (lib in library_list) {
  
  # Build path to corresponding metadata file
  meta_path <- glue("{annotated_dir}/tob_atac_{lib}_annotated_meta_new_ref_nonPBMC_excluded.csv")
  
  # Read metadata, clean cell type labels, and extract unique cell types
  celltypes <- fread(meta_path) %>%
    mutate(predicted.id = str_replace_all(predicted.id, " ", "_")) %>%
    pull(predicted.id) %>%
    unique()
  
  # Add new cell types to master list
  all_celltypes <- unique(c(all_celltypes, celltypes))
}

# Write the combined list of unique cell types to a file
output_path <- "/g/data/ei56/od8037/TenK10K/PeakCalling/celltype_names.txt"
writeLines(sort(all_celltypes), output_path)

# Generate list of libraries
rds_list <- list.files("/g/data/ei56/jf1058/TenK10K/Multiome/data/annotated/annotated", pattern = "*_annotated.rds")
library_list <- sub(".*(S\\d{4}_\\d).*", "\\1", rds_list)
writeLines(library_list, "/g/data/ei56/od8037/TenK10K/PeakCalling/library_names.txt")

# Create result folder for each celltype
for (celltype in all_celltypes) {
  dir.create(glue("/g/data/ei56/od8037/TenK10K/PeakCalling/CelltypeFragments/{celltype}"), recursive = TRUE, showWarnings = FALSE)
}
