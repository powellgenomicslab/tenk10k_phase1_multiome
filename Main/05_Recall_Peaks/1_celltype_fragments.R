# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Load libraries
library(Signac)
library(glue)
library(data.table)
library(tidyverse)

# Get library for current run
cur_library <- as.character(args[1])
print(paste("Processing library:", cur_library))

# Find the directory containing the fragments file
base_dirs <- list.dirs("/g/data/fy54/data/atac/atac_count_outs/force_cells_test", recursive = FALSE, full.names = TRUE)
for (dir in base_dirs) {
    subdirs <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
    if (cur_library %in% subdirs) {
        frag_dir <- file.path(dir, cur_library, "outs/fragments.tsv.gz")
        break
    }
}

# Read fragments file
print(paste("Reading fragments file:", frag_dir))
frag_df <- fread(frag_dir) 
colnames(frag_df) <- c("V1", "V2", "V3", "V4", "V5", "V6")
frag_df <- frag_df %>% 
    dplyr::filter(V1 %in% paste0("chr", 1:22))

# Read metadata file
print("Reading metadata file")
meta_df <- fread(glue(
    "/g/data/ei56/jf1058/TenK10K/Multiome/data/annotated/annotated/tob_atac_{cur_library}_annotated_meta_new_ref_nonPBMC_excluded.csv"
)) %>% 
    dplyr::mutate(predicted.id = stringr::str_replace_all(predicted.id, " ", "_"))

# Load list of celltypes
celltype_list <- readLines("/g/data/ei56/od8037/TenK10K/PeakCalling/celltype_names.txt")

for (celltype in celltype_list) {
    print(paste("Processing celltype:", celltype))
    barcode_list <- meta_df %>%
        dplyr::filter(predicted.id == celltype) %>%
        dplyr::pull(V1)
    
    # Filter fragments
    filtered_df <- frag_df %>% 
        dplyr::filter(V4 %in% barcode_list) 
    
    # Append filtered rows to the output file
    output_file <- glue(
        "/g/data/ei56/od8037/TenK10K/PeakCalling/CelltypeFragments/{celltype}/{cur_library}.tsv"
    )

    if (nrow(filtered_df) > 0) {
        filtered_df %>%
            dplyr::mutate(V4 = paste0(cur_library, "_", V4)) %>%
            fwrite(output_file, sep = "\t", col.names = FALSE)
        
        # Compress the output file using gzip
        system(paste("gzip", output_file))
    }
}

message(paste("Finished processing:", cur_library))
