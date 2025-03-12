#### Convert ATAC to Seurat object for all libraries in TOB cohrot ####
# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Load libraries
suppressPackageStartupMessages({
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
# library(tidyverse)
library(dplyr)
library(glue)
library(purrr)    
})

# Set random seed.
set.seed(860)

# Set up working directory
setwd(paste0("/g/data/ei56/ax3061/proj/tenk10k/caQTL"))
# Set up the data directory
data_dir = "/g/data/fy54/data/atac/atac_count_outs/force_cells_test/"

# Read in the pool names
pool <- as.character(args[1])
folder_names <- system(paste0("find /g/data/fy54/data/atac/atac_count_outs/force_cells_test/ -maxdepth 2 -type d -name '", pool, "'"), intern = TRUE)
# If no folder is found, search for flexible alternative patterns
if (length(folder_names) == 0) {
  base_id <- sub("_(\\d+)$", "", pool)  # Extracts "S0228" from "S0228_2"
  num_suffix <- sub("^.*_(\\d+)$", "\\1", pool)  # Extracts "2" from "S0228_2"
  
  # Generate alternative patterns
  alt_pattern_1 <- paste0(base_id, "_R_", num_suffix)  # "S0228_R_2"
  alt_pattern_2 <- paste0(base_id, "_", num_suffix, "_R")  # "S0228_2_R"
  
  folder_names <- system(paste0("find /g/data/fy54/data/atac/atac_count_outs/force_cells_test/ -maxdepth 2 -type d \\( -name '", alt_pattern_1, "' -o -name '", alt_pattern_2, "' \\)"), intern = TRUE)
}
pool_names <- basename(folder_names)
pool_modified <- sub("_[^_]+$", "", pool)
print(pool)

# Function to convert cellranger atac output to seurat object
atac_to_seurat <- function(pool) {
    print(pool)
    cellranger_outs_path <- glue(folder_names,"/outs/")

    # read the cellranger atac outs into seurat object
    counts <- Read10X_h5(
        filename = file.path(cellranger_outs_path, "filtered_peak_bc_matrix.h5")
    )

    metadata <- read.csv(
        file = file.path(cellranger_outs_path, "singlecell.csv"),
        header = TRUE,
        row.names = 1
    )

    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),
        fragments = file.path(cellranger_outs_path, "fragments.tsv.gz"),
        min.cells = 10,
        min.features = 200
    )

    seurat <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks",
        meta.data = metadata
    )
    rm(chrom_assay,counts,metadata)
    gc()
    # add pool to the metadata
    seurat$library <- pool
    seurat$pool <- sub("_[^_]+$", "", pool)
    # get souporcell demuliplexing results and add to the seurat object
    vireo_path <- glue("/g/data/ei56/ax3061/proj/tenk10k/caQTL/demultiplexing/output/Vireo/force_cell/{pool}/donor_ids.tsv")

    vireo <- read.table(vireo_path, header = T) %>%
        select(1:3)

    seurat <- AddMetaData(seurat, vireo[,-1])

    return(seurat)
}

# read in the library data + demuxafy info
pbmc <- atac_to_seurat(pool)
gc()
#### Calculate QC metrics ####
pbmc
pbmc[['peaks']]
Idents(pbmc) <- "library"

granges(pbmc)
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]

# get gene annotations for hg38
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# Add the gene information to the object
Annotation(pbmc) <- annotations

print("Computing QC metrics")
# Compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# Compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc)

# Add fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# Add blacklist ratio
pbmc$blacklist_ratio <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)
# Add nucleosome group
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# Save the object
pbmc
saveRDS(pbmc, paste0("./data/raw/tob_atac_",pool,"_raw.RDS"))
# Save meta data
meta <- pbmc@meta.data
write.csv(meta, paste0("./data/raw/tob_atac_", pool, "_meta.csv"))

