#### Convert ATAC to Seurat object for all libraries in TOB cohort ####
# 
# SCRIPT PURPOSE:
# This script converts Cell Ranger ATAC output files into Seurat objects and
# performs initial quality control for ATAC-seq data from the TOB cohort.
#
# KEY FUNCTIONS:
# 1. Load Cell Ranger ATAC outputs (peak-barcode matrix, metadata, fragments)
# 2. Create Seurat object with ChromatinAssay
# 3. Integrate Vireo demultiplexing results (donor assignments)
# 4. Calculate QC metrics (nucleosome signal, TSS enrichment, reads in peaks, blacklist ratio)
# 5. Save processed Seurat object and metadata
#
# USAGE:
# Rscript atac_01_atac_to_seurat.R [pool_name]
#
# INPUT:
# - pool_name: Pool identifier (e.g., "S0228_2")
# - Cell Ranger ATAC outputs from /g/data/fy54/data/atac/atac_count_outs/force_cells_test/
# - Vireo demultiplexing results
#
# OUTPUT:
# - ./data/raw/tob_atac_[pool]_raw.RDS - Seurat object with QC metrics
# - ./data/raw/tob_atac_[pool]_meta.csv - Cell metadata

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

# Read in the pool names (from command line argument)
pool <- as.character(args[1])
# Search for the pool directory in the Cell Ranger ATAC output folder
folder_names <- system(paste0("find /g/data/fy54/data/atac/atac_count_outs/force_cells_test/ -maxdepth 2 -type d -name '", pool, "'"), intern = TRUE)
# If no folder is found, search for flexible alternative naming patterns
# This handles different naming conventions used across batches
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
# Parameters:
#   pool: Pool identifier for locating Cell Ranger output files
# Returns:
#   Seurat object with ATAC data, metadata, and demultiplexing results
atac_to_seurat <- function(pool) {
    print(pool)
    cellranger_outs_path <- glue(folder_names,"/outs/")

    # Read the Cell Ranger ATAC outputs into Seurat object
    # 1. Load the filtered peak-barcode count matrix
    counts <- Read10X_h5(
        filename = file.path(cellranger_outs_path, "filtered_peak_bc_matrix.h5")
    )

    # 2. Load single-cell metadata (quality metrics from Cell Ranger)
    metadata <- read.csv(
        file = file.path(cellranger_outs_path, "singlecell.csv"),
        header = TRUE,
        row.names = 1
    )

    # 3. Create ChromatinAssay with quality filters
    # min.cells = 10: Only keep peaks detected in at least 10 cells
    # min.features = 200: Only keep cells with at least 200 features
    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),
        fragments = file.path(cellranger_outs_path, "fragments.tsv.gz"),
        min.cells = 10,
        min.features = 200
    )

    # 4. Create Seurat object with the ChromatinAssay
    seurat <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks",
        meta.data = metadata
    )
    # Clean up memory
    rm(chrom_assay,counts,metadata)
    gc()
    # Add pool and library identifiers to the metadata
    seurat$library <- pool
    seurat$pool <- sub("_[^_]+$", "", pool)
    
    # 5. Integrate Vireo demultiplexing results
    # Vireo assigns cells to individual donors within each pool
    vireo_path <- glue("/g/data/ei56/ax3061/proj/tenk10k/caQTL/demultiplexing/output/Vireo/force_cell/{pool}/donor_ids.tsv")

    vireo <- read.table(vireo_path, header = T) %>%
        select(1:3)

    # Add donor assignment metadata to Seurat object
    seurat <- AddMetaData(seurat, vireo[,-1])

    return(seurat)
}

# Execute the conversion function
# read in the library data + demuxafy info
pbmc <- atac_to_seurat(pool)
gc()

#### Calculate QC metrics ####
pbmc
pbmc[['peaks']]
Idents(pbmc) <- "library"

# Filter to standard chromosomes only
# This removes mitochondrial DNA and non-standard contigs
granges(pbmc)
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]

# Get gene annotations for hg38 genome build from Ensembl
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# Add the gene information to the object
# This enables downstream gene-centric analyses
Annotation(pbmc) <- annotations

print("Computing QC metrics")
# Compute nucleosome signal score per cell
# Nucleosome signal represents the ratio of mononucleosome to nucleosome-free fragments
# Higher values (>4) indicate lower quality cells
pbmc <- NucleosomeSignal(object = pbmc)

# Compute TSS enrichment score per cell
# TSS enrichment measures the enrichment of ATAC-seq signal at transcription start sites
# Higher values indicate better data quality
pbmc <- TSSEnrichment(object = pbmc)

# Add fraction of reads in peaks
# Higher percentages indicate better signal-to-noise ratio
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# Add blacklist ratio
# Blacklisted regions are known problematic genomic regions (e.g., repetitive elements)
# Lower values are better
pbmc$blacklist_ratio <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)
# Add nucleosome group classification
# Cells with NS > 4 are typically lower quality
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# Save the processed Seurat object and metadata
pbmc
saveRDS(pbmc, paste0("./data/raw/tob_atac_",pool,"_raw.RDS"))
# Save meta data
meta <- pbmc@meta.data
write.csv(meta, paste0("./data/raw/tob_atac_", pool, "_meta.csv"))

