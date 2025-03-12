#### Regular QC on the ATAC object ####
# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratDisk)
  library(patchwork)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
  library(dplyr)
  library(glue)
  library(purrr)
  library(data.table)
  library(presto)
})

# Set random seed.
set.seed(860)
# Set up working directory
setwd(paste0("/g/data/ei56/ax3061/proj/tenk10k/caQTL"))
# Set up the data directory
data_dir = "/g/data/fy54/data/atac/atac_count_outs/"

# Read in the lib names
lib <- as.character(args[1])
folder_names <- system(paste0("find /g/data/fy54/data/atac/atac_count_outs/ -maxdepth 2 -type d -name '", lib, "'"), intern = TRUE)
# If no folder is found, search for flexible alternative patterns
if (length(folder_names) == 0) {
  base_id <- sub("_(\\d+)$", "", lib)  # Extracts "S0228" from "S0228_2"
  num_suffix <- sub("^.*_(\\d+)$", "\\1", lib)  # Extracts "2" from "S0228_2"

  # Generate alternative patterns
  alt_pattern_1 <- paste0(base_id, "_R_", num_suffix)  # "S0228_R_2"
  alt_pattern_2 <- paste0(base_id, "_", num_suffix, "_R")  # "S0228_2_R"

  folder_names <- system(paste0("find /g/data/fy54/data/atac/atac_count_outs/ -maxdepth 2 -type d \\( -name '", alt_pattern_1, "' -o -name '", alt_pattern_2, "' \\)"), intern = TRUE)
}
lib_names <- basename(folder_names)
lib_modified <- sub("_[^_]+$", "", lib)
print(lib)


# Read merged Seurat object
pbmc <- readRDS(paste0("./data/raw/tob_atac_",lib,"_raw.RDS"))
pbmc

# Remove doublets and unassigned cells
Idents(pbmc) <- "donor_id"
print(paste0(ncol(pbmc)," cells before removing doublets"))
pbmc <- subset(pbmc, idents = c("doublet", "unassigned"), invert = TRUE)
print(paste0(ncol(pbmc)," cells after removing doublets"))
Idents(pbmc) <- "library"
# Distribution of each QC metric
p1 <- VlnPlot(
  object = pbmc,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.08,
  ncol = 5
)
ggsave(paste0("./output/plots/ATAC_QC_violin_", lib, ".pdf"), p1, width = 13, height = 5)

# nCount_peaks vs TSS.enrichment
p2 <- DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(paste0("./output/plots/ATAC_QC_DensityScatter_", lib, ".pdf"), p2, width = 8, height = 6)

# Fragment length periodicity for all the cells
p3 <- FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
ggsave(paste0("./output/plots/ATAC_QC_FragmentHistogram_", lib, ".pdf"), p3, width = 8, height = 6)

## QC for ATAC
# RNA data is annotated before any QCs
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 5000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)
print("After removing low quality cells:")
pbmc

# Save QCed object
saveRDS(pbmc,paste0("./data/QCed/tob_atac_",lib,"_QCed.rds"))

# Normalization and linear dimensional reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

# Visualize the top features
p4 <- DepthCor(pbmc)
ggsave(paste0("./output/plots/ATAC_QC_DepthCor_", lib, ".pdf"), p4, width = 8, height = 6)

# Non-linear dimension reduction and clustering
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

p5 <- DimPlot(object = pbmc, label = TRUE) + NoLegend()
ggsave(paste0("./output/plots/ATAC_QC_UMAP_", lib, ".pdf"), p5, width = 8, height = 6)

# Create a gene activity matrix
gene.activities <- GeneActivity(pbmc)
# Add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)

# Normliaze RNA data
# Shall I use SCTransform?
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
)

# Analyse RNA slot
DefaultAssay(pbmc) <- 'RNA'

p6 <- FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
ggsave(paste0("./output/plots/ATAC_QC_RNA_FeaturePlot_", lib, ".pdf"), p6, width = 12, height = 6)

# Save QCed object with imputed RNA data
# saveRDS(pbmc,paste0("./data/QCed/tob_atac_",lib,"_QCed_with_imputed_RNA.rds"))

# Save Session Info
sessionInfo()
