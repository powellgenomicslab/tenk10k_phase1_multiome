####
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)

# Read in pre-annotated data
# data = readRDS("/directflow/SCCGGroupShare/projects/josealq/multiome/results/2022-06-28_cell_type_annotation/annotation.RDS")

# load the RNA and ATAC data
wd = "/directflow/SCCGGroupShare/projects/josealq/multiome/data/ATAC/ES001_P1_1/outs/outs/"
counts <- Read10X_h5(paste0(wd,"filtered_feature_bc_matrix.h5"))
fragpath <- paste0(wd,"atac_fragments.tsv.gz")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg38"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

# Check
DefaultAssay(pbmc) <- "ATAC"
pbmc
pbmc[['ATAC']]
granges(pbmc)

# QC
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

p1 <- VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

saveRDS(p1, "multiome_plot_file.RDS")
ggsave(filename = "qc_metrics_ES001_P1_1.png", plot = p1, width = 20, height = 7)

p2 <- DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(filename = "qc_DensityScatter_ATAC_ES001_P1_1.png", plot = p2, width = 960*1.5, height = 720*1.5, units="px")

p3 <- DensityScatter(pbmc, x = 'nCount_RNA', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(filename = "qc_DensityScatter_RNA_ES001_P1_1.png", plot = p3, width = 960*1.5, height = 720*1.5, units="px")


# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 15000 &
    nCount_ATAC > 800 &
    nCount_RNA < 1400 &
    nCount_RNA > 80 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc




####
