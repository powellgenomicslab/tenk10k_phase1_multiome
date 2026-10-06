

#### Creat coverage plot for a gene of interest with a SNP highlighted
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratData)
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

# Set up working directory
setwd(paste0("/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/"))
# Set up the data directory
data_dir = "/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/data/"

# Read in a library
# load the pre-processed atac data
lib = "S0228_2"
pbmc = readRDS(paste0(data_dir,"QCed/first_56_libraries/tob_atac_",lib,"_annotated.rds"))

pbmc <- readRDS("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/SMR/coloc_SMR_compare/analysis/ETS2_plot/tob_atac_S0220_1_annotated.rds")

# Find differentially accessible peaks between cell types
# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'ATAC'
Idents(pbmc) <- "predicted.id"

# get gene annotations for hg38
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# Add the gene information to the object
Annotation(pbmc) <- annotations

pbmc@meta.data$predicted.id = gsub("_", " ", pbmc@meta.data$predicted.id)
pbmc <- subset(pbmc, idents = c("Eryth", "Platelet"), , invert = TRUE)
## Choose the gene of interest
colp = fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/colour_palette_table.tsv",header=T)
colp = fread("/g/data/ei56/jf1058/TenK10K/Multiome/scGLUE/full_cohort/colour_palette_table.tsv",header=T)

colp = colp[colp$cell_type %in% unique(pbmc@meta.data$predicted.id),]
# Sort the cell types
Idents(pbmc) <- "predicted.id"
Idents(pbmc) <- factor(pbmc@active.ident, levels = colp$cell_type)

custom_colors <- colp$color
names(custom_colors) = colp$cell_type

# Gene ETS2
gene = "MARCH8"
SNP_CHR = "chr10"
SNP_BP = 45592613
# roi = "chr10-45592479-45592785"
roi = "chr10-45400000-45650000"
window_size = - 45400000 + 45650000
SNP_highlight_size = window_size / 1000
# Create a GRanges object for the region of interest
regions_highlight <- GRanges( 
  seqnames = Rle(c(SNP_CHR)), 
  ranges = IRanges( start = c(SNP_BP - SNP_highlight_size), 
                    end = c(SNP_BP + SNP_highlight_size)), 
  strand = Rle("*") )
regions_highlight$color <- "black"

p3 <- CoveragePlot(
  object = subset(x = pbmc, idents = c('NK', 'cDC2')),
  region = roi,
  region.highlight = regions_highlight,
  annotation = FALSE,
  peaks = FALSE
) 

cov_plot3 <- p3 + scale_fill_manual(values = custom_colors, breaks = colp$cell_type) +
  theme(axis.title.y = element_text(size = 7))

gene_plot <- AnnotationPlot(
  object = pbmc,
  region = roi
)  + theme(axis.title.y = element_text(size = 9))
# gene_plot <- gene_plot + scale_color_manual(values = c("darkblue", "darkblue"))
gene_plot <- gene_plot + scale_color_manual(values = c("darkblue", "#006400"))

peak_plot <- PeakPlot(
  object = pbmc,
  region = roi
) + theme(axis.title.y = element_text(size = 9), axis.title.x = element_text(size = 9))

p_comb3 <- CombineTracks(
  plotlist = list(cov_plot3, gene_plot, peak_plot),
  heights = c(3, 2, 1),
  widths = c(10)
)

ggsave(paste0("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_opp_direction/atac_tracks/3.png"), p_comb3, height = 2.5, width = 9)
