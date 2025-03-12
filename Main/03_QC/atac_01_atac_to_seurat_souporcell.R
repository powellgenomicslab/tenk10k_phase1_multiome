#### Convert ATAC to Seurat object for all libraries in TOB cohrot ####
# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Load libraries
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
# library(tidyverse)
library(dplyr)
library(glue)
library(purrr)    

# Set random seed.
set.seed(860)

# Set up working directory
setwd(paste0("/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/"))
# Set up the data directory
data_dir = "/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/data/"

# Read in the pool names
# pool_names = read.table(paste0(data_dir, "tob_atac_pool_names.txt"), header = FALSE)
# pool_names = pool_names$V1
folder_names <- system(paste0("ls -d ", data_dir, "S0*/"), intern = TRUE)
pool_names <- basename(folder_names)
pool_index <- as.numeric(args[1])
pool = pool_names[pool_index]
pool_modified <- sub("_[^_]+$", "", pool)
print(pool)

# Function to convert cellranger atac output to seurat object
atac_to_seurat <- function(pool, pool_individuals) {
    print(pool)
    cellranger_outs_path <- glue("/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/data/{pool}/")

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

    # add pool to the metadata
    seurat$library <- pool
    seurat$pool <- sub("_[^_]+$", "", pool)
    # get souporcell demuliplexing results and add to the seurat object
    soup_path <- glue("/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/demultiplexing/output/Souporcell/ATAC_{pool}/clusters.tsv")

    soup <- read.table(soup_path,header=T) %>%
        select(1:3) %>%
        left_join(pool_individuals)

    seurat <- AddMetaData(seurat, soup[,-1])

    return(seurat)
}

# Generate pool-individual pair list
t <- read.table(
	glue("/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/data/vcf/{pool_modified}_sample_list.txt"),
	header = FALSE
    )
t = cbind(assignment = 0:7, individual = t$V1)
t = as.data.frame(t)

# read in the library data + demuxafy info
pbmc <- atac_to_seurat(pool, t)
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
saveRDS(pbmc, paste0("./data/QCed/tob_atac_",pool,"_raw.RDS"))
# Save meta data
meta <- pbmc@meta.data
write.csv(meta, paste0("./data/QCed/tob_atac_", pool, "_meta.csv"))

