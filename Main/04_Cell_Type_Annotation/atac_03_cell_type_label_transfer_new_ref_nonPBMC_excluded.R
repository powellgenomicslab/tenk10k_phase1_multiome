####
# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratData)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

# Set random seed.
set.seed(860)
# Set up working directory
setwd("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/ATAC_Final/")
# Set up the data directory
data_dir = "./data/QCed/"
          
file_name_list <- list.files(path = data_dir, 
                             pattern = "QCed\\.rds$", 
                             full.names = FALSE)
indices <- str_extract(file_name_list, "S\\d+_\\d+")
indices <- sort(indices)

lib <- indices[args[1]]

# load the pre-processed atac data
pbmc = readRDS(paste0(data_dir,"tob_atac_",lib,"_QCed.rds"))
# Rename the assay "peaks" to "ATAC"
pbmc <- RenameAssays(pbmc, peaks = "ATAC")

# load the pre-processed multiome data
pbmc.multi <- readRDS("/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/data/multiome_ref/pbmc_multiome_combined.rds")
DefaultAssay(pbmc.multi) <- "ATAC"
# Remove three non PBMC cell types
Idents(pbmc.multi) <- "predicted.id"
pbmc.multi <- subset(pbmc.multi, idents = c("Doublet", "Eryth", "Platelet"), invert = TRUE)

## ensure that the same features are measured in each dataset
# quantify multiome peaks in the scATAC-seq dataset
fragments <- CreateFragmentObject(
  path = paste0("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/ATAC_Final/data/fragments/", lib, "/fragments_nostrand.tsv.gz"),
  cells = colnames(pbmc), 
  validate.fragments = TRUE
)
Fragments(pbmc) <- NULL
Fragments(pbmc) <- fragments

counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = granges(pbmc.multi),
  cells = colnames(pbmc)
)

# add new assay with multiome peaks
pbmc[['ATAC']] <- CreateChromatinAssay(
  counts = counts,
  fragments = Fragments(pbmc)
)

rm(counts)
gc()

# compute LSI
DefaultAssay(pbmc) <- "ATAC"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 10)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

## merge the multiome and scATAC datasets together
# first add dataset-identifying metadata
pbmc$dataset <- "ATAC"
pbmc.multi$dataset <- "Multiome"

# merge
pbmc.combined <- merge(pbmc, pbmc.multi)

# process the combined dataset
pbmc.combined <- FindTopFeatures(pbmc.combined, min.cutoff = 10)
pbmc.combined <- RunTFIDF(pbmc.combined)
pbmc.combined <- RunSVD(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi", dims = 2:50)
p1 <- DimPlot(pbmc.combined, group.by = "dataset")

## integration anchors between the two datasets
# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(pbmc.multi, pbmc),
  anchor.features = rownames(pbmc.multi),
  reduction = "rlsi",
  dims = 2:50
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
)

rm(pbmc.combined)
gc()

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:50)
p2 <- DimPlot(integrated, group.by = "dataset")

p12 <- (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
print("Saving co-embedding plot")
ggsave(paste0("./output/plots/multiome_atac_integration_",lib , "_new_ref_nonPBMC_excluded.png"), p12, width = 12, height = 6)

## Reference mapping
# compute UMAP and store the UMAP model
# pbmc.multi <- RunUMAP(pbmc.multi, reduction = "lsi", dims = 2:30, return.model = TRUE)
pbmc.multi <- RunUMAP(pbmc.multi, reduction = "lsi", dims = 2:50, 
                      n.neighbors = 30, min.dist = 0.2, return.model = TRUE)

# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = pbmc.multi,
  query = pbmc,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:50
)

# map query onto the reference dataset
pbmc <- MapQuery(
  anchorset = transfer.anchors,
  reference = pbmc.multi,
  query = pbmc,
  refdata = pbmc.multi$predicted.id,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

## Visualisation
p3 <- DimPlot(pbmc.multi, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference")
p4 <- DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle(lib)

p34 <- p3 | p4
print("Saving annotation plot")
ggsave(paste0("./output/plots/ATAC_QC_Celltype_Annotations_", lib, "_new_ref_nonPBMC_excluded.png"), p34, width = 12, height = 6)

# Check annotation results
# table(pbmc$predicted.id)
# table(pbmc$predicted.id) / nrow(pbmc) * 100

# Clean 
rm(pbmc.multi, integration.anchors, transfer.anchors)
gc()

## save the annotated dataset
saveRDS(pbmc, paste0("./output/QCed/tob_atac_", lib, "_annotated.rds"))
meta = pbmc@meta.data
write.csv(meta,paste0("/output/QCed/tob_atac_", lib, "_annotated_meta_new_ref_nonPBMC_excluded.csv"))


#### END

