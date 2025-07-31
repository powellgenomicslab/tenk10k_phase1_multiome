args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide a celltype as an argument.")
}
celltype <- args[1]

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(reticulate)
library(Matrix)
reticulate::use_condaenv("/home/petall/miniforge3/envs/signac", required = T)

plan("multicore", workers = 6)
options(future.globals.maxSize = 200000 * 1024^2) # for 200 Gb RAM

# Read matrix
counts <- ReadMtx(
  mtx = file.path(paste0("output/20250630_celltype_subsets/", celltype, "/", celltype, "_subset_matrix.mtx")),
  features = file.path(paste0("output/20250630_celltype_subsets/", celltype, "/", celltype, "_subset_rep1_features.tsv")),
  cells = file.path(paste0("output/20250630_celltype_subsets/", celltype, "/", celltype, "_subset_rep1_barcodes.tsv")),
  feature.column = 1,
  mtx.transpose = TRUE
)

# Optionally read metadata
metadata <- read.csv(file.path(paste0("output/20250630_celltype_subsets/", celltype, "/", celltype, "_subset_rep1_metadata.csv")), row.names = 1)

# Create ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),  # adjust if your peak format is different
  genome = 'hg38',    # or 'mm10', etc.
  fragments = NULL    # you can add fragment files later
)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

library(EpiTrace)
source("src/EpiTrace.R")

peakset <- granges(seurat_obj@assays$peaks)
mtx <- as(seurat_obj@assays$peaks$counts, "sparseMatrix")

initiated_peaks <- Init_Peakset(peakset)
initiated_peaks_df <- as.data.frame(initiated_peaks, row.names = NULL)

colnames(seurat_obj) -> cellname_vec
paste0(initiated_peaks_df$seqnames, "_", initiated_peaks_df$start, "_", initiated_peaks_df$end) -> initiated_peaks_df$peakName

initiated_mm <- Init_Matrix(cellname = cellname_vec, peakname = initiated_peaks_df$peakName, matrix = mtx)
rm(counts, chrom_assay, metadata)
gc()


tryCatch({
  print(seurat_obj)
  epitrace_obj_age_estimated <- EpiTraceAge_Convergence(
    initiated_peaks,
    initiated_mm,
    celltype = NULL,
    qualnum = 10,
    Z_cutoff = 2.5,
    mean_error_limit = 0.01,
    iterative_time = 20,
    parallel = T,
    ncore_lim = 3,
    ref_genome = "hg38",
    non_standard_clock = T
  )
  write.csv(epitrace_obj_age_estimated@meta.data, file = paste0("output/20250630_celltype_subsets/", celltype, "_epitrace_age.csv"), row.names = T)
  message("EpiTraceAge_Convergence completed successfully.")
}, error = function(e) {
  message("EpiTraceAge_Convergence failed. Restarting by removing first feature.")
  message("Previous Seurat Object:")
  print(seurat_obj)
  seurat_obj <- seurat_obj[-1,]
  message("New Seurat Object:")
  print(seurat_obj)

  peakset <- granges(seurat_obj@assays$peaks)
  mtx <- as(seurat_obj@assays$peaks$counts, "sparseMatrix")

  initiated_peaks <- Init_Peakset(peakset)
  initiated_peaks_df <- as.data.frame(initiated_peaks, row.names = NULL)

  colnames(seurat_obj) -> cellname_vec
  paste0(initiated_peaks_df$seqnames, "_", initiated_peaks_df$start, "_", initiated_peaks_df$end) -> initiated_peaks_df$peakName

  initiated_mm <- Init_Matrix(cellname = cellname_vec, peakname = initiated_peaks_df$peakName, matrix = mtx)

  epitrace_obj_age_estimated <- EpiTraceAge_Convergence(
    initiated_peaks,
    initiated_mm,
    celltype = NULL,
    qualnum = 10,
    Z_cutoff = 2.5,
    mean_error_limit = 0.01,
    iterative_time = 20,
    parallel = T,
    ncore_lim = 3,
    ref_genome = "hg38",
    non_standard_clock = T
  )
  write.csv(epitrace_obj_age_estimated@meta.data, file = paste0("output/20250630_celltype_subsets/", celltype, "_epitrace_age.csv"), row.names = T)
  message("EpiTraceAge_Convergence completed successfully after retry.")
})

message("Job Complete.")
