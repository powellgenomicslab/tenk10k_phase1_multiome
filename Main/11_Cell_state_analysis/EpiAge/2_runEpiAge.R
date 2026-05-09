#!/usr/bin/env Rscript

#' Epigenetic Age Estimation Using EpiTrace
#' 
#' This script computes epigenetic age estimates for single-cell ATAC-seq data
#' using the EpiTrace algorithm (Xiao, 2024). It processes cell type-specific data subsets
#' and generates epigenetic age predictions for individual cells.
#' 
#' Usage: Rscript 2_runEpiAge.R <celltype>
#' 
#' Author: Peter C Allen

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 2_runEpiAge.R <celltype>")
}
celltype <- args[1]

cat("Processing cell type:", celltype, "\n")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(reticulate)
  library(Matrix)
  library(EpiTrace)
})

# Configure Python environment (adjust path as needed)
tryCatch({
  reticulate::use_condaenv("signac", required = TRUE)
}, error = function(e) {
  cat("Warning: Could not activate conda environment, using system Python\n")
})

# Set up parallel processing
plan("multicore", workers = 6)
options(future.globals.maxSize = 200 * 1024^3)  # 200 GB memory limit

# Define file paths
base_dir <- "output/20250630_celltype_subsets"
celltype_dir <- file.path(base_dir, celltype)

cat("Loading data files...\n")

# Load count matrix
counts <- ReadMtx(
  mtx = file.path(celltype_dir, paste0(celltype, "_subset_matrix.mtx")),
  features = file.path(celltype_dir, paste0(celltype, "_subset_rep1_features.tsv")),
  cells = file.path(celltype_dir, paste0(celltype, "_subset_rep1_barcodes.tsv")),
  feature.column = 1,
  mtx.transpose = TRUE
)

# Load cell metadata
metadata <- read.csv(file.path(celltype_dir, paste0(celltype, "_subset_rep1_metadata.csv")), 
                     row.names = 1)

# Create Seurat object with chromatin accessibility data
cat("Creating Seurat object...\n")

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = NULL
)

seurat_obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

cat("Created Seurat object:", ncol(seurat_obj), "cells ×", nrow(seurat_obj), "peaks\n")

# Load EpiTrace functions
source("src/EpiTrace.R")

# Prepare data for EpiTrace analysis
cat("Preparing data for EpiTrace...\n")

peakset <- granges(seurat_obj@assays$peaks)
mtx <- as(seurat_obj@assays$peaks$counts, "sparseMatrix")

# Initialize EpiTrace components
initiated_peaks <- Init_Peakset(peakset)
initiated_peaks_df <- as.data.frame(initiated_peaks, row.names = NULL)

cellname_vec <- colnames(seurat_obj)
initiated_peaks_df$peakName <- paste0(initiated_peaks_df$seqnames, "_", 
                                     initiated_peaks_df$start, "_", 
                                     initiated_peaks_df$end)

initiated_mm <- Init_Matrix(cellname = cellname_vec, 
                           peakname = initiated_peaks_df$peakName, 
                           matrix = mtx)

# Free up memory
rm(counts, chrom_assay, metadata)
gc()


# Run EpiTrace epigenetic age estimation
cat("Running EpiTrace age estimation...\n")

# Helper function for EpiTrace computation
compute_epitrace <- function() {
  EpiTraceAge_Convergence(
    initiated_peaks,
    initiated_mm,
    celltype = NULL,
    qualnum = 10,
    Z_cutoff = 2.5,
    mean_error_limit = 0.01,
    iterative_time = 20,
    parallel = TRUE,
    ncore_lim = 3,
    ref_genome = "hg38",
    non_standard_clock = TRUE
  )
}

# Attempt EpiTrace computation with error recovery
tryCatch({
  epitrace_obj <- compute_epitrace()
  output_file <- file.path(base_dir, paste0(celltype, "_epitrace_age.csv"))
  write.csv(epitrace_obj@meta.data, file = output_file, row.names = TRUE)
  cat("EpiTrace completed successfully\n")
  cat("Results saved to:", output_file, "\n")
  
}, error = function(e) {
  cat("EpiTrace failed, attempting recovery by removing first peak...\n")
  cat("Error:", e$message, "\n")
  
  # Remove first peak and retry
  cat("Original object:", ncol(seurat_obj), "cells ×", nrow(seurat_obj), "peaks\n")
  seurat_obj <<- seurat_obj[-1, ]
  cat("Modified object:", ncol(seurat_obj), "cells ×", nrow(seurat_obj), "peaks\n")

  # Recompute EpiTrace components
  peakset <- granges(seurat_obj@assays$peaks)
  mtx <- as(seurat_obj@assays$peaks$counts, "sparseMatrix")

  initiated_peaks <<- Init_Peakset(peakset)
  initiated_peaks_df <- as.data.frame(initiated_peaks, row.names = NULL)

  cellname_vec <- colnames(seurat_obj)
  initiated_peaks_df$peakName <- paste0(initiated_peaks_df$seqnames, "_", 
                                       initiated_peaks_df$start, "_", 
                                       initiated_peaks_df$end)

  initiated_mm <<- Init_Matrix(cellname = cellname_vec, 
                              peakname = initiated_peaks_df$peakName, 
                              matrix = mtx)

  # Retry EpiTrace computation
  epitrace_obj <- compute_epitrace()
  output_file <- file.path(base_dir, paste0(celltype, "_epitrace_age.csv"))
  write.csv(epitrace_obj@meta.data, file = output_file, row.names = TRUE)
  cat("EpiTrace completed after retry\n")
  cat("Results saved to:", output_file, "\n")
})

cat("Epigenetic age computation completed for", celltype, "\n")
