# Cell Type-Specific Epigenetic Age Analysis Pipeline

This directory contains scripts for analyzing epigenetic age patterns and genotype interactions in single-cell ATAC-seq data across different cell types.

**Author:** Peter C Allen  
**Project:** TenK10K Epigenetic Age and caQTL Analysis

## Overview

The pipeline performs cell type-specific analysis of epigenetic aging patterns and their interactions with genetic variants. It processes single-cell ATAC-seq data to:

1. Create cell type-specific data subsets
2. Compute epigenetic age estimates using EpiTrace
3. Test for genotype × epigenetic age interactions
4. Visualize results using UMAP and other plots

## Workflow

### 1. Cell Type Subsetting
**Script:** `1_celltype_subsetting.py`

Creates cell type-specific subsets from merged single-cell data with:
- Donor-stratified sampling for balanced representation
- Selection of highly variable peaks
- Export to MTX format for downstream analysis

```bash
python 1_celltype_subsetting.py --input merged_dataset.h5ad --output output/celltype_subsets
```

### 2. Epigenetic Age Computation
**Scripts:** `2_runEpiAge.R`, `2c_computeEpiage.sh`

Computes epigenetic age estimates for each cell type using EpiTrace:
- Processes MTX format data
- Applies EpiTrace algorithm with convergence criteria
- Outputs per-cell epigenetic age estimates

```bash
# Single cell type
Rscript 2_runEpiAge.R <celltype>

# All cell types via SGE array job
qsub 2c_computeEpiage.sh
```

### 3. Interaction Data Preparation
**Script:** `3_interactionCreation.py`

Processes epigenetic age estimates to create summary statistics:
- Computes median and standard deviation per donor
- Prepares data for interaction testing

```bash
python 3_interactionCreation.py --merged-file merged_dataset.h5ad
```

### 4. Genotype × Epigenetic Age Interaction Analysis
**Scripts:** `4_interactionRegression.R`, `4b_regressionJob.sh`

Tests for interactions between genetic variants and epigenetic age:
- Linear regression with covariates
- Tests SNP × epigenetic age interactions
- Multiple testing correction

```bash
# Single analysis
Rscript 4_interactionRegression.R <celltype> <chromosome>

# All combinations via SGE jobs
bash 4b_regressionJob.sh
```

### 5. Results Visualization
**Script:** `5-interactionPlotting.py`

Creates comprehensive visualizations of interaction results:
- Summary plots across cell types
- Significant interaction distributions
- Effect size comparisons

```bash
python 5-interactionPlotting.py
```

### 6. UMAP Visualizations
**Scripts:** `6-umap-epiagePeaks.py`, `6a-umap-epiageNK.py`

Generates UMAP plots colored by epigenetic age:
- Cell type-specific embeddings
- Genotype-stratified visualizations (NK cells)
- Multi-panel interaction plots

```bash
python 6-umap-epiagePeaks.py
python 6a-umap-epiageNK.py
```

## Input Data Requirements

- **Merged H5AD file:** Single-cell ATAC-seq data with cell type annotations
- **caQTL results:** Significant chromatin accessibility QTLs
- **Genotype data:** Individual-level SNP genotypes
- **Clinical data:** Disease status and demographic information
- **Covariate files:** Technical and biological covariates

## Output Structure

```
output/
├── celltype_subsets/           # Cell type-specific MTX files
├── regression-results/         # Interaction analysis results
└── merged_NK_cells.h5ad       # Processed NK cell data

figures/
├── 4-epigeneticAge-Celltype/  # Basic UMAP plots
├── 5-epigeneticAge-Genotype-UMAPs/  # Interaction UMAPs

data/
└── epiages/                   # Epigenetic age summaries per donor
```

## Dependencies

### R packages
- Seurat
- Signac
- EpiTrace
- lme4, lmerTest
- data.table, dplyr
- tibble, purrr

### Python packages
- scanpy
- anndata
- pandas, numpy
- matplotlib, seaborn
- scipy

## Notes

- Scripts are designed for SGE cluster environments
- Memory requirements vary by cell type (100-500GB)
- Runtime varies from minutes to hours per cell type
- Update file paths in scripts to match your environment
- Ensure conda environment 'signac' is available
