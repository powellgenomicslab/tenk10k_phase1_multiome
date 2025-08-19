# Cell Type-Specific Epigenetic Age Analysis Pipeline

This directory contains the complete analysis pipeline for studying epigenetic aging patterns and genetic interactions in single-cell ATAC-seq data across different cell types. The scripts listed below were used for the publication analysis.

**Author:** Peter C Allen  
**Project:** TenK10K Epigenetic Age and caQTL Analysis

## Analysis Scripts Overview

This pipeline consists of 12 core scripts that process single-cell ATAC-seq data to identify and visualize genetic variants that interact with epigenetic aging patterns:

### Data Processing Pipeline (Steps 1-6)
1. **`1_celltype_subsetting.py`** - Cell type-specific data subset creation
2. **`2_runEpiAge.R`** - Epigenetic age computation using EpiTrace
3. **`2c_computeEpiage.sh`** - SGE job submission for epigenetic age computation
4. **`3_interactionCreation.py`** - Interaction analysis data preparation
5. **`4_interactionRegression.R`** - Statistical testing of SNP × epigenetic age interactions
6. **`4b_regressionJob.sh`** - SGE job submission for interaction regression

### Visualization and Analysis (Steps 7-12)
7. **`plotAllCelltypes.R`** - Comprehensive cell type visualization
8. **`5-interactionPlotting.py`** - Comprehensive interaction results visualization
9. **`6-umap-epiagePeaks.py`** - UMAP analysis with epigenetic age overlays
10. **`6a-plotNK-interaction.py`** - NK cell-specific interaction plotting

#### Supporting Scripts
11. **`functions_5.py`** - Modular analysis functions library
12. **`nk_umap_plotting.py`** - Modular functions to plot epigenetic age + chromatin accessibility x Genotype|EpiAge


## Detailed Workflow

### Step 1: Cell Type Subsetting
**Script:** `1_celltype_subsetting.py`

Creates cell type-specific subsets from merged single-cell data:
- Donor-stratified sampling for balanced representation
- Selection of highly variable peaks
- Export to MTX format for downstream analysis

```bash
python 1_celltype_subsetting.py --input merged_dataset.h5ad --output output/celltype_subsets
```

### Step 2: Epigenetic Age Computation
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

### Step 3: Interaction Data Preparation
**Script:** `3_interactionCreation.py`

Processes epigenetic age estimates to create summary statistics:
- Computes median and standard deviation per donor
- Prepares data for interaction testing
- Creates binary age classifications (Young/Old)

```bash
python 3_interactionCreation.py --merged-file merged_dataset.h5ad
```

### Step 4: Statistical Interaction Testing
**Scripts:** `4_interactionRegression.R`, `4b_regressionJob.sh`

Tests for interactions between genetic variants and epigenetic age:
- Linear regression with technical and biological covariates
- Tests SNP × epigenetic age interactions on chromatin accessibility
- Applies multiple testing correction across variants

```bash
# Single analysis
Rscript 4_interactionRegression.R <celltype> <chromosome>

# All combinations via SGE jobs
bash 4b_regressionJob.sh
```

### Step 5: Results Visualization
**Script:** `5-interactionPlotting.py`
**Support Script:** `functions_5.py`

Creates comprehensive visualizations of interaction results:
- Summary plots across cell types showing significant interactions
- Effect size distributions and directional analysis
- Statistical comparisons between cell types

```bash
python 5-interactionPlotting.py
```

### Step 6: UMAP Analysis
**Scripts:** `6-umap-epiagePeaks.py`, `6a-plotNK-interaction.py`
**Support Script:** `nk_umap_plotting.py`

Generates UMAP plots with epigenetic age and genotype overlays:
- Cell type-specific embeddings colored by epigenetic age
- NK cell-specific genotype-stratified visualizations
- Multi-panel interaction effect plots

```bash
python 6-umap-epiagePeaks.py
python 6a-plotNK-interaction.py
```

### Step 7: Comprehensive Visualization
**Script:** `plotAllCelltypes.R`

Creates plot of the epigenetic ages of all cell-types -- EpiTrace Age calculated for each pool then combined

```bash
Rscript plotAllCelltypes.R
```

## Results

The analysis pipeline identified:
- **174,045 total tests** performed across 24 cell types
- **3,080 significant interactions** (FDR < 0.05)
- **Directional bias:** 54.2% positive vs 45.8% negative effects (p = 1.8×10⁻⁶)
- **Multiple GWAS hits** with significant genotype × epigenetic age interactions

## Data Structure

### Input Requirements
```
data/
├── ExpressionBeds/chr{N}.bed.gz  # Pseudobulked Accessibility data
├── TenK10K_*_chr{N}_common_variants.raw      # Genotype data
├── {celltype}_epiage.txt                 # Epigenetic ages
└── covariates.txt                                         # Sample metadata
```

### Output Structure
```
output/
├── celltype_subsets/           # Cell type-specific MTX files
├── regression-results/         # Interaction analysis results
└── merged_NK_cells.h5ad       # Processed NK cell data

figures/
├── 3-regression/              # Interaction analysis plots
├── 4-epigeneticAge-Celltype/  # Basic UMAP plots
└── 5-epigeneticAge-Genotype-UMAPs/  # Interaction UMAPs
```

## Usage Notes

- Scripts are optimized for SGE cluster environments
- Update file paths in scripts to match your data structure
- Ensure conda environment with required packages is available
- Some scripts require substantial memory (up to 32GB for large cell types)
- Pipeline can be run in parallel for different cell types

