# ATAC-seq Quality Control Pipeline

This directory contains scripts for processing and quality control of ATAC-seq data from the TenK10K phase 1 multiome study.

## Scripts Overview

### atac_01_atac_to_seurat.R

**Purpose:** Converts Cell Ranger ATAC output to Seurat objects and performs initial quality control for TOB cohort libraries.

**What it does:**
1. **Data Conversion:**
   - Reads Cell Ranger ATAC outputs (filtered peak-barcode matrix, single-cell metadata, and fragment files)
   - Creates Chromatin Assay with quality filters (minimum 10 cells per peak, minimum 200 features per cell)
   - Builds a Seurat object with the "peaks" assay

2. **Demultiplexing Integration:**
   - Integrates Vireo demultiplexing results to assign cells to individual donors
   - Adds donor assignment information to the Seurat object metadata

3. **Quality Control Metrics:**
   - Filters to standard chromosomes only (removes non-standard contigs)
   - Adds gene annotations from Ensembl (hg38 genome build)
   - Computes nucleosome signal score per cell
   - Calculates TSS (Transcription Start Site) enrichment score
   - Computes percentage of reads in peaks
   - Calculates blacklist ratio (proportion of reads in ENCODE blacklisted regions)
   - Classifies cells by nucleosome signal (NS > 4 vs NS < 4)

4. **Output:**
   - Saves processed Seurat object as `tob_atac_[pool]_raw.RDS`
   - Exports metadata as `tob_atac_[pool]_meta.csv`

**Usage:**
```bash
Rscript atac_01_atac_to_seurat.R [pool_name]
```

**Input:**
- `pool_name`: Name of the pool to process (e.g., "S0228_2")
- Cell Ranger ATAC output directory containing:
  - `filtered_peak_bc_matrix.h5`
  - `singlecell.csv`
  - `fragments.tsv.gz`
- Vireo demultiplexing results from `/g/data/ei56/ax3061/proj/tenk10k/caQTL/demultiplexing/output/Vireo/force_cell/[pool]/donor_ids.tsv`

**Output:**
- `./data/raw/tob_atac_[pool]_raw.RDS` - Seurat object with QC metrics
- `./data/raw/tob_atac_[pool]_meta.csv` - Cell metadata CSV file

**Dependencies:**
- Seurat
- Signac
- EnsDb.Hsapiens.v86
- BSgenome.Hsapiens.UCSC.hg38
- ggplot2
- dplyr
- glue
- purrr

### atac_02_qc.R

Additional quality control and filtering steps for ATAC-seq data.

## Job Submission

The directory includes helper scripts for batch job submission:
- `create_atac_01_atac_to_seurat_jobs.sh` - Generate job scripts for all pools
- `atac_01_atac_to_seurat_pool_*.qsub.sh` - Example job submission scripts

## Notes

- Each pool contains 8 donors and was sequenced twice (2 technical replicates)
- The script handles flexible naming patterns for pool directories (e.g., "S0228_2", "S0228_R_2", "S0228_2_R")
- Standard chromosomes filtering removes mitochondrial DNA and non-standard contigs
- Nucleosome signal > 4 typically indicates lower quality cells
