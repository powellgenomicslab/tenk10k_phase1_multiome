# 10_Fine_Mapping

This folder contains scripts for fine-mapping caQTL signals in the TenK10K phase 1 multiome dataset using SuSiE and mvSuSiE approaches.

## Overview

The workflow consists of two fine-mapping approaches:

1. **SuSiE fine-mapping** (`SuSiE/`)
   - Fine-mapping of caQTL signals per cell type using univariate SuSiE
   - Workflow: Process caQTL results → Run SuSiE fine-mapping

2. **mvSuSiE fine-mapping** (`mvSuSiE/`)
   - Multivariate fine-mapping combining caQTL with colocalized GWAS/eQTL signals
   - Aggregates fine-mapped credible sets across multiple traits

## Files and Structure

### SuSiE/

- `1_process_caQTL.R` / `1_process_caQTL.sh`  
  Processes and formats caQTL summary statistics (extract effect sizes, standard errors, and supporting data) for fine-mapping input.

- `2_finemapping.R` / `2_finemapping.sh`  
  Runs univariate SuSiE fine-mapping on processed caQTL data to identify credible sets of causal variants per locus.

### mvSuSiE/

- `collect_genes_for_ld.R`  
  Collects gene/peak information and linkage disequilibrium (LD) data for use in multivariate fine-mapping.

- `run_mvsusie_test_v9_siggwas_colocsmrQCed_1mb.R`  
  Runs multivariate SuSiE fine-mapping combining caQTL, colocalized GWAS, and QC'd signals within 1 Mb windows.

## Output

The outputs include credible sets of fine-mapped variants with posterior probabilities, indicating likely causal variants for each caQTL signal.