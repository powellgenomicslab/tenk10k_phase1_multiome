# Multi-Trait Colocalization Analysis

This folder contains scripts for performing multi-trait colocalization analysis between caQTL and multiple GWAS/eQTL traits in the TenK10K phase 1 multiome dataset.

## Overview

The workflow combines GWAS and eQTL data with caQTL results to identify shared causal signals across multiple traits simultaneously using multi-trait colocalization.

Workflow:
1. Process and format GWAS summary statistics
2. Generate sample sizes and standardize data across traits
3. Run multi-trait colocalization analysis

## Files

- `1_process_GWAS.R`  
  Processes and formats GWAS summary statistics for colocalization (extract relevant columns, filter, standardize format).

- `2_generate_N.R`  
  Generates and standardizes sample sizes for GWAS and other traits used in multi-trait colocalization.

- `3_coloc.R` / `3_coloc.sh`  
  Runs multi-trait colocalization analysis combining caQTL peaks with processed GWAS/eQTL data using multi-trait coloc methods.

## Output

The outputs include multi-trait colocalization results indicating shared causal variants between caQTL and multiple GWAS/eQTL loci simultaneously.