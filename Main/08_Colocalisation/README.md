# 08_Colocalisation

This folder contains scripts for performing colocalisation (coloc) analysis between caQTL and GWAS/eQTL traits in the TenK10K phase 1 multiome dataset.

## Overview

The workflow consists of two main colocalisation analyses:

1. **caQTL to GWAS colocalisation** (`caQTL2GWAS/`)
   - Tests for shared causal variants between caQTL peaks and GWAS signals
   - Organized by trait category (BloodTraits, DiseaseTraits, NewDiseases)
   - Workflow: Process GWAS summary statistics → Run colocalisation

2. **caQTL to eQTL colocalisation** (`caQTL2eQTL/`)
   - Tests for shared causal variants between caQTL peaks and eQTL signals
   - Workflow: Run colocalisation on preprocessed eQTL data

## Files and Structure

### caQTL2GWAS/

This directory contains subdirectories for different trait categories, each with the same workflow:

#### BloodTraits/, DiseaseTraits/, NewDiseases/

- `1_process_GWAS.R` / `1_process_GWAS.sh`  
  Processes and formats GWAS summary statistics for colocalisation (extract relevant columns, filter, standardize format).

- `2_coloc.R` / `2_coloc.sh`  
  Runs colocalisation analysis between caQTL peaks and processed GWAS data using the `coloc` R package.

### caQTL2eQTL/

- `1_coloc.R` / `1_coloc.sh`  
  Runs colocalisation analysis between caQTL peaks and eQTL signals.

## Output

The outputs include colocalisation results (posterior probabilities of colocalisation, hypothesis support) indicating shared causal variants between caQTL and GWAS/eQTL loci.

