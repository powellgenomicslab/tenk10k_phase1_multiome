# WASP Analysis

This folder contains scripts for allele-specific expression (ASE) analysis using the WASP pipeline in the TenK10K phase 1 multiome dataset.

## Overview

The workflow performs allele-specific expression calling and analysis by remapping reads to a personalized reference genome and quantifying allele-specific counts for downstream QTL mapping.

Workflow:
1. Generate and assign cell barcodes
2. Subset BAM files by cell type/library
3. Convert BAM to FASTQ and realign with WASP
4. Convert genotypes and BAM files
5. Prepare CHT input and adjust read counts
6. Calculate overdispersions and run statistical tests

## Files

- `1_generate_barcodes.R`  
  Generates and assigns cell barcodes from metadata for downstream sample tracking.

- `2_subset_bams.sh`  
  Subsets BAM files by cell type or library based on barcode assignments.

- `3_bam_to_fastq.sh`  
  Converts BAM files to FASTQ format for realignment.

- `3.5_prepare_STAR.sh`  
  Prepares STAR genome indices and reference files for WASP-aware alignment.

- `4_star_wasp.sh`  
  Runs STAR aligner with WASP mode enabled to map reads while tracking alleles.

- `5_convert_geno.sh`  
  Converts and formats genotype data (VCF format) for WASP processing.

- `6_convert_BAM.sh`  
  Converts aligned BAMs into standardized format for ASE analysis.

- `7_get_regions.sh`  
  Extracts genomic regions of interest for CHT input preparation.

- `8_make_CHT_input.sh`  
  Generates input files for CHT (comparative hypothesis testing or similar analysis tool).

- `9_adjust_read_counts.sh`  
  Adjusts allele-specific read counts for bias and other artifacts.

- `10_adjust_geno_errors.sh`  
  Adjusts genotyping errors in allele-specific counts.

- `11_overdispersions.sh`  
  Calculates overdispersion parameters for statistical modeling.

- `12_run_CHT.sh`  
  Runs statistical tests on allele-specific counts.

## Output

The outputs include allele-specific expression counts and statistical results for downstream allele-specific QTL mapping.