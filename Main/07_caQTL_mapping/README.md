# caQTL mapping using TensorQTL and SAIGE-QTL

This folder contains the scripts for caQTL mapping with two software packages

## TensorQTL
The caQTL results from TensorQTL were used as the primary results in our paper.

`1_correct_GC.R`

This script processes pseudobulk ATAC-seq matrices for a specified cell type by reading the raw counts, filtering peaks, creating genomic ranges, calculating GC content, and normalizing read counts based on GC bias. It then saves the resulting normalized matrix to a dedicated output directory for downstream analyses.

`2_preprocessing.R`

This script processes normalized pseudobulk ATAC-seq matrices from two technical repeats for a given cell type by filtering low-count peaks, averaging overlapping donors, standardizing counts, and generating per-chromosome expression BED files. It then performs PCA to select top components, combines them with donor covariates, and saves the final standardized matrix, expression files, and covariate table for downstream QTL analyses.

`3_tensorQTL.py`

This script performs chromosome-specific cis-caQTL mapping for a given cell type by loading genotype and phenotype data, filtering donors with missing values, and computing nominal and permutation-based association p-values using TensorQTL. It then calculates FDR-adjusted q-values, identifies significant and conditionally independent QTLs, and saves all results to structured output files for downstream analysis.

<br>

## SAIGE-QTL

We also performed SAIGE-QTL analysis for caQTL mapping. 





The caQTL summary statistics are publicly available and can be downloaded at Zenodo (#### link to be updated)



