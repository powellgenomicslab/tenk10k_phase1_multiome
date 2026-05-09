# caQTL mapping using SAIGEQTL

This folder contains the scripts and helper functions to perform caQTL mapping using SAIGE-based association tests. The workflow is split into three main parts:

- `Phenotypes/` — generate and prepare phenotype (pseudobulk) matrices and peak lists. **Run this first.**
- `VRE/` — prepare variant-region extraction inputs used for downstream association testing. **Run this after `Phenotypes/`, but before `Unconditional/` or `Conditional/`.**
- `Unconditional/` — run unconditional (single-variant) association tests and combine results.
- `Conditional/` — run conditional association tests (if conditional analyses are required) and combine results.

There is also a helper `qval_functions.R` used to compute q-values/FDR across results.

## Overview of key scripts:

`Phenotypes/`

- `0_filter_peaks.py`
  - Filter raw peaks and apply initial peak-level filters prior to matrix creation.
- `1_process_objects.py`
  - Process pseudobulk objects and metadata; prepares donor/sample-level inputs.
- `1.5_combine_metadata.R`
  - Helper to merge metadata tables used across steps.
- `2_make_matrices.py` 
  - Build per-donor (pseudobulk) count matrices ready for association testing.
- `3_combine_files.py` 
  - Combine per-chromosome or per-chunk files into unified phenotype files.
- `4_get_peak_lists.R`
  - Generate final peak lists/bed files used to define phenotype regions.

`VRE/`

- `1_get_VRE_plink.sh`
  - Generate PLINK files for variant-region extraction workflows.
- `2_combine_variants.sh`
  - Combine variant files for downstream use.
- `3_get_regions.R`
  - Build region lists used by later association testing steps.

`Unconditional/`

- `1.5_get_regions.R`
  - Generate region lists/bed files for genome chunks to test.
- `1_index_VCF.sh`
  - Index VCFs for fast access (bgzip/tabix or equivalent) before running SAIGE.
- `2_run_SAIGE.sh`
  - Chromosome/chunk-specific runner script that executes SAIGE for unconditional tests.
- `3_combine_raw_results.sh`
  - Combine raw per-chunk SAIGE outputs into a single file per phenotype/celltype.
- `4_combine_single_assoc.sh`
  - Post-process combined results into single-variant association tables.
- `5_compute_q_values.R`
  - Compute q-values/FDR and annotate significant hits.

`Conditional/`

- `0_convert_VCF.sh`
  - Prepare/convert VCFs for conditional analyses (subset/format as needed).
- `1_prepare_SNPs.R`
  - Select and format SNP lists used as covariates in conditional models.
- `2_run_SAIGE_cond.sh`
  - Runner script to execute conditional SAIGE analyses given conditioning SNPs.
- `3_combine_raw_results.sh`
  - Combine raw conditional outputs across chunks.
- `4_combine_single_assoc.sh`
  - Post-process combined conditional results.
- `5_compute_q_values.R`
  - Compute q-values/FDR for conditional results.

## Recommended running order:

1. Run the entire `Phenotypes/` pipeline to produce phenotype matrices and peak lists. This step is required before any association testing can proceed.

2. Run `VRE/` to prepare the variant-region extraction inputs needed for downstream analysis.

3. Run either `Unconditional/` or `Conditional/` depending on your analysis plan:
   - For a standard single-variant scan, run `Unconditional/`.
   - If you need conditional analyses, run `Conditional/` after producing the SNP lists to condition on