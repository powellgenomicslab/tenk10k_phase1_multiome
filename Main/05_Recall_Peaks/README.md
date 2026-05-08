# 05_Recall_Peaks

This folder contains scripts for recalling ATAC peaks across cell types in the TenK10K phase 1 multiome dataset.

## Overview

The workflow:
1. Generates lists of libraries and cell types from annotated metadata.
2. Extracts fragment files for each cell type within each library.
3. Combines fragments across libraries for each cell type.
4. Calls peaks on the combined fragments using MACS3 via `snapatac2`.
5. Merges recalled peaks across cell types into a final combined peak set.

## Files

- `0_generate_lists.R`  
  Creates library and cell type lists and initializes output directories.

- `1_celltype_fragments.R`  
  Extracts cell type-specific fragments from each library using metadata annotations.

- `2_combine_fragments.sh`  
  Concatenates per-library fragment files into one combined fragment file per cell type.

- `3_MACS3.py`  
  Runs MACS3 peak calling on combined cell type fragment files.

- `4_merge_peaks.R`  
  Merges peaks across cell types and exports a final BED file.

## Output

The main output of this workflow is a merged BED file of recalled peaks across cell types.
