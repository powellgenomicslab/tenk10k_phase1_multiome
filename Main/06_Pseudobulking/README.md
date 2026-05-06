# 06_Pseudobulking

This folder contains scripts for generating pseudobulk ATAC-seq count matrices from single-cell multiome data.

## Overview

The main script aggregates ATAC peak counts by **donor** within a selected **cell type**. For each run, it:

1. reads a cell type name from `celltype_names.txt` using a command-line index,
2. loads input `.h5ad` files,
3. filters cells by `predicted.id`,
4. aggregates counts by `donor_id`,
5. concatenates results across files, and
6. saves the final pseudobulk matrix as a CSV file in `PseudobulkMatrices/`.

## Main script

- `1_pseudobulking.py` — generate donor-level pseudobulk matrices for a specified cell type.

## Inputs

- `celltype_names.txt`
- ATAC `.h5ad` files from the configured matrix directory

## Output

- `PseudobulkMatrices/<celltype>.csv`

## Notes

- Only files ending in `_1.h5ad` are processed.
- Cell types are defined using the `predicted.id` metadata field.
- Aggregation is performed using `donor_id`.
