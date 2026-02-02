# Quality control of scATAC-seq datasets

Nuclei that failed the following quality control thresholds were excluded from analyses: number of peak counts > 9,000 and < 100,000, percentage of reads overlapping with peaks > 40, blacklist ratio < 0.01, nucleosome signal < 4, TSS enrichment score > 4. Seven sequencing libraries with low number of nuclei after quality control (< 300) were excluded, which resulted in the loss of one pool (8 donors excluded). After matching with whole genome sequencing data, 20 donors could not be matched and were removed. One pool (containing 8 donors) was excluded due to the extremely low number of nuclei. The final data includes 3,472,552 nuclei from 922 TOB donors. On average, each donor had 3,766 nuclei from scATAC-seq and 2,825 cells from matching scRNA-seq data.

`atac_01_atac_to_seurat.R`

This script converts 10X Genomics ATAC-seq outputs for all TOB cohort libraries into Seurat objects, integrating fragment and metadata files while adding demultiplexing results from Vireo. It then calculates standard QC metrics including nucleosome signal, TSS enrichment, fraction of reads in peaks, and blacklist ratios, and saves both the Seurat object and its metadata for downstream analysis.


`atac_02_qc.R`

This script performs regular QC on previously generated ATAC-seq Seurat objects by removing doublets and low-quality cells based on metrics such as peak counts, TSS enrichment, nucleosome signal, and fraction of reads in peaks, and saves the filtered objects. It then applies normalization, dimensionality reduction, clustering, UMAP visualization, and generates a gene activity matrix to enable RNA-based feature analysis and downstream interpretation.

The `*.qsub.sh` scripts are bash scripts to submit array jobs in HPC.
