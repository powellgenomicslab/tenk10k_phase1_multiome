# TenK10K phase1 multiome data processing

TenK10K is a population cohort that will profile single-cell RNA-seq data of ~50 million peripheral blood mononuclear cells (PBMCs) from 10,000 individuals in Australia. This cohort will collect
- 5,000 healthy donors, an extended version of the [OneK1K cohort](https://www.science.org/doi/10.1126/science.abf3041)
- 2,000 individuals with coronary artery disease, the [BioHeart cohort](https://www.heartresearch.com.au/research/bioheart/)
- 1,500 individuals from a pan-autoimmune disease cohort, the [AIM study](https://www.mrc.unsw.edu.au/australian-inflammatory-bowel-disease-microbiome-study)
- 1,500 individuals from a pan-cancer cohort, the LBIO study
- 600 individuals from a long-covid cohort, the [ADAPT study](https://www.svhs.org.au/research-education/participating-in-research-trials/adapt-study)

<!-- potentially also from [HOPE Research Program](https://www.garvan.org.au/research/collaboration/hope-research) -->
This repository will contain the analysis code for computational and statistical analyses of TenK10K phase 1 ATAC/multiome data, with a focus on caQTL (i.e., SNPs associated with chromatin accessibility levels) and multi-omics integration.
We will include 1,000 TOB individuals, 900 for ATAC-seq and 100 for muliome data.

<br>

![Tenk10K](https://github.com/powellgenomicslab/tenk10k_phase1_multiome/blob/main/Figures/TenK10K_icon_poster.png)

# Pilot data processing

We obtained two sets of pilot data from BioHeart and LBIO cohorts.

**BioHeart**: 4 pools and 16 individuals per pool (generated back in 2022)

**LIBIO**: 4 pools and 8 individuals per pool (currently undergoing)

## Preprocessing

- Reads alignment (cell ranger ARC)
- Ambient RNA detection
- Demultiplexing and doublet detection
- Cell type annotation
- Batch correction
- Multi-omics layers integration
- QC & Normalisation

## ATAC-seq processing

- Peak calling by MACS2
- DNA accessibility process by latent semantic indexing (LSI)
- Create a gene activity matrix
- Clustering using Azimuth / label transferring
- Verify with scRNA-seq data
- Generate Peak Count matrix

## caQTL mapping

- To be updated

# Data generation tracking (ATAC-seq)
1st batch: Received on Aug 12nd 2024 (`240807`). 8 libraries, 4 pools, and 32 TOB individuals (each pool repeated twice).







