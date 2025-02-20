# TenK10K phase1 multiome data processing

TenK10K is a population cohort that will profile single-cell RNA-seq data of ~50 million peripheral blood mononuclear cells (PBMCs) from 10,000 individuals in Australia. This cohort will collect
- 5,000 healthy donors, an extended version of the [OneK1K cohort](https://www.science.org/doi/10.1126/science.abf3041)
- 2,000 individuals with coronary artery disease, the [BioHeart cohort](https://bmjopen.bmj.com/content/9/9/e028649)
- 1,500 samples from a pan-cancer cohort, the LBIO study
- 1,000 individuals from a pan-autoimmune disease cohort, the [AIM study](https://bmjopen.bmj.com/content/11/2/e042493.long)
- 600 samples from a long-covid cohort, the [ADAPT study](https://www.svhs.org.au/research-education/participating-in-research-trials/adapt-study)

<!-- potentially also from [HOPE Research Program](https://www.garvan.org.au/research/collaboration/hope-research) -->
This repository will contain the analysis code for computational and statistical analyses of TenK10K phase 1 ATAC/multiome data, with a focus on caQTL (i.e., SNPs associated with chromatin accessibility levels) and multi-omics integration.

We will include 952 TOB individuals for scATAC-seq, 48 BioHear Individuals and 27 LIBIO individuals for muliome data.

<br>

![Tenk10K](https://github.com/powellgenomicslab/tenk10k_phase1_multiome/blob/main/Figures/TenK10K_icon_poster.png)

# Main data for this project

We obtained three sets of scATAC/Multiome data from TOB, BioHeart and LBIO cohorts.

**TOB**: 119 pools (238 libraries) and 8 individuals per pool (generated in 2024)

**BioHeart**: 4 pools and 16 individuals per pool (generated back in 2022)

**LIBIO**: 4 pools and 6-8 individuals per pool (generated back in 2023)

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

## Gene regulatory network inference

- To be updated

# Data generation tracking (scATAC-seq for TOB cohort)
1st batch: Received on Aug 12nd 2024 (`240807`). 8 libraries from 4 pools.

2nd batch: Received on Sep 05th 2024 (`240905`). 8 libraries from 4 pools.

3rd batch: Received on Oct 02nd 2024 (`241002`). 20 libraries, from 11 pools.

4th batch: Received on Oct 08th 2024 (`241008`). 11 libraries, from 7 pools.

5th batch: Received on Oct 10th 2024 (`241010`). 10 libraries, from 6 pools.

6th batch: Received on Oct 15th 2024 (`241015`). 28 libraries, from 13 pools.

7th batch: Received on Oct 18th 2024 (`241018`). 6 libraries from 4 pools.

8th batch: Received on Nov 01st 2024 (`241101`). 7 libraries, from 5 pools.

9th batch: Received on Nov 04th 2024 (`241104`). 7 libraries from 4 pools.

10th batch: Received on Nov 06th 2024 (`241106`). 14 libraries, from 8 pools.

11st batch: Received on Nov 08nd 2024 (`241108`). 21 libraries, from 11 pools.

12nd batch: Received on Nov 11st 2024 (`241111`). 21 libraries, from 11 pools.

13rd batch: Received on Dec 02nd 2024 (`241202`). 31 libraries, from 17 pools.

14th batch: Received on Dec 04th 2024 (`241204`). 30 libraries, from 15 pools.

15th batch: Received on Dec 06th 2024 (`241206`). 16 libraries, from 8 pools.


In total, there are 238 scATAC-seq libraries from 119 pools

Note: each pool contains 8 donors and was captured/sequenced twice





