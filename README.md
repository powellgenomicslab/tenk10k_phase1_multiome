# TenK10K phase1 multiome project

TenK10K is a population cohort that will profile single-cell RNA-seq data of ~50 million peripheral blood mononuclear cells (PBMCs) from 10,000 individuals in Australia. This cohort will collect
- 5,000 healthy donors, an extended version of the [OneK1K cohort](https://www.science.org/doi/10.1126/science.abf3041)
- 2,000 individuals undergoing clinically indicated CT coronary angiogram, the [BioHeart cohort](https://bmjopen.bmj.com/content/9/9/e028649)
- 1,500 samples from a pan-cancer cohort, the LBIO study (500 donors, 1~5 timepoints for each)
- 1,000 individuals from a pan-autoimmune disease cohort, the [AIM study](https://bmjopen.bmj.com/content/11/2/e042493.long)
- 600 samples from a long-COVID cohort, the [ADAPT study](https://www.svhs.org.au/research-education/participating-in-research-trials/adapt-study)

<!-- potentially also from [HOPE Research Program](https://www.garvan.org.au/research/collaboration/hope-research) -->
This repository will contain the analysis code for computational and statistical analyses of TenK10K phase 1 ATAC/multiome data, with a focus on caQTL (i.e., SNPs associated with chromatin accessibility levels) and multi-omics integration.

In the phase 1 multiome project, we processed 952 TOB individuals for scATAC-seq, 64 BioHeart individuals, and 26 LIBIO individuals for muliome data.

<br>

![Tenk10K](https://github.com/powellgenomicslab/tenk10k_phase1_multiome/blob/main/Figures/TenK10K%20-%20Phase%201%20cohort.png)

# Main data for this project

We obtained PBMCs of 1042 individuals from three sets of scATAC/Multiome data from TOB, BioHeart and LBIO cohorts.

**TOB**: 952 donors, 119 pools (238 libraries) and 8 individuals per pool (generated in 2024)

**BioHeart**: 64 donors, 4 pools and 16 individuals per pool (generated in 2022)

**LBIO**: 26 donors, 4 pools and 6-8 individuals per pool (generated in 2023)

After quality controls, 922 TOB donors were used for caQTL mapping and 60 donors (BioHeart + LBIO) were used for replication.

# Data availability
We are currently preparing the data sharing, including full caQTL summary statistics (both common and rare variants), fine-mapping, coloc, SMR, peak-gene links, and GRN inference results.
The download link will be posted here in May 2026.

# Citation
Xue et al. Genetic regulation of cell type–specific chromatin accessibility shapes immune function and disease risk. [medRxiv](https://www.medrxiv.org/content/10.1101/2025.08.27.25334533v1). 2025.

# Contact
Please send all enquiries regarding this study to the corresponding authors:
Angli Xue (a.xue@garvan.org.au) and Joseph Powell (j.powell@garvan.org.au)


