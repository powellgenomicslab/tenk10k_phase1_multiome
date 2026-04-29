# SMR Inference Pipeline (caQTL & eQTL → GWAS)

This directory contains scripts and workflows for performing **Summary-based Mendelian Randomisation (SMR)** analyses within the TenK10K Phase 1 multiome project.

The SMR framework integrates **molecular QTL signals (e.g., eQTL / caQTL)** with **GWAS summary statistics** to identify putative causal relationships between genetic variants, molecular traits, and complex diseases.

---

## Overview

The SMR analysis aims to:

- Test whether the effect of a genetic variant on a complex trait is mediated through:
  - gene expression (eQTL)
  - chromatin accessibility (caQTL)
- Distinguish **pleiotropy / causality** from **linkage** using the HEIDI test
- Prioritise candidate genes and regulatory elements underlying GWAS loci

This module operates downstream of:

- QTL mapping (eQTL / caQTL)
- fine-mapping
- GWAS summary statistics harmonisation

---

## File Structure

### Data Preparation

**1_tensorcaQTL2MatrixQTL.R**
- Converts tensor-based caQTL results into Matrix eQTL-compatible format

**2_convert2besd.sh**
- Converts QTL summary statistics into BESD format required by SMR

**3_createEPI.R**
- Generates `.epi` file containing probe (gene/peak) information

**4_createESI.R**
- Generates `.esi` file containing SNP information

---

### caQTL → eQTL SMR Analyses

**5_run_SMR_caQTL_eQTL.sh**
- Performs SMR between caQTL and eQTL
- Identifies peak-to-gene regulatory relationships

---

### caQTL → GWAS Integration

**6_1_run_SMR_GWAS_diseasetraits.sh**
- SMR analysis between QTL and disease GWAS traits

**6_2_run_SMR_GWAS_bloodtraits.sh**
- SMR analysis between QTL and blood trait GWAS

---

### eQTL → GWAS Integration

**7_1_run_SMR_GWAS_eQTL_diseasetraits.sh**
- SMR analysis between eQTL and disease GWAS traits

**7_2_run_SMR_GWAS_eQTL_bloodtraits.sh**
- SMR analysis between eQTL and blood trait GWAS

---

## Workflow
---

## Inputs

- caQTL summary statistics (tensor format)
- eQTL summary statistics
- GWAS summary statistics (disease and blood traits)
- LD reference panel

---

## Outputs

- SMR association results (effect sizes, p-values)
- HEIDI test results (heterogeneity p-values)
- Prioritised:
  - peak–gene links
  - gene–trait associations
  - regulatory mechanisms

---

## Key Dependencies

- SMR software (e.g., GCTA-SMR)
- R
- Bash

---

## Interpretation

Results from different steps can be combined to support causal chains such as:

SNP → chromatin accessibility → gene expression → disease

---

## References

- Zhu et al., 2016 — SMR method
- TenK10K Phase 1 multiome project

---

## 🚀 Usage

Scripts can be run individually or as part of a larger pipeline (e.g., HPC batch jobs). It is recommended to execute them sequentially following the workflow above.
