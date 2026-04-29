# Causal Inference Pipeline

This directory contains workflows for **causal inference analysis** in the TenK10K Phase 1 multiome project, integrating molecular QTL data with GWAS summary statistics.

The pipeline is organised into two main modules:

- **SMR-based inference** for testing causal relationships between molecular traits and complex traits  
- **GWAS clumping and colocalisation** for identifying shared genetic signals between QTLs and GWAS loci  

---

## Directory Structure

---

## Modules

### 1. SMR_inference

Implements **Summary-based Mendelian Randomisation (SMR)** analyses to:

- Link genetic variants to:
  - chromatin accessibility (caQTL)
  - gene expression (eQTL)
  - complex traits (GWAS)
- Identify putative **causal relationships**
- Apply HEIDI tests to distinguish causality from linkage

This module supports analyses such as:

- caQTL ↔ eQTL (regulatory inference)  
- caQTL ↔ GWAS (chromatin-mediated effects)  
- eQTL ↔ GWAS (expression-mediated effects)  

See `SMR_inference/README.md` for details.

---

### 2. GWAS_clumping_colocalization

Performs **GWAS signal refinement and colocalisation analysis** to:

- Identify independent GWAS loci using PLINK clumping  
- Evaluate colocalisation between QTL and GWAS signals  
- Extract and summarise high-confidence colocalised signals (e.g., PP.H4 and p_SMR)  
- Visualise results across traits  

This module includes:

- Clumping of GWAS summary statistics  
- Colocalisation result comparison and aggregation  
- Visualisation of colocalisation signals  
---

## Analysis Strategy

These two modules are complementary:

- **Colocalisation analysis** prioritises loci with shared genetic signals  
- **SMR analysis** tests causal relationships at those loci  

Together, they provide a robust framework for identifying:

> Genetic variants → molecular traits → complex traits

---

## Requirements

- R  
- Bash  
- PLINK  
- SMR software (GCTA-SMR)  

---
