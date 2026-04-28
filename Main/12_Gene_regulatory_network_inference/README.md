# Gene Regulatory Network (GRN) Inference Pipeline

This directory contains scripts for constructing **gene regulatory networks (GRNs)** from single-cell multiome data (scRNA-seq + scATAC-seq) using an integrative framework based on **scGLUE** and **SCENIC**.

The pipeline combines:

- Multi-omics integration (RNA + ATAC)
- Prior regulatory graph construction (distance, overlap, eQTL)
- Feature embedding via GLUE
- Gene–peak and TF–gene inference
- Cis-regulatory ranking and network pruning

---

## 🔬 Overview

The goal of this pipeline is to infer **cell-type-specific regulatory networks** by integrating:

- Chromatin accessibility (ATAC)
- Gene expression (RNA)
- Prior biological knowledge:
  - genomic distance
  - eQTL links
  - TF binding (ChIP-seq)

This enables reconstruction of regulatory relationships:

> TF → regulatory element (peak) → target gene

---

## 📁 File Structure
---

## ⚙️ Pipeline Steps

### 1. Preprocessing and Prior Construction  
**`s01_preprocessing.py`**

- Loads scRNA-seq and scATAC-seq data  
- Performs:
  - quality control and filtering  
  - normalization and dimensionality reduction  
  - batch correction (Harmony)  
- Annotates genomic coordinates  
- Constructs regulatory graphs:
  - peak–gene overlap  
  - genomic distance  
  - eQTL links  
- Outputs:
  - processed AnnData objects  
  - graph structures (GraphML)

---

### 2. Multi-omics Integration (scGLUE)  
**`s02_glue.py`**

- Trains **scGLUE model** using RNA + ATAC data  
- Incorporates prior regulatory graph  
- Generates:
  - integrated cell embeddings  
  - feature (gene/peak) embeddings  
- Performs:
  - pretraining + fine-tuning  
  - integration consistency evaluation  

---

### 3. Gene–TF Network Inference  
**`s03_infer_gene_tf.py`**

- Infers regulatory relationships using multiple signals:

#### Gene–Peak Links
- Based on:
  - genomic distance  
  - eQTL evidence  
  - GLUE feature similarity  

#### TF Binding
- Uses ChIP-seq data to map TF binding:
  - promoter flanks  
  - peaks  

#### Cis-regulatory Ranking
- Computes TF–gene regulatory strength using:
  - distance-based links  
  - eQTL-supported links  
  - GLUE-inferred links  

#### SCENIC Integration
- Builds co-expression network (GRN)  
- Applies cisTarget pruning  
- Produces final TF–target networks  

Outputs:
- gene–peak connections  
- TF–gene rankings  
- final GRNs (GraphML, CSV, plots)

---

### 4. Comparison and Evaluation  
**`s04_compare.py`**

- Compares GRNs derived from different strategies:
  - distance-based  
  - eQTL-based  
  - GLUE-based  
- Evaluates consistency and differences across methods  

---

## 🔄 Workflow

---

## 📊 Inputs

- scRNA-seq data (`.h5ad`)
- scATAC-seq data (`.h5ad`)
- Gene annotation (GTF)
- eQTL summary data
- TF binding data (e.g., ChIP-seq BED)
- Optional: precomputed priors

---

## 📤 Outputs

- Processed multiome datasets  
- Prior regulatory graphs  
- Integrated embeddings (cells + features)  
- Gene–peak regulatory links  
- TF–gene interaction networks  
- Final GRNs (GraphML / CSV / visualisations)

---

## ⚙️ Dependencies

- Python (3.8+)
- scanpy  
- anndata  
- scglue  
- muon  
- networkx  
- numpy / pandas / scipy  
- seaborn / matplotlib  
- pySCENIC (via Singularity)

---

## 🚀 Usage

Scripts are designed to be run sequentially:

1. `s01_preprocessing.py`  
2. `s02_glue.py`  
3. `s03_infer_gene_tf.py`  
4. `s04_compare.py`  

Each step produces intermediate outputs used by downstream analyses.

---

