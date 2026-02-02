# Cell type annotation for scATAC-seq dataset

This script processes single-cell ATAC-seq data by loading pre-processed scATAC and multiome datasets, quantifying multiome peaks in the scATAC cells, performing LSI and UMAP dimensionality reduction, and integrating the datasets to generate co-embeddings. It then maps the scATAC query onto the multiome reference for cell type annotation, visualizes the results, and saves the annotated dataset and metadata for downstream analysis.
