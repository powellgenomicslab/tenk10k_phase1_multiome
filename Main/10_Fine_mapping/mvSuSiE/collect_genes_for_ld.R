#!/usr/bin/env Rscript
library(data.table)

# Load gene locations and eGene info
gene_loc <- fread("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/coloc_1mb/gene_pos.csv")
eGene_info <- fread("/directflow/SCCGGroupShare/projects/angxue/proj/tenk10k/eQTL/coloc_gwas_threshold/1Mb/NewDiseases/ibd.csv")

anchor_trait <- "ibd_liu2023"
celltypes <- list.dirs("coloc_1mb/eQTL1MbRuns/caQTL2eQTL/Coloc_Results", full.names = FALSE, recursive = FALSE)

# Process each chromosome
for (chr_num in 1:22) {
    all_genes <- c()
    
    # Collect genes from all cell types
    for (anchor_celltype in celltypes) {
        files <- c(
            paste0("coloc_1mb/eQTL1MbRuns/caQTL2eQTL/Coloc_Results/", anchor_celltype, "/chr", chr_num, ".csv"),
            paste0("coloc_1mb/caQTL2GWAS/DiseaseTraits/Coloc_Results/", anchor_trait, "/", anchor_celltype, "/chr", chr_num, ".csv"),
            paste0("coloc_1mb/eQTL1MbRuns/DiseaseTraits/Coloc_Results/", anchor_trait, "/", anchor_celltype, "/chr", chr_num, ".csv")
        )
        
        if (all(file.exists(files))) {
            # Load coloc results
            coloc_eQTL_caQTL <- fread(files[1])
            coloc_caQTL_GWAS <- fread(files[2])
            coloc_eQTL_GWAS <- fread(files[3])
            
            # Filter eGene info for this chromosome and cell type
            eGene_info_sub <- eGene_info[chrom == paste0("chr", chr_num) & celltype == anchor_celltype]
            # Add GWAS p-values to eQTL-GWAS coloc results
            coloc_eQTL_GWAS[, min_p_gwas := eGene_info_sub[match(coloc_eQTL_GWAS$gene, eGene_info_sub$gene), min_p_gwas]]
            # Filter for genome-wide significant GWAS hits
            coloc_eQTL_GWAS <- coloc_eQTL_GWAS[min_p_gwas < 1e-5]
            # Filter for H3+H4 > 0.8
            coloc_eQTL_GWAS_filter <- coloc_eQTL_GWAS[(PP.H3.abf + PP.H4.abf) > 0.8]
            
            # Collect genes
            all_genes <- c(all_genes, coloc_eQTL_GWAS_filter$gene)
        }
    }
    
    # Get unique genes
    unique_genes <- unique(all_genes)
    if (length(unique_genes) == 0) next
    
    # Get gene centers from gene_loc
    gene_loc_chr <- gene_loc[gene_id %in% unique_genes]
    
    # Calculate windows (1MB on each side of gene center)
    gene_df <- data.frame(
        gene = gene_loc_chr$gene_id,
        chr = paste0("chr", chr_num),
        centre = gene_loc_chr$centre,
        window_start = pmax(gene_loc_chr$centre - 1000000, 1),
        window_end = gene_loc_chr$centre + 1000000
    )
    
    # Save gene coordinates
    dir.create(paste0("ld_matrix/LD_chr", chr_num, "/genes"), recursive = TRUE, showWarnings = FALSE)
    write.table(gene_df, paste0("ld_matrix/LD_chr", chr_num, "/genes/chr", chr_num, "_gene_coords.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    print(paste("Chr", chr_num, ":", nrow(gene_df), "genes with genome-wide significant GWAS hits and H3+H4 > 0.8"))
}
