#!/usr/bin/env Rscript
# Simplified mvSuSiE analysis pipeline
library(mvsusieR)
library(data.table)
library(dplyr)
library(qvalue)
library(mashr)
library(ggplot2)


# Set working directory
setwd("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie")

# ========== MAIN ANALYSIS FUNCTION ==========
doFinemapping <- function(anchor_gene) {
    print(paste("Processing gene:", anchor_gene))

    # Load SNP/variant info for anchor_gene (eGene)
    ld_base <- "/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/ld_matrix"
    variants_file <- paste0(ld_base, "/LD_chr", chr_num, "/ld_matrices_1000000bp/", anchor_gene, ".variants")
    ld_file <- paste0(ld_base, "/LD_chr", chr_num, "/ld_matrices_1000000bp/", anchor_gene, ".ld")
    if (!file.exists(variants_file)) {
        stop(data.frame(gene = anchor_gene, status = "Variants file not found"))
    }
    if (!file.exists(ld_file)) {
        stop(data.frame(gene = anchor_gene, status = "LD file not found"))
    }
    variant_ids <- fread(variants_file, header=FALSE)[[1]]
    ld_matrix <- as.matrix(fread(ld_file, header=FALSE))
    print(paste("Found", length(variant_ids), "variants in LD matrix"))
    
    # get overlapping snps across 3 layers (gene, caPeaks, gwas) - to do: refactor
    eqtl_region <- eqtl_data[MarkerID %in% variant_ids]
    caqtl_region <- caqtl_data[variant_id %in% variant_ids]
    gwas_region <- gwas_data[SNP %in% variant_ids]
    caqtl_region_sig <- caqtl_region[fdr < 0.05] # only caPeaks
    caqtl_region_sig_snps <- unique(caqtl_region_sig$variant_id)
    gwas_region_sig <- gwas_region[P < 1e-5]
    overlapping_snps <- Reduce(intersect, list(unique(gwas_region_sig$SNP), unique(eqtl_region$MarkerID), caqtl_region_sig_snps))
    if (nrow(gwas_region_sig) == 0) {
        print("No significant GWAS SNPs found in this region")
        stop(data.frame(gene = anchor_gene, status = "No significant GWAS SNPs"))
    } else {
        print(paste("Found", nrow(gwas_region_sig), "significant GWAS SNPs in this region - proceeding to fine-mapping"))
    }
    if (length(overlapping_snps) < 5) {
        print("No enough SNPs found in the region")
        stop(data.frame(gene = anchor_gene, status = "No overlapping SNPs"))
    } else {
        print(paste("Found", length(overlapping_snps), "overlapping SNPs in the region"))
    }
    
    eqtl_region_sub <- eqtl_region[(MarkerID %in% overlapping_snps) & (gene == anchor_gene)] # all snp gene expression pairs
    caqtl_region_sub <- caqtl_region[(variant_id %in% overlapping_snps) & (phenotype_id %in% caqtl_region_sig$phenotype_id) & (phenotype_id %in% coloc_eQTL_caQTL[gene == anchor_gene]$peak) ] # only caQTL sig pairs
    gwas_region_sub <- gwas_region[gwas_region$SNP %in% overlapping_snps] # only gwas sig pairs
    # Keep only peaks with complete SNP data
    n_required <- length(overlapping_snps)
    peaks_complete <- caqtl_region_sub[variant_id %in% overlapping_snps,
                                     .(n_snps = uniqueN(variant_id)),
                                     by = phenotype_id][n_snps == n_required]$phenotype_id
    caqtl_region_sub <- caqtl_region_sub[phenotype_id %in% peaks_complete]
    if (nrow(eqtl_region_sub) == 0 || nrow(caqtl_region_sub) == 0) {
        print("No eQTL or caQTL data available for this peak after selecting top traits.")
        stop(data.frame(gene = anchor_gene, status = "No eQTL or caQTL data"))
    }



    # filter and match LD matrix to the overlapping SNPs
    snp_indices <- match(overlapping_snps, variant_ids)
    ld_matrix_filtered <- ld_matrix[snp_indices, snp_indices]
    print(paste("Match LD matrix to", nrow(ld_matrix_filtered), "SNPs"))
    # Get traits and check completeness
    eqtl_traits <- unique(eqtl_region_sub$gene)
    caqtl_traits <- unique(caqtl_region_sub$phenotype_id)
    # Build trait list
    all_traits <- c(paste0("eQTL_", eqtl_traits), 
                    paste0("caQTL_", caqtl_traits),
                    "GWAS")


    #Z matrix
    Z <- matrix(0, nrow = length(overlapping_snps), ncol = length(all_traits))
    rownames(Z) <- overlapping_snps
    colnames(Z) <- all_traits
    col_idx <- 1
    # eQTL z-scores
    for (gene_name in eqtl_traits) {
        for (i in 1:length(overlapping_snps)) {
            snp <- overlapping_snps[i]
            eqtl_row <- eqtl_region_sub[MarkerID == snp & gene == gene_name]#####
            if (nrow(eqtl_row) == 1) {
                Z[i, col_idx] <- eqtl_row$z[1]
            } else {
                stop(data.frame(gene = anchor_gene, status = "Z matrix error check"))
            }
        }
        col_idx <- col_idx + 1
    }
    # caQTL z-scores
    for (peak_id in caqtl_traits) {
        for (i in 1:length(overlapping_snps)) {
            snp <- overlapping_snps[i]
            caqtl_row <- caqtl_region_sub[variant_id == snp & phenotype_id == peak_id]
            if (nrow(caqtl_row) == 1) {
                Z[i, col_idx] <- caqtl_row$z[1]
            } else {
                stop(data.frame(gene = anchor_gene, status = "Z matrix error check"))
            }
        }
        col_idx <- col_idx + 1
    }
    # GWAS z-scores
    for (i in 1:length(overlapping_snps)) {
        snp <- overlapping_snps[i]
        gwas_row <- gwas_region_sub[SNP == snp]
        if (nrow(gwas_row) == 1) {
            Z[i, col_idx] <- gwas_row$z[1]
        } else {
            stop(data.frame(gene = anchor_gene, status = "Z matrix error check"))
        }
    }
    
    # Run mvSuSiE
    print("Running mvSuSiE...")
    prior <- create_mixture_prior(R=ncol(Z), null_weight=0.0)
    n_samples <- 922  # Conservative: use caQTL sample size
    colnames(ld_matrix_filtered) <- rownames(ld_matrix_filtered) <- overlapping_snps
    set.seed(2333)
    fit <- mvsusie_rss(
        Z = Z,
        R = ld_matrix_filtered,
        N = n_samples,
        prior_variance = prior,
        L = 10,
        tol = 0.01,
        max_iter = 100
    )

    # Collect results
    if (length(fit$sets$cs) > 0) {
        # PIP plot
        out <- mvsusie_plot(fit, chr=chr_num)
        pip_plot <- out$pip_plot
        dir_path <- paste0("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/results_1mb/plots/", anchor_trait, "/", anchor_celltype, "/chr", chr_num)
        if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
        ggsave(paste0("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/results_1mb/plots/", anchor_trait, "/", anchor_celltype, "/chr", chr_num, "/", anchor_gene, "_pip_plot.pdf"), pip_plot, width = 6, height = 3)
        print("PIP plot saved")

        # Effect plot
        out <- mvsusie_plot(fit, chr=chr_num, add_cs=TRUE)
        effect_plot <- out$effect_plot
        dir_path <- paste0("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/results_1mb/plots/", anchor_trait, "/", anchor_celltype, "/chr", chr_num)
        if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
        ggsave(paste0("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/results_1mb/plots/", anchor_trait, "/", anchor_celltype, "/chr", chr_num, "/", anchor_gene, "_effect_plot.pdf"), effect_plot, width = 5, height = 6)
        print("effect plot saved")

        results_list <- list()
        for (i in 1:length(fit$sets$cs)) {
            cs_name <- names(fit$sets$cs)[i]
            cs_indices <- fit$sets$cs[[i]]
            
            snp_pips <- fit$pip[cs_indices]
            names(snp_pips) <- overlapping_snps[cs_indices]
            lfsr <- fit$single_effect_lfsr[as.numeric(substring(cs_name,2)), ]
            sig_traits <- colnames(Z)[lfsr < 0.01]
            results_list[[cs_name]] <- list(
                snps = snp_pips,  # Named vector with PIPs
                top_snp = names(snp_pips)[which.max(snp_pips)],
                coverage = round(fit$sets$coverage[i], 3),
                purity = round(fit$sets$purity[i, "min.abs.corr"], 3),
                lfsr = lfsr,
                traits_affected = sig_traits
            )
        }
        return(results_list)
    } else {
        stop(data.frame(gene = anchor_gene, status = "No credible sets found"))
    }
}

# ========== ARGUMENTS ==========
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
    chr_num <- as.numeric(args[1])
    anchor_celltype <- args[2]
    anchor_trait <- args[3]
} else {
    chr_num <- 11
    anchor_celltype <- "CD14_Mono"
    anchor_trait <- "ibd_liu2023"
}

# ========== LOAD DATA ==========
print("Loading data...")
# Load eQTL data and match chromosome
eqtl_file <- paste0("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/gwas/saige-1mb-summ-stats/",anchor_celltype,"_common_all_cis_raw_pvalues_1000000bp.tsv")
eqtl_data <- fread(eqtl_file)
eqtl_data <- eqtl_data[eqtl_data$CHR == chr_num] 
eqtl_data$fdr <- qvalue(eqtl_data$'p.value')$qvalues
eqtl_data$z <- eqtl_data$BETA / eqtl_data$SE 

# Load caQTL data for the chromosome
caqtl_file <- paste0("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/Final_caQTL/NewGenotypes/Runs/", 
                     anchor_celltype, "/Results/TenK10K.cis_qtl_pairs.chr", chr_num, ".csv")
print("Loading caQTL data...")
caqtl_data <- fread(caqtl_file)
caqtl_data$fdr <- qvalue(caqtl_data$'pval_nominal')$qvalues
caqtl_data$z <- caqtl_data$slope / caqtl_data$slope_se

# Load GWAS data 
print("Loading GWAS data...")
gwas_file <- paste0("/directflow/SCCGGroupShare/projects/angxue/data/GWAS_summary/IBD_Liu_et_al_2023/updated_02AUG2025/IBD_EAS_EUR_meta_Liu_et_al_QCed.txt.gz")
gwas_data <- fread(gwas_file)
gwas_data[, SNP := sub("^chr", "", SNP)]
gwas_data <- gwas_data[gwas_data$CHR == chr_num]
gwas_data$z <- gwas_data$b / gwas_data$se
print("Data loaded successfully")

# ========== READ PEAK ==========
coloc_eQTL_GWAS <- fread(paste0("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/coloc_1mb/eQTL1MbRuns/DiseaseTraits/Coloc_Results/", anchor_trait,
"/",anchor_celltype, "/chr", chr_num, ".csv"))
eGene_info <- fread("/directflow/SCCGGroupShare/projects/angxue/proj/tenk10k/eQTL/coloc_gwas_threshold/1Mb/NewDiseases/ibd.csv")
eGene_info <- eGene_info[chrom==paste0("chr",{chr_num}) & celltype==anchor_celltype] 
coloc_eQTL_GWAS[, min_p_gwas := eGene_info[match(coloc_eQTL_GWAS$gene, eGene_info$gene), min_p_gwas]]
coloc_eQTL_GWAS <- coloc_eQTL_GWAS[min_p_gwas < 1e-5]
coloc_eQTL_GWAS_filter <- coloc_eQTL_GWAS[PP.H4.abf > 0.8]
gene_loc <- fread("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/coloc_1mb/gene_pos.csv")
gene_loc <- gene_loc[gene_id %in% coloc_eQTL_GWAS_filter$gene]
gene_loc[, window_start := pmax(centre - 1000000, 1)]
gene_loc[, window_end := centre + 1000000]

coloc_eQTL_caQTL <- fread(paste0("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/coloc_1mb/eQTL1MbRuns/caQTL2eQTL/Coloc_Results/",
anchor_celltype,"/chr",chr_num,".csv"))
coloc_eQTL_caQTL <- coloc_eQTL_caQTL[PP.H4.abf > 0.8]
coloc_caQTL_GWAS <- fread(paste0("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/coloc_1mb/caQTL2GWAS/DiseaseTraits/Coloc_Results/", anchor_trait,
"/", anchor_celltype, "/chr", chr_num, ".csv"))

print(paste0("Performing mvsusie finemapping for: ", nrow(gene_loc), " genes..."))
genes <- unique(gene_loc$gene_id)


all_results <- list()
summary_df <- list()

for (i in genes){
    result <- try(doFinemapping(i))
    if (!inherits(result, "try-error")) {
        all_results[[i]] <- result  # Comprehensive results by gene
        
        # Build concise summary
        if (length(result) > 0) {
            for (cs_name in names(result)) {
                cs <- result[[cs_name]]
                summary_df[[length(summary_df) + 1]] <- data.frame(
                    gene = i,
                    cs = cs_name,
                    top_snp = cs$top_snp,
                    pip = round(max(cs$snps), 3),
                    coverage = cs$coverage,
                    purity = cs$purity,
                    n_snps = length(cs$snps),
                    affects_all = any(grepl("^eQTL_", cs$traits_affected)) &
                    any(grepl("^caQTL_", cs$traits_affected)) &
                    ("GWAS" %in% cs$traits_affected)
                )
            }
        }
    }
}

# Save outputs
dir_path <- paste0("/directflow/SCCGGroupShare/projects/lawhua/projects/proj_ATAC/mvsusie/results_1mb/summary/", anchor_trait, "/", anchor_celltype)
if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)

# Save comprehensive RDS
saveRDS(all_results, paste0(dir_path, "/chr", chr_num, "_comprehensive.rds"))
print(paste("Saved comprehensive RDS for", length(all_results), "genes"))

# Save concise CSV
if (length(summary_df) > 0) {
    final_table <- do.call(rbind, summary_df)
    write.csv(final_table, paste0(dir_path, "/chr", chr_num, "_summary.csv"), row.names = FALSE)
    print(paste("Saved concise CSV with", nrow(final_table), "credible sets"))
} else {
    print("No credible sets found")
}
