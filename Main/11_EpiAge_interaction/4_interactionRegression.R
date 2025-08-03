#!/usr/bin/env Rscript

#' Genotype × Epigenetic Age Interaction Analysis
#' 
#' This script tests for interactions between genetic variants and epigenetic age
#' on chromatin accessibility. It performs linear regression analysis using 
#' previously identified caQTLs and epigenetic age estimates.
#' 
#' Usage: Rscript 4_interactionRegression.R <cell_type> <chromosome>
#' 
#' Author: Peter C Allen

# Load required libraries
suppressPackageStartupMessages({
  library(lme4)
  library(readr)
  library(lmerTest)
  library(dplyr)
  library(data.table)
  library(tibble)
  library(purrr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 4_interactionRegression.R <cell_type> <chromosome>")
}

cell_type <- args[1]
chr <- args[2]

cat("Analyzing", cell_type, "on", chr, "\n")

# Load caQTL pairs for this cell type
cat("Loading caQTL data...\n")
caQTL_pairs <- read.csv("data/TenK10K.caQTL.1Mb.final_summary_qval005_15jun2025.csv") %>%
  filter(cell_types == cell_type)

cat("Found", nrow(caQTL_pairs), "caQTL pairs for", cell_type, "\n")

# Load chromatin accessibility data
cat("Loading accessibility data...\n")
expression_matrix <- fread(paste0("data/expressionBeds/", cell_type, "/ExpressionBeds/", chr, ".bed.gz")) %>% 
  as.data.frame()

# Format accessibility matrix
rownames(expression_matrix) <- expression_matrix$phenotype_id
expression_matrix <- expression_matrix[, -(1:4)]  # Remove coordinate columns
expression_matrix <- expression_matrix[na.omit(match(caQTL_pairs$phenotype_id, rownames(expression_matrix))), ]

cat("Processing", nrow(expression_matrix), "peaks\n")

# Load and format covariates
cat("Loading covariates...\n")
covariates_df <- read_csv(paste0("data/expressionBeds/", cell_type, "/covariates.txt"), trim_ws = TRUE) %>% 
  as.data.frame() %>% t()

colnames(covariates_df) <- covariates_df[1, ]
covariates_df <- covariates_df[-1, ]
rownames(covariates_df) <- gsub("X", "", rownames(covariates_df))

# Process covariates
covariates_df <- as.data.frame(covariates_df) %>%
  mutate(across(-sex, as.numeric)) %>%  
  mutate(sex = as.factor(ifelse(as.numeric(as.character(sex)) == 1, 1, 2)))

# Load epigenetic age standard deviation (used as covariate)
cat("Loading epigenetic age variability data...\n")
epiage_sd <- read.delim(paste0("data/20250630_epiages/", cell_type, "_epiageSD.txt")) %>%
  filter(!is.na(sd_epiage))

# Merge epigenetic age SD with covariates
covariates_df <- covariates_df %>%
  rownames_to_column(var = "ID") %>%
  left_join(epiage_sd, by = c("ID" = "individual")) %>%
  column_to_rownames(var = "ID") %>%
  filter(!is.na(sd_epiage))

# Load median epigenetic age data (the main effect of interest)
cat("Loading median epigenetic age data...\n")
epiage <- read.delim(paste0("data/20250630_epiages/", cell_type, "_epiage.txt")) %>%
  filter(!is.na(median_epiage))

# Merge median epigenetic age with covariates
covariates_df <- covariates_df %>%
  rownames_to_column(var = "ID") %>%
  left_join(epiage, by = c("ID" = "individual")) %>%
  column_to_rownames(var = "ID")

cat("Covariates prepared for", nrow(covariates_df), "individuals\n")

# Import genotypes
# geno_df <- fread(paste0("data/genotype/hg38/genotypes_chr", chr, ".raw"), header = TRUE, showProgress = TRUE, data.table = FALSE)
geno_df <- fread(paste0("data/genotype/wgs/TenK10K_TOB_ATAC_renamed_", chr, "_common_variants.raw"), header = TRUE, showProgress = TRUE, data.table = FALSE)

geno_df$ID <- geno_df$IID

# Ensure unique column names -- duplication error
colnames(geno_df) <- make.unique(colnames(geno_df))

# Select relevant columns
geno_df <- geno_df %>% select(-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE)
colnames(geno_df) <- gsub("_[GCAT]$", "", colnames(geno_df))
geno_df <- geno_df %>% select(ID, na.omit(match(caQTL_pairs$variant_id, colnames(geno_df)))) # 10631/11666 (91% - ASDC Chr1)

# Count duplicate column names in geno_df
duplicate_colnames <- colnames(geno_df)[duplicated(colnames(geno_df))]
cat("Number of duplicate Variants in Genotype DF:", length(duplicate_colnames), "\n")

# Remove duplicate column names from geno_df
geno_df <- geno_df[, !duplicated(colnames(geno_df))]
geno_df <- geno_df %>% column_to_rownames(var = "ID")

# Only keep individuals with data
geno_df <- geno_df[na.omit(match(rownames(covariates_df), rownames(geno_df))), , drop = FALSE]
# rownames(geno_df) <- rownames(covariates_df)

peaks_keep <- caQTL_pairs$phenotype_id[na.omit(match(colnames(geno_df), caQTL_pairs$variant_id))]
variants_keep <- caQTL_pairs$variant_id[na.omit(match(colnames(geno_df), caQTL_pairs$variant_id))]

expression_matrix <- expression_matrix[peaks_keep,rownames(covariates_df)]
geno_df <- geno_df[rownames(covariates_df), variants_keep]

# Transpose expression matrix
expr_matrix_t <- as.data.frame(t(expression_matrix))
expr_matrix_t$ID <- rownames(expr_matrix_t)

# Combine expression and genotype matrix into one data frame
geno_df$ID <- rownames(geno_df)
covariates_df$ID <- rownames(covariates_df)

merged_df <- reduce(
  list(covariates_df, geno_df, expr_matrix_t),
  function(x, y) merge(x, y, by = "ID")
)

# ---- Peak to Variant Mapping ----
peak_names <- rownames(expression_matrix)  # from original expr_matrix
snp_names <- colnames(geno_df)[-which(colnames(geno_df) == "ID")]
peak_snp_map <- data.frame(peak = peak_names, snp = snp_names)

# ---- Loop through peaks and run regression ----
results_list <- list()

for (i in seq_len(nrow(peak_snp_map))) {
  peak <- peak_snp_map$peak[i]
  snp <- peak_snp_map$snp[i]

  merged_df$expression <- merged_df[[peak]]
  merged_df$genotype_snp <- merged_df[[snp]]

  nGenotypes <- length(table(merged_df$genotype_snp))
  if (nGenotypes == 3) {
    T<- table(merged_df$genotype_snp)
    MAF <- (sum(T[2])+sum(T[3])+sum(T[3]))/(2*sum(T))
  }

  model <- tryCatch({
    lm(expression ~ sex + age + sd_epiage + genotype_pc1 + genotype_pc2 + genotype_pc3 +
         genotype_pc4 + genotype_pc5 + genotype_pc6 +
         genotype_snp*median_epiage, data = merged_df)
    # lm(expression ~ sex + age + genotype_pc1 + genotype_pc2 + genotype_pc3 +
    #      genotype_pc4 + genotype_pc5 + genotype_pc6 +
    #      genotype_snp, data = merged_df)
  }, error = function(e) NULL)

  if (!is.null(model)) {
    coef_df <- as.data.frame(coef(summary(model)))
    coef_df$peak <- peak
    coef_df$snp <- snp
    coef_df$term <- rownames(coef_df)
    coef_df$nGenotypes <- nGenotypes
    coef_df$MAF <- MAF
    results_list[[peak]] <- coef_df
  }
}

# Combine
final_results <- bind_rows(results_list)

# Filter for the interaction term
interaction_results <- final_results %>% filter(term == "genotype_snp:median_epiage")
interaction_results <- interaction_results[order(interaction_results$`Pr(>|t|)`), ]
interaction_results$adj_pval <- p.adjust(interaction_results$`Pr(>|t|)`, method = "BH")
# interaction_results <- final_results %>% filter(term == "genotype_snp")
# interaction_results <- interaction_results[order(interaction_results$`Pr(>|t|)`), ]
# interaction_results$adj_pval <- p.adjust(interaction_results$`Pr(>|t|)`, method = "fdr")

# Print number of significant results
cat("\n\n")
cat("Number of significant results (p.adj < 0.05):", sum(interaction_results$adj_pval < 0.05), "\n")

# Save results
write.csv(final_results, paste0("output/4-regression-results/", "regression_results_", cell_type, "_", chr, ".csv"), row.names = FALSE)
write.csv(interaction_results, paste0("output/4-regression-results/", "interaction_results_", cell_type, "_", chr, ".csv"), row.names = FALSE)

