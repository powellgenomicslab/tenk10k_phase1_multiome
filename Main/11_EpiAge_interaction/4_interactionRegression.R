# Load packages
library(lme4)
library(readr)
library(lmerTest)  # for p-values
library(dplyr)
library(data.table)
library(tibble)  # for column_to_rownames
library(purrr)


# Define variables
args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[1]
chr <- args[2]

# Import data
# caQTL_pairs <- read.table(paste0('/directflow/SCCGGroupShare/projects/oscdon/TensorQTL/OldPeaks_NewDemultiplexing/Runs/', cell_type, '/Results/TenK10K.sig_cis_qtl_pairs.', chr, '.csv'), header = TRUE)
# caQTL_pairs <- read.table(paste0('data/caQTL_results/', cell_type, '/Results_1000000/TenK10K.sig_cis_qtl_pairs.', chr, '.csv'), header = TRUE)
# caQTL_pairs <- caQTL_pairs %>% filter(qval<0.05)
caQTL_pairs <- read.csv("data/TenK10K.caQTL.1Mb.final_summary_qval005_15jun2025.csv", header = TRUE) %>%
  filter(cell_types == cell_type)

# Load and format expression data
# expression_matrix <- fread(paste0("/directflow/SCCGGroupShare/projects/oscdon/TensorQTL/OldPeaks_NewDemultiplexing/Runs/", cell_type, "/ExpressionBeds/", chr, ".bed.gz")) %>% as.data.frame()
expression_matrix <- fread(paste0("data/expressionBeds/", cell_type, "/ExpressionBeds/", chr, ".bed.gz")) %>% as.data.frame()
rownames(expression_matrix) <- expression_matrix$phenotype_id
expression_matrix <- expression_matrix[,-(1:4)]
expression_matrix <- expression_matrix[na.omit(match(caQTL_pairs$phenotype_id, rownames(expression_matrix))),]

# Load covariate data
# covariates_df <- read.csv(paste0("/directflow/SCCGGroupShare/projects/oscdon/TensorQTL/OldPeaks_NewDemultiplexing/Runs/", cell_type, "/covariates.txt")) %>% t()
covariates_df <- read_csv(paste0("data/expressionBeds/", cell_type, "/covariates.txt"), trim_ws = TRUE) %>% as.data.frame() %>% t()
colnames(covariates_df) <- covariates_df[1,]
covariates_df <- covariates_df[-1,]
rownames(covariates_df) <- gsub("X", "", rownames(covariates_df))

# Convert columns to appropriate types
covariates_df <- as.data.frame(covariates_df) %>%
  mutate(across(-sex, as.numeric)) %>%  # Convert all columns except 'sex' to numeric
  mutate(sex = as.factor(ifelse(as.numeric(as.character(sex)) == 1, 1, 2)))  # Ensure 'sex' is either 1 or 2

# Import additional covariate, epiage_sd
epiage_sd <- read.delim(paste0("data/20250630_epiages/", cell_type, "_epiageSD.txt")) %>%
  filter(!is.na(sd_epiage))

# Merge sd_epiage from epiage_sd with covariates_df by individual
covariates_df <- covariates_df %>%
  rownames_to_column(var = "ID") %>%
  left_join(epiage_sd, by = c("ID" = "individual")) %>%
  column_to_rownames(var = "ID")

covariates_df <- covariates_df %>% filter(!is.na(sd_epiage))

# Import and merge EpiAge data
epiage <- read.delim(paste0("data/20250630_epiages/", cell_type, "_epiage.txt")) %>%
  filter(!is.na(median_epiage))

covariates_df <- covariates_df %>%
  rownames_to_column(var = "ID") %>%
  left_join(epiage, by = c("ID" = "individual")) %>%
  column_to_rownames(var = "ID")

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
# write.csv(final_results, paste0("output/3_regression/", "SNP_regression_results_", cell_type, "_", chr, ".csv"), row.names = FALSE)
# write.csv(interaction_results, paste0("output/3_regression/", "SNP_specific_results_", cell_type, "_", chr, ".csv"), row.names = FALSE)
