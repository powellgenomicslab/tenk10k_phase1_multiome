# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(patchwork)

trait_name <- "asthma"
cell_type <- "NK_CD56bright"
chrNumber <- 17

# Record SNPs for LD matrix
min_pos <- 39700000
max_pos <- 40100000

# Load eQTL data
# GSDMB - ENSG00000073605 (chr17:39912224)
old_eqtl_df <- fread(glue(
  "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/december24_freeze/{cell_type}/{cell_type}_common_all_cis_raw_pvalues.tsv"
))[gene == "ENSG00000073605"]
old_eqtl_df <- old_eqtl_df[nchar(Allele1) == 1 & nchar(Allele2) == 1]
old_eqtl_df <- old_eqtl_df[POS >= min_pos & POS <= max_pos]

# Load caQTL data
caqtl_raw <- fread(glue(
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/Final_caQTL/NewGenotypes/Runs/{cell_type}/Results/TenK10K.cis_qtl_pairs.chr{chrNumber}.csv"
), select = c("phenotype_id", "variant_id", "af", "slope", "slope_se", "pval_nominal"))
caqtl_raw <- caqtl_raw[phenotype_id == "chr17:39872745-39873334"]
caqtl_raw[, c("chr_num", "position", "ref", "effect") := tstrsplit(variant_id, ":", fixed = TRUE)]
caqtl_raw <- caqtl_raw[nchar(ref) == 1 & nchar(effect) == 1]
caqtl_raw[, position := as.numeric(position)]
caqtl_raw <- caqtl_raw[position >= min_pos & position <= max_pos]

# Load GWAS data
gwas_df <- fread(glue("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/GWAS/Additional_Set/ma_no_INDEL_RARE_REPEAT/{trait_name}/chr{chrNumber}.ma"))
gwas_df[, c("chr_num", "position", "ref", "effect") := tstrsplit(SNP, ":", fixed = TRUE)]
gwas_df[, position := as.numeric(position)]
gwas_df <- gwas_df[position >= min_pos & position <= max_pos]

# Save SNPs and run PLINK LD calculation (comment out once executed)
writeLines(
  unique(old_eqtl_df$MarkerID),
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/old_eqtl.txt"
)

writeLines(
  unique(caqtl_raw$variant_id),
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/caqtl.txt"
)

writeLines(
  unique(gwas_df$SNP),
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/gwas.txt"
)

system("/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink \
        --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr17_common_variants_qced \
        --r2 inter-chr \
        --ld-snp-list /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/caqtl.txt \
        --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/caqtl")

system("/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink \
        --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr17_common_variants_qced \
        --r2 inter-chr \
        --ld-snp-list /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/old_eqtl.txt \
        --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/old_eqtl")

system("/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink \
        --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr17_common_variants_qced \
        --r2 inter-chr \
        --ld-snp-list /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/gwas.txt \
        --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/gwas")


# Load LD matrix
caqtl_ld_matrix <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/caqtl.ld")
old_eqtl_ld_matrix <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/old_eqtl.ld")
gwas_ld_matrix <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/gwas.ld")

ld_colours <- c("1" = "#1c529e", "2" = "#81c6ea", "3" = "#6bae9a", "4" = "#e8a76b", "5" = "#c74a57", "variant" = "black")
# Generate old_eQTL locus plot
top_variant <- old_eqtl_df[which.min(p.value)]$MarkerID
old_eqtl_ld <- old_eqtl_ld_matrix[SNP_A == top_variant] %>%
  select(SNP_B, R2) %>%
  rename(MarkerID = SNP_B) %>%
  mutate(R2_bin = as.numeric(cut(R2, breaks = c(0,0.2,0.4,0.6,0.8,1))))
old_eqtl_plot <- old_eqtl_df %>%
  left_join(old_eqtl_ld, by = "MarkerID") %>%
  mutate(
    R2_bin = case_when(
      MarkerID == top_variant ~ "variant",
      is.na(R2_bin) ~ "1",
      TRUE ~ as.character(R2_bin)
    ),
    point_size = ifelse(MarkerID == top_variant, 2, 1),
    log_p = -log10(p.value)
  )

p1 <- ggplot(old_eqtl_plot, aes(x = POS/1e6, y = log_p, size = point_size)) +
  geom_point(aes(colour = factor(R2_bin))) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = "eQTL: GSDMB (H4: 0.926) - NK CD56bright",
    x = " ",
    y = " "
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 20) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )

# Generate caQTL locus plot
peak <- "chr17:39872745-39873334"
peak_coords <- gsub("chr17:", "", peak)
peak_start <- as.numeric(strsplit(peak_coords, "-")[[1]][1])
peak_end <- as.numeric(strsplit(peak_coords, "-")[[1]][2])
caqtl_df <- caqtl_raw[phenotype_id == peak]

top_variant <- caqtl_df[which.min(pval_nominal)]$variant_id
caqtl_ld <- caqtl_ld_matrix[SNP_A == top_variant] %>%
  select(SNP_B, R2) %>%
  rename(variant_id = SNP_B) %>%
  mutate(R2_bin = as.numeric(cut(R2, breaks = c(0,0.2,0.4,0.6,0.8,1))))
caqtl_df <- caqtl_df %>%
  left_join(caqtl_ld, by = "variant_id") %>%
  mutate(
    R2_bin = case_when(
      variant_id == top_variant ~ "variant",
      is.na(R2_bin) ~ "1",
      TRUE ~ as.character(R2_bin)
    ),
    point_size = ifelse(variant_id == top_variant, 2, 1),
    log_p = case_when(
      pval_nominal == 0 ~ -log10(2.225074e-308),
      TRUE ~ -log10(pval_nominal)
    )
  )

p2 <- ggplot(caqtl_df, aes(x = position/1e6, y = log_p)) +
  geom_point(aes(colour = factor(R2_bin), size = point_size)) +
  geom_rect(aes(xmin = 39872745/1e6, xmax = 39873334/1e6, ymin = -Inf, ymax = Inf),
            fill = "#8B4513", alpha = 0.1) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = glue("caQTL: {peak} (H4: 5e-5) - NK CD56bright"),
    x = " ",
    y = " "
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 10) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 5, 10), labels = c("0", "5", "10")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )


#####
# Do all the stuff for another cell type

cell_type <- "NK_Proliferating"

# Load eQTL data
# GSDMB - ENSG00000073605 (chr17:39912224), ORMDL3 - ENSG00000172057
old_eqtl_df <- fread(glue(
  "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/december24_freeze/{cell_type}/{cell_type}_common_all_cis_raw_pvalues.tsv"
))[gene == "ENSG00000073605"]
old_eqtl_df <- old_eqtl_df[nchar(Allele1) == 1 & nchar(Allele2) == 1]
old_eqtl_df <- old_eqtl_df[POS >= min_pos & POS <= max_pos]

# Load caQTL data
caqtl_raw <- fread(glue(
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/Final_caQTL/NewGenotypes/Runs/{cell_type}/Results/TenK10K.cis_qtl_pairs.chr{chrNumber}.csv"
), select = c("phenotype_id", "variant_id", "af", "slope", "slope_se", "pval_nominal"))
caqtl_raw <- caqtl_raw[phenotype_id == "chr17:39864796-39866052"]
caqtl_raw[, c("chr_num", "position", "ref", "effect") := tstrsplit(variant_id, ":", fixed = TRUE)]
caqtl_raw <- caqtl_raw[nchar(ref) == 1 & nchar(effect) == 1]
caqtl_raw[, position := as.numeric(position)]
caqtl_raw <- caqtl_raw[position >= min_pos & position <= max_pos]

## Save SNPs and run PLINK LD calculation (comment out once executed)
writeLines(
  unique(old_eqtl_df$MarkerID),
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/old_eqtl_ct2.txt"
)

writeLines(
  unique(caqtl_raw$variant_id),
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/caqtl_ct2.txt"
)

system("/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink \
        --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr17_common_variants_qced \
        --r2 inter-chr \
        --ld-snp-list /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/caqtl_ct2.txt \
        --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/caqtl_ct2")

system("/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink \
        --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr17_common_variants_qced \
        --r2 inter-chr \
        --ld-snp-list /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/old_eqtl_ct2.txt \
        --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/old_eqtl_ct2")

# Load LD matrix
caqtl_ld_matrix <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/caqtl_ct2.ld")
old_eqtl_ld_matrix <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/LD_calculation/old_eqtl_ct2.ld")

ld_colours <- c("1" = "#1c529e", "2" = "#81c6ea", "3" = "#6bae9a", "4" = "#e8a76b", "5" = "#c74a57", "variant" = "black")
# Generate old_eQTL locus plot
top_variant <- old_eqtl_df[which.min(p.value)]$MarkerID
old_eqtl_ld <- old_eqtl_ld_matrix[SNP_A == top_variant] %>%
  select(SNP_B, R2) %>%
  rename(MarkerID = SNP_B) %>%
  mutate(R2_bin = as.numeric(cut(R2, breaks = c(0,0.2,0.4,0.6,0.8,1))))
old_eqtl_plot <- old_eqtl_df %>%
  left_join(old_eqtl_ld, by = "MarkerID") %>%
  mutate(
    R2_bin = case_when(
      MarkerID == top_variant ~ "variant",
      is.na(R2_bin) ~ "1",
      TRUE ~ as.character(R2_bin)
    ),
    point_size = ifelse(MarkerID == top_variant, 2, 1),
    log_p = -log10(p.value)
  )

p3 <- ggplot(old_eqtl_plot, aes(x = POS/1e6, y = log_p, size = point_size)) +
  geom_point(aes(colour = factor(R2_bin))) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = "eQTL: GSDMB (H4: 0.890) - NK Proliferating",
    x = " ",
    y = "-log10(p-value)"
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 15) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )

# Generate caQTL locus plot
peak <- "chr17:39864796-39866052"
peak_coords <- gsub("chr17:", "", peak)
peak_start <- as.numeric(strsplit(peak_coords, "-")[[1]][1])
peak_end <- as.numeric(strsplit(peak_coords, "-")[[1]][2])
caqtl_df <- caqtl_raw[phenotype_id == peak]

top_variant <- caqtl_df[which.min(pval_nominal)]$variant_id
caqtl_ld <- caqtl_ld_matrix[SNP_A == top_variant] %>%
  select(SNP_B, R2) %>%
  rename(variant_id = SNP_B) %>%
  mutate(R2_bin = as.numeric(cut(R2, breaks = c(0,0.2,0.4,0.6,0.8,1))))
caqtl_df <- caqtl_df %>%
  left_join(caqtl_ld, by = "variant_id") %>%
  mutate(
    R2_bin = case_when(
      variant_id == top_variant ~ "variant",
      is.na(R2_bin) ~ "1",
      TRUE ~ as.character(R2_bin)
    ),
    point_size = ifelse(variant_id == top_variant, 2, 1),
    log_p = case_when(
      pval_nominal == 0 ~ -log10(2.225074e-308),
      TRUE ~ -log10(pval_nominal)
    )
  )

p4 <- ggplot(caqtl_df, aes(x = position/1e6, y = log_p)) +
  geom_point(aes(colour = factor(R2_bin), size = point_size)) +
  geom_rect(aes(xmin = peak_start/1e6, xmax = peak_end/1e6, ymin = -Inf, ymax = Inf),
            fill = "#CD853F", alpha = 0.1) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = glue("caQTL: {peak} (H4: Not tested) - NK Proliferating"),
    x = " ",
    y = " "
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 10) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 5, 10), labels = c("0", "5", "10")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )


#####
# Generate GWAS locus plot
top_variant <- gwas_df[which.min(P)]$SNP
gwas_ld <- gwas_ld_matrix[SNP_A == top_variant] %>%
  select(SNP_B, R2) %>%
  rename(SNP = SNP_B) %>%
  mutate(R2_bin = as.numeric(cut(R2, breaks = c(0,0.2,0.4,0.6,0.8,1))))
gwas_plot <- gwas_df %>%
  left_join(gwas_ld, by = "SNP") %>%
  mutate(
    R2_bin = case_when(
      SNP == top_variant ~ "variant",
      is.na(R2_bin) ~ "1",
      TRUE ~ as.character(R2_bin)
    ),
    point_size = ifelse(SNP == top_variant, 2, 1),
    log_p = -log10(P)
  )

p5 <- ggplot(gwas_plot, aes(x = position/1e6, y = log_p, size = point_size)) +
  geom_point(aes(colour = factor(R2_bin))) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = "GWAS: asthma",
    # x = "Position on chr21 (Mbp)",
    x = " ",
    y = " "
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 80) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )

# Draw genomic tracks
library(GenomeInfoDb)
library(Signac)
library(scales)

# Exclude the following genes from annotation for better visual appearance

annotations <- readRDS("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/annotations.rds")
pbmc <- readRDS("/directflow/SCCGGroupShare/projects/angxue/proj/multiome/TOB_ATAC/data/QCed/first_56_libraries/tob_atac_S0220_1_annotated.rds")
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]

Annotation(pbmc) <- annotations

region_string <- paste0("chr17-", min_pos, "-", max_pos)
gene_plot <- AnnotationPlot(object = pbmc, region = region_string)
gene_plot$layers[[5]]$aes_params$size <- 4
gene_plot <- gene_plot +
  scale_x_continuous(
    name = "Position on Chromosome 17 (Mbp)",
    breaks = seq(round(min_pos / 1e5) * 1e5, max_pos, by = 200000),
    labels = number_format(scale = 1e-6, accuracy = 0.01),
    limits = c(min_pos, max_pos)
  ) +
  labs(x = "Position on Chromosome 17 (Mbp)") +
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )


# Merge figures and save
p <- p1 / p2 / p3 / p4 / p5 / gene_plot

ggsave(p, filename = "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Find_only1_QTL/plot/locus1.png", width = 7.5, height = 7)
