# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(patchwork)

peak_list <- c("chr21:39094096-39094839")
# Load caQTL data
caqtl_raw <- fread(glue(
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/Final_caQTL/NewGenotypes/Runs/CD14_Mono/Results/TenK10K.cis_qtl_pairs.chr21.csv"
), select = c("phenotype_id", "variant_id", "af", "slope", "slope_se", "pval_nominal"))
caqtl_raw <- caqtl_raw[phenotype_id %in% peak_list]
caqtl_raw[, c("chr_num", "position", "ref", "effect") := tstrsplit(variant_id, ":", fixed = TRUE)]
caqtl_raw <- caqtl_raw[nchar(ref) == 1 & nchar(effect) == 1]
caqtl_raw[, position := as.numeric(position)]

# Load old eQTL data
old_eqtl_df <- fread(
  "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/december24_freeze/CD14_Mono/CD14_Mono_common_all_cis_raw_pvalues.tsv"
)[gene == "ENSG00000157557"]
old_eqtl_df <- old_eqtl_df[nchar(Allele1) == 1 & nchar(Allele2) == 1]

# Load new eQTL data
new_eqtl_df <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/src_NewGenotypes/ETS2_1Mb_eQTL_summary/saige-qtl_tenk10k-genome-2-3-eur_output_files_241210_CD14_Mono_chr21_CD14_Mono_ENSG00000157557_cis_1000000bp.txt")
new_eqtl_df <- new_eqtl_df[nchar(Allele1) == 1 & nchar(Allele2) == 1]

# Load GWAS data
gwas_df <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/GWAS/Flagship_disease_trait/ma_no_INDEL_RARE_REPEAT/ibd_EAS_EUR_SiKJEF_meta_IBD/chr21.ma")
gwas_df[, c("chr_num", "position", "ref", "effect") := tstrsplit(SNP, ":", fixed = TRUE)]
gwas_df[, position := as.numeric(position)]

# Record SNPs for LD matrix
min_pos <- 38500000
max_pos <- 39500000

caqtl_raw <- caqtl_raw[position >= min_pos & position <= max_pos]
writeLines(
  unique(caqtl_raw$variant_id),
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/caqtl.txt"
)

old_eqtl_df <- old_eqtl_df[POS >= min_pos & POS <= max_pos]
writeLines(
  unique(old_eqtl_df$MarkerID),
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/old_eqtl.txt"
)

new_eqtl_df <- new_eqtl_df[POS >= min_pos & POS <= max_pos]
writeLines(
  unique(new_eqtl_df$MarkerID),
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/eqtl.txt"
)

gwas_df <- gwas_df[position >= min_pos & position <= max_pos]
writeLines(
  unique(gwas_df$SNP),
  "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/gwas.txt"
)

system("/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink \
        --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr21_common_variants_qced \
        --r2 inter-chr \
        --ld-snp-list /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/caqtl.txt \
        --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/caqtl")

system("/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink \
        --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr21_common_variants_qced \
        --r2 inter-chr \
        --ld-snp-list /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/old_eqtl.txt \
        --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/old_eqtl")

system("/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink \
        --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr21_common_variants_qced \
        --r2 inter-chr \
        --ld-snp-list /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/new_eqtl.txt \
        --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/new_eqtl")

system("/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink \
        --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr21_common_variants_qced \
        --r2 inter-chr \
        --ld-snp-list /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/gwas.txt \
        --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/gwas")

# Load LD matrix
caqtl_ld_matrix <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/caqtl.ld")
old_eqtl_ld_matrix <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/old_eqtl.ld")
new_eqtl_ld_matrix <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/new_eqtl.ld")
gwas_ld_matrix <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/LD_calculation/gwas.ld")

ld_colours <- c("1" = "#1c529e", "2" = "#81c6ea", "3" = "#6bae9a", "4" = "#e8a76b", "5" = "#c74a57", "variant" = "black")
gene_start <- 38805183
gene_end <- 38824955
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
    title = "100kb eQTL: ETS2 (H4: 5e-6)",
    x = " ",
    y = "-log10(p-value)"
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 75) +
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

# Generate new_eQTL locus plot
top_variant <- new_eqtl_df[which.min(p.value)]$MarkerID
new_eqtl_ld <- new_eqtl_ld_matrix[SNP_A == top_variant] %>%
  select(SNP_B, R2) %>%
  rename(MarkerID = SNP_B) %>%
  mutate(R2_bin = as.numeric(cut(R2, breaks = c(0,0.2,0.4,0.6,0.8,1))))
new_eqtl_plot <- new_eqtl_df %>%
  left_join(new_eqtl_ld, by = "MarkerID") %>%
  mutate(
    R2_bin = case_when(
      MarkerID == top_variant ~ "variant",
      is.na(R2_bin) ~ "1",
      TRUE ~ as.character(R2_bin)
    ),
    point_size = ifelse(MarkerID == top_variant, 2, 1),
    log_p = -log10(p.value)
  )

p2 <- ggplot(new_eqtl_plot, aes(x = POS/1e6, y = log_p, size = point_size)) +
  geom_point(aes(colour = factor(R2_bin))) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = "1Mb eQTL: ETS2 (H4: 0.9959)",
    x = " ",
    y = "-log10(p-value)"
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 75) +
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
peak <- "chr21:39094096-39094839"
peak_coords <- gsub("chr21:", "", peak)
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

p3 <- ggplot(caqtl_df, aes(x = position/1e6, y = log_p)) +
  geom_point(aes(colour = factor(R2_bin), size = point_size)) +
  geom_rect(aes(xmin = peak_start/1e6, xmax = peak_end/1e6, ymin = -Inf, ymax = Inf),
            fill = "#8B4513", alpha = 0.1) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = glue("caQTL: {peak} (H4: 0.9953)"),
    x = " ",
    y = "-log10(p-value)"
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, max(caqtl_df$log_p, na.rm = TRUE)+10) +
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

# Generate GWAS locus plot
top_variant <- gwas_df[which.min(p)]$SNP
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
    log_p = -log10(p)
  )

p4 <- ggplot(gwas_plot, aes(x = position/1e6, y = log_p, size = point_size)) +
  geom_point(aes(colour = factor(R2_bin))) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = "GWAS: IBD",
    # x = "Position on chr21 (Mbp)",
    x = " ",
    y = "-log10(p-value)"
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 75) +
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

region_string <- paste0("chr21-", min_pos, "-", max_pos)
gene_plot <- AnnotationPlot(object = pbmc, region = region_string)
gene_plot$layers[[5]]$aes_params$size <- 4
gene_plot <- gene_plot +
  scale_x_continuous(
    name = "Position on Chromosome 21 (Mbp)",
    breaks = seq(round(min_pos / 1e5) * 1e5, max_pos, by = 200000),
    labels = number_format(scale = 1e-6, accuracy = 0.01),
    limits = c(min_pos, max_pos)
  ) +
  labs(x = "Position on Chromosome 21 (Mbp)") +
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )


# Merge figures and save
p <- p1 / p2 / p3 / p4 / gene_plot

ggsave(p, filename = "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/ETS2_plot/repo/figure/p.png", width = 7.5, height = 7)
