
# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(patchwork)
library(rtracklayer)
library(GenomicRanges)
library(scales)

trait_name <- "IBD_EAS_EUR"
cell_type <- "CD4_TCM"

chrNumber <- 16
peak_to_look <- "chr16:50326863-50327752"
gene_to_look <- "ENSG00000121281"

# Record SNPs for LD matrix
min_pos <- 50270000
max_pos <- 50420000

#
caqtl_df_2 <- fread(glue(
  "/g/data/ei56/od8037/NewGenotypes/caQTL/Runs/{cell_type}/Results/TenK10K.cis_qtl_pairs.chr{chrNumber}.csv"
), select = c("phenotype_id", "variant_id", "af", "slope", "slope_se", "pval_nominal"))
caqtl_df_2 <- caqtl_df_2[phenotype_id == peak_to_look]
caqtl_df_2[, c("chr_num", "position", "ref", "effect") := tstrsplit(variant_id, ":", fixed = TRUE)]
caqtl_df_2 <- caqtl_df_2[nchar(ref) == 1 & nchar(effect) == 1]
caqtl_df_2[, position := as.numeric(position)]
caqtl_df_2 <- caqtl_df_2[position >= min_pos & position <= max_pos]

#
eqtl_df_1Mb_raw <- fread(glue(
  "/g/data/fy54/results/eqtl_1Mb/{cell_type}_common_all_cis_raw_pvalues_1000000bp.tsv"
))
eqtl_df_1Mb_2 <- eqtl_df_1Mb_raw[gene == "ENSG00000121281"]
eqtl_df_1Mb_2 <- eqtl_df_1Mb_2[nchar(Allele1) == 1 & nchar(Allele2) == 1]
eqtl_df_1Mb_2 <- eqtl_df_1Mb_2[POS >= min_pos & POS <= max_pos]

eQTL_finemap_path <- glue("/g/data/fy54/results/eqtl/susie/susie_outputs/{gene_to_look}_{cell_type}_susie_summary.tsv")
eQTL_finemap_df <- fread(eQTL_finemap_path)
eQTL_finemap_df[, position := as.numeric(tstrsplit(SNP, ":")[[2]])]

caQTL_finemap_path <- glue("/g/data/fy54/od8037/TenK10K/FineMapping/SuSiEResults/{cell_type}/chr{chrNumber}.csv")
caQTL_finemap_df <- fread(caQTL_finemap_path)[peak == peak_to_look]
caQTL_finemap_df[, position := as.numeric(tstrsplit(SNP, ":")[[2]])]

# Load GWAS data
gwas_df <- fread(glue("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/SMR/GWAS/Additional_Set/ma_no_INDEL_RARE_REPEAT/{trait_name}/chr{chrNumber}.ma"))
gwas_df[, c("chr_num", "position", "ref", "effect") := tstrsplit(SNP, ":", fixed = TRUE)]
gwas_df[, position := as.numeric(position)]
gwas_df <- gwas_df[position >= min_pos & position <= max_pos]


caqtl_ld_matrix_2 <- fread("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/locus_plot/LD_calculation/caqtl_2.ld")
eqtl_ld_matrix_2 <- fread("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/locus_plot/LD_calculation/eqtl_2.ld")
gwas_ld_matrix <- fread("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/locus_plot/LD_calculation/gwas.ld")

ld_colours <- c("1" = "#1c529e", "2" = "#81c6ea", "3" = "#6bae9a", "4" = "#e8a76b", "5" = "#c74a57", "variant" = "black")
# Generate caQTL locus plot
peak <- peak_to_look
peak_coords <- gsub("chr16:", "", peak)
peak_start <- as.numeric(strsplit(peak_coords, "-")[[1]][1])
peak_end <- as.numeric(strsplit(peak_coords, "-")[[1]][2])

#
top_variant <- caqtl_df_2[which.min(pval_nominal)]$variant_id
caqtl_ld <- caqtl_ld_matrix_2[SNP_A == top_variant] %>%
  select(SNP_B, R2) %>%
  mutate(R2_bin = as.numeric(cut(R2, breaks = c(0,0.2,0.4,0.6,0.8,1))))
setnames(caqtl_ld, "SNP_B", "variant_id")
caqtl_df_2 <- caqtl_df_2 %>%
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

p2 <- ggplot(caqtl_df_2, aes(x = position/1e6, y = log_p)) +
  geom_point(aes(colour = factor(R2_bin), size = point_size)) +
  geom_rect(aes(xmin = peak_start/1e6, xmax = peak_end/1e6, ymin = -Inf, ymax = Inf),
            fill = "#8B4513", alpha = 0.1) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = glue("caQTL: {peak} (p_SMR: 5.65e-10) - CD4 TCM"),
    x = " ",
    y = " "
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, max(caqtl_df_2$log_p, na.rm = TRUE) + 3) +
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

# Generate finemapping caQTL plot

p20 <- ggplot(caQTL_finemap_df, aes(x = position/1e6, y = PIP, 
                                     shape = factor(Credible_Set),
                                     colour = factor(Credible_Set))) +
  geom_point() +
  scale_shape_manual(values = c("1" = 17, "2" = 15)) +
  scale_colour_manual(values = c("1" = "#9B59B6", "2" = "#2ECC71")) +
  labs(
    title = glue("caQTL: {peak} - CD4 TCM (Fine-mapping)"),
    x = " ",
    y = " "
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 0.4) +
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

# Generate eQTL locus plot
top_variant <- eqtl_df_1Mb_2[which.min(p.value)]$MarkerID
eqtl_ld <- eqtl_ld_matrix_2[SNP_A == top_variant] %>%
  select(SNP_B, R2) %>%
  mutate(R2_bin = as.numeric(cut(R2, breaks = c(0,0.2,0.4,0.6,0.8,1))))
setnames(eqtl_ld, "SNP_B", "MarkerID")
eqtl_plot_2 <- eqtl_df_1Mb_2 %>%
  left_join(eqtl_ld, by = "MarkerID") %>%
  mutate(
    R2_bin = case_when(
      MarkerID == top_variant ~ "variant",
      is.na(R2_bin) ~ "1",
      TRUE ~ as.character(R2_bin)
    ),
    point_size = ifelse(MarkerID == top_variant, 2, 1),
    log_p = -log10(p.value)
  )

p12 <- ggplot(eqtl_plot_2, aes(x = POS/1e6, y = log_p, size = point_size)) +
  geom_point(aes(colour = factor(R2_bin))) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +  
  labs(
    title = "eQTL: ADCY7 (p_SMR: 1.10e-9) - CD4 TCM",
    x = " ",
    y = " "
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, max(eqtl_plot_2$log_p, na.rm = TRUE) + 3) +
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

# Generate finemapping eQTL plot

p120 <- ggplot(eQTL_finemap_df, aes(x = position/1e6, y = PIP, 
                                     shape = factor(Credible_Set),
                                     colour = factor(Credible_Set))) +
  geom_point() +
  scale_shape_manual(values = c("1" = 17, "2" = 15)) +
  scale_colour_manual(values = c("1" = "#9B59B6", "2" = "#2ECC71")) +
  labs(
    title = glue("eQTL: ADCY7 - CD4 TCM (Fine-mapping)"),
    x = " ",
    y = " "
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 1) +
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
  mutate(R2_bin = as.numeric(cut(R2, breaks = c(0,0.2,0.4,0.6,0.8,1))))
setnames(gwas_ld, "SNP_B", "SNP")
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

p9 <- ggplot(gwas_plot, aes(x = position/1e6, y = log_p, size = point_size)) +
  geom_point(aes(colour = factor(R2_bin))) +
  scale_colour_manual(values = ld_colours) +
  scale_size_identity() +
  labs(
    title = "GWAS: IBD",
    x = " ",
    y = " "
  ) +
  theme_classic() +
  ylim(0, max(gwas_plot$log_p, na.rm = TRUE) + 3) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(size = 12)
  )

# Draw genomic tracks
library(GenomeInfoDb)
library(Signac)
library(scales)

# Exclude the following genes from annotation for better visual appearance

annotations <- readRDS("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/SMR/coloc_SMR_compare/analysis/ETS2_plot/annotations.rds")
pbmc <- readRDS("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/SMR/coloc_SMR_compare/analysis/ETS2_plot/tob_atac_S0220_1_annotated.rds")
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]

Annotation(pbmc) <- annotations

region_string <- paste0("chr16-", min_pos, "-", max_pos)
gene_plot <- AnnotationPlot(object = pbmc, region = region_string)
gene_plot$layers[[5]]$aes_params$size <- 4
gene_plot <- gene_plot +
  scale_x_continuous(
    name = "Position on Chromosome 16 (Mbp)",
    breaks = seq(round(min_pos / 1e5) * 1e5, max_pos, by = 40000),
    labels = number_format(scale = 1e-6, accuracy = 0.01),
    limits = c(min_pos, max_pos)
  ) +
  labs(x = "Position on Chromosome 16 (Mbp)") +
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )

p <- p2 / p20 / p12 / p120 / p9 / gene_plot + plot_layout(heights = c(1, 0.8, 1, 0.8, 1, 0.8))
ggsave(p, filename = "/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/locus_plot/plot/locus5_finemap_CD4_TCM.png", width = 7.5, height = 6.5, dpi = 300)
