# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(patchwork)
library(scales)
library(ggplot2)

peak_list <- c("chr10:45592479-45592785")
cell_type_1 <- "NK"
# Load caQTL data
caqtl_raw_1 <- fread(glue(
  "/g/data/ei56/od8037/TenK10K/caQTLNew/Results/{cell_type_1}/TenK10K.cis_qtl_pairs.chr10.csv"
), select = c("phenotype_id", "variant_id", "af", "slope", "slope_se", "pval_nominal"))
caqtl_raw_1 <- caqtl_raw_1[phenotype_id %in% peak_list]
caqtl_raw_1[, c("chr_num", "position", "ref", "effect") := tstrsplit(variant_id, ":", fixed = TRUE)]
caqtl_raw_1 <- caqtl_raw_1[nchar(ref) == 1 & nchar(effect) == 1]
caqtl_raw_1[, position := as.numeric(position)]

# Load old eQTL data
eqtl_df_1 <- fread(glue(
  "/g/data/fy54/results/eqtl_1Mb/{cell_type_1}_common_all_cis_raw_pvalues_1000000bp.tsv"
))[gene == "ENSG00000165406"]
eqtl_df_1 <- eqtl_df_1[nchar(Allele1) == 1 & nchar(Allele2) == 1]

# cell type 2
cell_type_2 <- "cDC2"
# Load caQTL data
caqtl_raw_2 <- fread(glue(
  "/g/data/ei56/od8037/TenK10K/caQTLNew/Results/{cell_type_2}/TenK10K.cis_qtl_pairs.chr10.csv"
), select = c("phenotype_id", "variant_id", "af", "slope", "slope_se", "pval_nominal"))
caqtl_raw_2 <- caqtl_raw_2[phenotype_id %in% peak_list]
caqtl_raw_2[, c("chr_num", "position", "ref", "effect") := tstrsplit(variant_id, ":", fixed = TRUE)]
caqtl_raw_2 <- caqtl_raw_2[nchar(ref) == 1 & nchar(effect) == 1]
caqtl_raw_2[, position := as.numeric(position)]

# Load old eQTL data
eqtl_df_2 <- fread(glue(
  "/g/data/fy54/results/eqtl_1Mb/{cell_type_2}_common_all_cis_raw_pvalues_1000000bp.tsv"
))[gene == "ENSG00000165406"]
eqtl_df_2 <- eqtl_df_2[nchar(Allele1) == 1 & nchar(Allele2) == 1]

# Record SNPs for LD matrix
min_pos <- 44900000
max_pos <- 46300000

caqtl_raw_1_cut <- caqtl_raw_1[position >= min_pos & position <= max_pos]
eqtl_df_1_cut <- eqtl_df_1[POS >= min_pos & POS <= max_pos]
caqtl_raw_2_cut <- caqtl_raw_2[position >= min_pos & position <= max_pos]
eqtl_df_2_cut <- eqtl_df_2[POS >= min_pos & POS <= max_pos]

caqtl_raw_1_cut[, log_p := -log10(pval_nominal)]
caqtl_raw_1_cut[, celltype := cell_type_1]
caqtl_raw_2_cut[, log_p := -log10(pval_nominal)]
caqtl_raw_2_cut[, celltype := cell_type_2]

eqtl_df_1_cut[, log_p := -log10(p.value)]
eqtl_df_1_cut[, position := POS]
eqtl_df_1_cut[, celltype := cell_type_1]
eqtl_df_2_cut[, log_p := -log10(p.value)]
eqtl_df_2_cut[, position := POS]
eqtl_df_2_cut[, celltype := cell_type_2]

plot_df <- rbind(
  caqtl_raw_1_cut[, .(position, log_p, celltype)], 
  caqtl_raw_2_cut[, .(position, log_p, celltype)]
) %>%
  dplyr::mutate(
    outline = case_when(
      celltype == cell_type_1 & position == 45592613 ~ "black",
      celltype == cell_type_2 & position == 45429520 ~ "black",
      celltype == cell_type_1 & position != 45592613 ~ "transparent",
      celltype == cell_type_2 & position != 45429520 ~ "transparent"
    )
  )

plot_df_eqtl <- rbind(
  eqtl_df_1_cut[, .(position, log_p, celltype)], 
  eqtl_df_2_cut[, .(position, log_p, celltype)]
) %>%
  dplyr::mutate(
    outline = case_when(
      celltype == cell_type_1 & position == 45592613 ~ "black",
      celltype == cell_type_2 & position == 45429520 ~ "black",
      celltype == cell_type_1 & position != 45592613 ~ "transparent",
      celltype == cell_type_2 & position != 45429520 ~ "transparent"
    )
  )
start_bp <- min_pos
end_bp <- max_pos

p1 <- plot_df %>%
  ggplot(aes(x = position, y = log_p)) +
  geom_point(
    data = filter(plot_df, outline != "black"),
    aes(fill = celltype),
    shape = 21,
    size = 1,
    stroke = 0,
    colour = "transparent"
  ) +
  geom_point(
    data = filter(plot_df, outline == "black"),
    aes(fill = celltype),
    shape = 21,
    size = 2,
    stroke = 0.5,
    colour = "black"
  ) +
  geom_rect(aes(xmin = 45591000, xmax = 45594000, ymin = -Inf, ymax = Inf),
            fill = "#b8d8f9", alpha = 0.1) +
  scale_fill_manual(values = c("NK" = "#c5c5c5", "cDC2" = "#cb3a27")) +
  labs(
    x = "Position on Chromosome 10",
    y = expression(-log["10"](italic(P)["caQTL"]))
  ) +
  xlim(start_bp, end_bp) +
  scale_x_continuous(
    name = "Position on Chromosome 10 (Mbp)",
    breaks = seq(round(start_bp / 1e5) * 1e5, end_bp, by = 500000),
    labels = number_format(scale = 1e-6, accuracy = 0.01),
    limits = c(start_bp, end_bp)
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank()
  )

p2 <- plot_df_eqtl %>%
  ggplot(aes(x = position, y = log_p)) +
  geom_point(
    data = filter(plot_df_eqtl, outline != "black"),
    aes(fill = celltype),
    shape = 21,
    size = 1,
    stroke = 0,
    colour = "transparent"
  ) +
  geom_point(
    data = filter(plot_df_eqtl, outline == "black"),
    aes(fill = celltype),
    shape = 21,
    size = 2,
    stroke = 0.5,
    colour = "black"
  ) +
  scale_fill_manual(values = c("NK" = "#c5c5c5", "cDC2" = "#cb3a27")) +
  labs(
    x = "Position on Chromosome 10",
    y = expression(-log["10"](italic(P)["eQTL"]))
  ) +
  xlim(start_bp, end_bp) +
  scale_x_continuous(
    name = "Position on Chromosome 10 (Mbp)",
    breaks = seq(round(start_bp / 1e5) * 1e5, end_bp, by = 500000),
    labels = number_format(scale = 1e-6, accuracy = 0.01),
    limits = c(start_bp, end_bp)
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank()
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

# annotations <- annotations[annotations$gene_name != "TCL1B"]
Annotation(pbmc) <- annotations

region_string <- paste0("chr10-", min_pos, "-", max_pos)

gene_plot <- AnnotationPlot(object = pbmc, region = region_string)   # fresh object

track_map <- c(
  TMEM72 = 1, OR13A1 = 1, WASHC2C = 1, TIMM23 = 1, ANTXRL = 1,
  RASSF4 = 2, ALOX5  = 2, AGAP4   = 2, NCOA4  = 2,
  DEPP1  = 3, MARCH8 = 3, MSMB    = 3,
  ZNF22  = 4, ZFAND4 = 4
)
label_offset <- 0.3

fix_dodge <- function(d, is_label) {
  if (!is.data.frame(d) || !all(c("gene_name", "dodge") %in% names(d))) return(d)
  hit <- d$gene_name %in% names(track_map)
  d$dodge[hit] <- unname(track_map[d$gene_name[hit]]) + if (is_label) label_offset else 0
  d
}

# plot-level data (used by layers that inherit it)
gene_plot$data <- fix_dodge(gene_plot$data, FALSE)

for (i in seq_along(gene_plot$layers)) {
  is_label <- inherits(gene_plot$layers[[i]]$geom, "GeomText")
  gene_plot$layers[[i]]$data <- fix_dodge(gene_plot$layers[[i]]$data, is_label)

  # make sure gene bodies / arrows read y from the dodge column
  if (inherits(gene_plot$layers[[i]]$geom, "GeomSegment")) {
    gene_plot$layers[[i]]$mapping$y    <- aes(y = dodge)$y
    gene_plot$layers[[i]]$mapping$yend <- aes(yend = dodge)$yend
  }
}

for (i in seq_along(gene_plot$layers)) {
  L <- gene_plot$layers[[i]]

  # gene labels
  if (inherits(L$geom, "GeomText")) {
    gene_plot$layers[[i]]$aes_params$size <- 3.3          # text size (default ~3)
  }

  # gene bodies, exons, and arrows
  if (inherits(L$geom, "GeomSegment")) {
    # thicker lines (ggplot2 >= 3.4 uses linewidth; older versions use size)
    lw_name <- if ("linewidth" %in% names(L$aes_params)) "linewidth" else "size"
    old_lw  <- L$aes_params[[lw_name]]
    if (is.null(old_lw)) old_lw <- 0.5
    gene_plot$layers[[i]]$aes_params[[lw_name]] <- old_lw * 1.5   # double thickness

    # bigger arrowheads, only on layers that draw arrows
    if (!is.null(L$geom_params$arrow)) {
      gene_plot$layers[[i]]$geom_params$arrow$length <- unit(0.08, "inches")
    }
  }
}

gene_plot <- gene_plot +
  scale_x_continuous(
    name   = "Position on Chromosome 10 (Mbp)",
    breaks = seq(round(min_pos / 1e5) * 1e5, max_pos, by = 500000),
    labels = number_format(scale = 1e-6, accuracy = 0.01),
    limits = c(min_pos, max_pos)
  ) +
  scale_y_continuous(limits = c(0.7, 4.5)) +   # make room for 4 tracks
  theme(
    axis.line.y  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )

p <- p1 / p2 / gene_plot + plot_layout(heights = c(1, 1, 0.75))
ggsave("/g/data/fy54/jf1058/ei56_nomem/tenk10k_phase1/SMR/coloc_SMR_compare/Find_opp_direction/eQTL_plot/Plots/MARCH8_all.svg", plot = p, width = 10, height = 5, dpi = 300)
