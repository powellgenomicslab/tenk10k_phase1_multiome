# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(patchwork)
library(ggrepel)

trait_name <- "IBD_EAS_EUR"
cell_type_1 <- "CD4_Naive"
cell_type_2 <- "CD16_Mono"
chrNumber <- 16

# Record SNPs for LD matrix
min_pos <- 50200000
max_pos <- 50800000

# Load 1Mb eQTL data
eqtl_df_1Mb_raw_ct_1 <- fread(glue(
  "/g/data/fy54/results/eqtl_1Mb/{cell_type_1}_common_all_cis_raw_pvalues_1000000bp.tsv"
))
eqtl_df_1Mb_raw_ct_2 <- fread(glue(
  "/g/data/fy54/results/eqtl_1Mb/{cell_type_2}_common_all_cis_raw_pvalues_1000000bp.tsv"
))
eqtl_df_1Mb_1 <- eqtl_df_1Mb_raw_ct_1[gene == "ENSG00000121281"]
eqtl_df_1Mb_1 <- eqtl_df_1Mb_1[nchar(Allele1) == 1 & nchar(Allele2) == 1]
eqtl_df_1Mb_1 <- eqtl_df_1Mb_1[POS >= min_pos & POS <= max_pos]

eqtl_df_1Mb_2 <- eqtl_df_1Mb_raw_ct_1[gene == "ENSG00000121274"]
eqtl_df_1Mb_2 <- eqtl_df_1Mb_2[nchar(Allele1) == 1 & nchar(Allele2) == 1]
eqtl_df_1Mb_2 <- eqtl_df_1Mb_2[POS >= min_pos & POS <= max_pos]

eqtl_df_1Mb_3 <- eqtl_df_1Mb_raw_ct_2[gene == "ENSG00000166164"]
eqtl_df_1Mb_3 <- eqtl_df_1Mb_3[nchar(Allele1) == 1 & nchar(Allele2) == 1]
eqtl_df_1Mb_3 <- eqtl_df_1Mb_3[POS >= min_pos & POS <= max_pos]

eqtl_df_1Mb_4 <- eqtl_df_1Mb_raw_ct_2[gene == "ENSG00000167207"]
eqtl_df_1Mb_4 <- eqtl_df_1Mb_4[nchar(Allele1) == 1 & nchar(Allele2) == 1]
eqtl_df_1Mb_4 <- eqtl_df_1Mb_4[POS >= min_pos & POS <= max_pos]

# Load GWAS data
gwas_df <- fread(glue("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/SMR/GWAS/Additional_Set/ma_no_INDEL_RARE_REPEAT/{trait_name}/chr{chrNumber}.ma"))
gwas_df[, c("chr_num", "position", "ref", "effect") := tstrsplit(SNP, ":", fixed = TRUE)]
gwas_df[, position := as.numeric(position)]
gwas_df <- gwas_df[position >= min_pos & position <= max_pos]


###
track_colours <- c(eQTL = "#1c529e", caQTL = "#c74a57", GWAS = "#6bae9a")

# ADCY7
eqtl_cs_leads <- c(
  CS1 = "16:50349166:C:T",
  CS2 = "16:50357685:C:T",
  CS3 = "16:50401107:T:C",
  CS4 = "16:50409496:C:T"
)

eqtl_coloc_cs <- "CS1"

top_variant <- eqtl_df_1Mb_1[which.min(p.value)]$MarkerID
eqtl_plot <- eqtl_df_1Mb_1 %>%
  mutate(
    point_size = ifelse(MarkerID == top_variant, 2, 1),
    log_p = -log10(p.value)
  )

eqtl_cs_df <- eqtl_plot %>%
  filter(MarkerID %in% eqtl_cs_leads) %>%
  mutate(
    cs_label   = names(eqtl_cs_leads)[match(MarkerID, eqtl_cs_leads)],
    is_coloc   = cs_label == eqtl_coloc_cs,
    cs_colour  = ifelse(is_coloc, "red", "black")   # red for coloc, black otherwise
  )

p1 <- ggplot(eqtl_plot, aes(x = POS/1e6, y = log_p, size = point_size)) +
  geom_point(colour = track_colours["eQTL"]) +
  scale_size_identity() +
  labs(title = "eQTL: ADCY7 - CD4 TCM", x = " ", y = " ") +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  # ylim(0, 4) +
  theme(
    axis.ticks.x = element_blank(), axis.text.x = element_blank(),
    axis.title.x = element_blank(), axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )

p11 <- p1 +
  geom_point(
    data = eqtl_cs_df,
    aes(x = POS/1e6, y = log_p, colour = cs_colour),
    shape = 18, size = 3,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = eqtl_cs_df,
    aes(x = POS/1e6, y = log_p, label = cs_label, colour = cs_colour),
    size = 3.5,
    fontface = "bold",
    nudge_y = 0.4,
    box.padding = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.4,
    arrow = arrow(length = unit(0.015, "npc"), type = "closed"),
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_identity()

# TENT4B
eqtl_cs_leads <- c(
  CS1 = "16:50331225:C:T"
)

eqtl_coloc_cs <- "CS1"

top_variant <- eqtl_df_1Mb_2[which.min(p.value)]$MarkerID
eqtl_plot <- eqtl_df_1Mb_2 %>%
  mutate(
    point_size = ifelse(MarkerID == top_variant, 2, 1),
    log_p = -log10(p.value)
  )

eqtl_cs_df <- eqtl_plot %>%
  filter(MarkerID %in% eqtl_cs_leads) %>%
  mutate(
    cs_label   = names(eqtl_cs_leads)[match(MarkerID, eqtl_cs_leads)],
    is_coloc   = cs_label == eqtl_coloc_cs,
    cs_colour  = ifelse(is_coloc, "purple", "black")   # red for coloc, black otherwise
  )

p2 <- ggplot(eqtl_plot, aes(x = POS/1e6, y = log_p, size = point_size)) +
  geom_point(colour = track_colours["eQTL"]) +
  scale_size_identity() +
  labs(title = "eQTL: TENT4B - CD4 TCM", x = " ", y = " ") +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  # ylim(0, 4) +
  theme(
    axis.ticks.x = element_blank(), axis.text.x = element_blank(),
    axis.title.x = element_blank(), axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )

p22 <- p2 +
  geom_point(
    data = eqtl_cs_df,
    aes(x = POS/1e6, y = log_p, colour = cs_colour),
    shape = 18, size = 3,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = eqtl_cs_df,
    aes(x = POS/1e6, y = log_p, label = cs_label, colour = cs_colour),
    size = 3.5,
    fontface = "bold",
    nudge_y = 0.4,
    box.padding = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.4,
    arrow = arrow(length = unit(0.015, "npc"), type = "closed"),
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_identity()

# BRD7
eqtl_cs_leads <- c(
  CS1 = "16:50349166:C:T",
  CS2 = "16:50357685:C:T",
  CS3 = "16:50423085:G:T"
)

eqtl_coloc_cs <- "CS1"

top_variant <- eqtl_df_1Mb_3[which.min(p.value)]$MarkerID
eqtl_plot <- eqtl_df_1Mb_3 %>%
  mutate(
    point_size = ifelse(MarkerID == top_variant, 2, 1),
    log_p = -log10(p.value)
  )

eqtl_cs_df <- eqtl_plot %>%
  filter(MarkerID %in% eqtl_cs_leads) %>%
  mutate(
    cs_label   = names(eqtl_cs_leads)[match(MarkerID, eqtl_cs_leads)],
    is_coloc   = cs_label == eqtl_coloc_cs,
    cs_colour  = ifelse(is_coloc, "orange", "black")   # red for coloc, black otherwise
  )

p3 <- ggplot(eqtl_plot, aes(x = POS/1e6, y = log_p, size = point_size)) +
  geom_point(colour = track_colours["eQTL"]) +
  scale_size_identity() +
  labs(title = "eQTL: BRD7 - CD16 Mono", x = " ", y = " ") +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  # ylim(0, 4) +
  theme(
    axis.ticks.x = element_blank(), axis.text.x = element_blank(),
    axis.title.x = element_blank(), axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )

p33 <- p3 +
  geom_point(
    data = eqtl_cs_df,
    aes(x = POS/1e6, y = log_p, colour = cs_colour),
    shape = 18, size = 3,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = eqtl_cs_df,
    aes(x = POS/1e6, y = log_p, label = cs_label, colour = cs_colour),
    size = 3.5,
    fontface = "bold",
    nudge_y = 0.4,
    box.padding = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.4,
    arrow = arrow(length = unit(0.015, "npc"), type = "closed"),
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_identity()

# NOD2
eqtl_cs_leads <- c(
  CS1 = "16:50685832:G:A",
  CS2 = "16:50712015:C:T",
  CS3 = "16:50722863:C:T",
  CS4 = "16:50833006:T:C"
)

eqtl_coloc_cs <- "CS1"

top_variant <- eqtl_df_1Mb_4[which.min(p.value)]$MarkerID
eqtl_plot <- eqtl_df_1Mb_4 %>%
  mutate(
    point_size = ifelse(MarkerID == top_variant, 2, 1),
    log_p = -log10(p.value)
  )

eqtl_cs_df <- eqtl_plot %>%
  filter(MarkerID %in% eqtl_cs_leads) %>%
  mutate(
    cs_label   = names(eqtl_cs_leads)[match(MarkerID, eqtl_cs_leads)],
    is_coloc   = cs_label == eqtl_coloc_cs,
    cs_colour  = ifelse(is_coloc, "pink", "black")   # red for coloc, black otherwise
  )

p4 <- ggplot(eqtl_plot, aes(x = POS/1e6, y = log_p, size = point_size)) +
  geom_point(colour = track_colours["eQTL"]) +
  scale_size_identity() +
  labs(title = "eQTL: NOD2 - CD16 Mono", x = " ", y = " ") +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  ylim(0, 11) +
  theme(
    axis.ticks.x = element_blank(), axis.text.x = element_blank(),
    axis.title.x = element_blank(), axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )

p44 <- p4 +
  geom_point(
    data = eqtl_cs_df,
    aes(x = POS/1e6, y = log_p, colour = cs_colour),
    shape = 18, size = 3,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = eqtl_cs_df,
    aes(x = POS/1e6, y = log_p, label = cs_label, colour = cs_colour),
    size = 3.5,
    fontface = "bold",
    nudge_y = 0.4,
    box.padding = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.4,
    arrow = arrow(length = unit(0.015, "npc"), type = "closed"),
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_identity()

# Generate caQTL_1 locus plot
peak <- "chr2:144659432-144660322"
peak_coords <- gsub("chr2:", "", peak)
peak_start <- as.numeric(strsplit(peak_coords, "-")[[1]][1])
peak_end <- as.numeric(strsplit(peak_coords, "-")[[1]][2])

top_variant <- caqtl_df_1[which.min(pval_nominal)]$variant_id
caqtl_df_1 <- caqtl_df_1 %>%
  mutate(
    point_size = ifelse(variant_id == top_variant, 2, 1),
    log_p = ifelse(pval_nominal == 0, -log10(2.225074e-308), -log10(pval_nominal))
  )

p2 <- ggplot(caqtl_df_1, aes(x = position/1e6, y = log_p)) +
  geom_point(aes(size = point_size), colour = track_colours["caQTL"]) +
  geom_rect(aes(xmin = 144659432/1e6, xmax = 144660322/1e6,
                ymin = -Inf, ymax = Inf),
            fill = "#8B4513", alpha = 0.1) +
  scale_size_identity() +
  labs(title = glue("caQTL: {peak} (H4: 0.942) - CD14 Mono"), x = " ", y = " ") +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  # ylim(0, 21) +
  theme(
    axis.ticks.x = element_blank(), axis.text.x = element_blank(),
    axis.title.x = element_blank(), axis.line.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 12)
  )

#####
# Generate GWAS locus plot
gwas_cs_leads <- c(
  CS1 = "16:50311867:C:T",
  CS2 = "16:50325014:G:C",
  CS3 = "16:50335977:C:T",
  CS4 = "16:50367408:G:A",
  CS5 = "16:50685832:G:A",
  CS6 = "16:50718198:C:T"
)

gwas_coloc_cs_1 <- "CS3"
gwas_coloc_cs_2 <- "CS2"
gwas_coloc_cs_3 <- "CS4"
gwas_coloc_cs_4 <- "CS5"

# Map each coloc CS to its highlight color
coloc_palette <- c("red", "purple", "orange", "pink")
names(coloc_palette) <- c(gwas_coloc_cs_1, gwas_coloc_cs_2,
                          gwas_coloc_cs_3, gwas_coloc_cs_4)

top_variant <- gwas_df[which.min(p)]$SNP
gwas_plot <- gwas_df %>%
  mutate(
    point_size = ifelse(SNP == top_variant, 2, 1),
    log_p = -log10(p)
  )

gwas_cs_df <- gwas_plot %>%
  filter(SNP %in% gwas_cs_leads) %>%
  mutate(
    cs_label  = names(gwas_cs_leads)[match(SNP, gwas_cs_leads)],
    is_coloc  = cs_label %in% names(coloc_palette),
    cs_colour = ifelse(is_coloc, unname(coloc_palette[cs_label]), "black")
  )

p5 <- ggplot(gwas_plot, aes(x = position/1e6, y = log_p, size = point_size)) +
  geom_point(colour = track_colours["GWAS"]) +
  scale_size_identity() +
  labs(
    title = "GWAS: IBD",
    x = " ",
    y = " "
  ) +
  theme_classic() +
  xlim(min_pos/1e6, max_pos/1e6) +
  # ylim(0, 11) +
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

p55 <- p5 +
  geom_point(
    data = gwas_cs_df,
    aes(x = position/1e6, y = log_p, colour = cs_colour),
    shape = 18, size = 3,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = gwas_cs_df,
    aes(x = position/1e6, y = log_p, label = cs_label, colour = cs_colour),
    size = 3.5,
    fontface = "bold",
    nudge_y = 0.8,
    box.padding = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.4,
    arrow = arrow(length = unit(0.015, "npc"), type = "closed"),
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_identity()

# Draw genomic tracks
library(GenomeInfoDb)
library(Signac)
library(scales)

# Exclude the following genes from annotation for better visual appearance

annotations <- readRDS("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/SMR/coloc_SMR_compare/analysis/ETS2_plot/annotations.rds")
pbmc <- readRDS("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/SMR/coloc_SMR_compare/analysis/ETS2_plot/tob_atac_S0220_1_annotated.rds")
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]

annotations <- annotations[annotations$gene_name != "SNX20"]
Annotation(pbmc) <- annotations

# take a look at gene_plot$layers[[2/3/4]]$data
# layer2's data shows the overall position of all genes, can change dodge
region_string <- paste0("chr16-", min_pos, "-", max_pos)
gene_plot <- AnnotationPlot(object = pbmc, region = region_string)

genes_to_move <- c("NKD1", "NOD2", "TENT4B")

# # Update layers 2, 3, 4 (gene bodies and arrows) - set dodge to 2
# for (i in 1:4) {
#   gene_plot$layers[[i]]$data$dodge[gene_plot$layers[[i]]$data$gene_name %in% genes_to_move] <- 2
#   gene_plot$layers[[i]]$mapping$y <- y_values <- gene_plot$layers[[i]]$data$dodge
#   gene_plot$layers[[i]]$mapping$yend <- y_values <- gene_plot$layers[[i]]$data$dodge
# }

# # Update layer 5 (text labels) - set dodge to 2.2
# gene_plot$layers[[5]]$data$dodge[gene_plot$layers[[5]]$data$gene_name %in% genes_to_move] <- 2.2

gene_plot$layers[[5]]$aes_params$size <- 4
gene_plot <- gene_plot +
  scale_x_continuous(
    name = "Position on Chromosome 16 (Mbp)",
    breaks = seq(round(min_pos / 1e5) * 1e5, max_pos, by = 100000),
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


# Merge figures and save
p <- p11 / p22 / p33 / p44 / p55 / gene_plot

ggsave(p, filename = "/home/913/jf1058/tmp/save/locus.png", width = 7.5, height = 7)
