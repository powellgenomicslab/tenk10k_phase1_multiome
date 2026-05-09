# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(patchwork)
library(scales)

# Set parameters
peak_name <- "chr13:24670806-24672096"
img_name <- gsub("-", "_", gsub(":", "_", peak_name))
chr_num <- gsub("chr", "", strsplit(peak_name, ":")[[1]][1])

snp_cd14_mono <- "13:24570579:C:A"
snp_cd4_tcm <- "13:24671328:T:C"
ld_value <- 4.51791e-07

# Load files
cd4_tcm_df <- fread(glue(
    "/g/data/ei56/od8037/Plotting/QTL_Boxplots/Example_2/CD4_TCM_{img_name}.csv"
)) %>%
    column_to_rownames("phenotype_id") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("donor_id")

cd14_mono_df <- fread(glue(
    "/g/data/ei56/od8037/Plotting/QTL_Boxplots/Example_2/CD14_Mono_{img_name}.csv"
)) %>%
    column_to_rownames("phenotype_id") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("donor_id")

geno_df <- fread(glue(
    "/g/data/ei56/od8037/Plotting/QTL_Boxplots/Example_2/genotype_{img_name}.csv"
)) %>%
    column_to_rownames("snp") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("donor_id")

# Generate dataframes for plotting
plot_df_cd4_tcm <- cd4_tcm_df %>%
    left_join(geno_df, by = "donor_id") %>%
    dplyr::select(all_of(c(peak_name, snp_cd4_tcm))) %>%
    drop_na()

plot_df_cd14_mono <- cd14_mono_df %>%
    left_join(geno_df, by = "donor_id") %>%
    dplyr::select(all_of(c(peak_name, snp_cd14_mono))) %>%
    drop_na()

# Load raw caQTL results
cd4_tcm_caqtl_raw <- fread(
    cmd = glue(
        "head -n 1 /g/data/ei56/od8037/NewGenotypes/caQTL/Runs/CD4_TCM/Results/TenK10K.cis_qtl_pairs.chr{chr_num}.csv && grep '{peak_name}' /g/data/ei56/od8037/NewGenotypes/caQTL/Runs/CD4_TCM/Results/TenK10K.cis_qtl_pairs.chr{chr_num}.csv"
    )
)

cd4_tcm_caqtl_raw[, position := as.integer(sub("^[^:]+:([^:]+):.*", "\\1", variant_id))]
cd4_tcm_caqtl_raw[, log_p := -log10(pval_nominal)]
cd4_tcm_caqtl_raw[, celltype := "CD4 TCM"]

cd14_mono_caqtl_raw <- fread(
    cmd = glue(
        "head -n 1 /g/data/ei56/od8037/NewGenotypes/caQTL/Runs/CD14_Mono/Results/TenK10K.cis_qtl_pairs.chr{chr_num}.csv && grep '{peak_name}' /g/data/ei56/od8037/NewGenotypes/caQTL/Runs/CD14_Mono/Results/TenK10K.cis_qtl_pairs.chr{chr_num}.csv"
    )
)

cd14_mono_caqtl_raw[, position := as.integer(sub("^[^:]+:([^:]+):.*", "\\1", variant_id))]
cd14_mono_caqtl_raw[, log_p := -log10(pval_nominal)]
cd14_mono_caqtl_raw[, celltype := "CD14 Mono"]

# Obtain nominal p-values
p_cd4_tcm <- cd4_tcm_caqtl_raw[
    phenotype_id == peak_name & variant_id == snp_cd4_tcm,
    signif(pval_nominal, 2)
]

p_cd14_mono <- cd14_mono_caqtl_raw[
    phenotype_id == peak_name & variant_id == snp_cd14_mono,
    signif(pval_nominal, 2)
]

format_p <- function(p) {
    sci <- formatC(p, format = "e", digits = 2)
    parts <- strsplit(sci, "e")[[1]]
    mantissa <- parts[1]
    exponent <- as.integer(parts[2])
    list(m = mantissa, e = exponent)
}

p_parts_cd4 <- format_p(p_cd4_tcm)
p_parts_cd14 <- format_p(p_cd14_mono)

# Generate violin plots
p1 <- ggplot(plot_df_cd4_tcm, aes(x = .data[[snp_cd4_tcm]], y = .data[[peak_name]])) +
    geom_violin(aes(group = .data[[snp_cd4_tcm]]), fill = "#d8d8d8", colour = "#c5c5c5") +
    geom_boxplot(aes(group = .data[[snp_cd4_tcm]]), width = 0.2, outlier.shape = NA, fill = "#c5c5c5", colour = "#c5c5c5") +
    geom_smooth(method = "lm", colour = "#c5c5c5", se = FALSE) +
    labs(
        title = bquote(atop("CD4"["TCM"], "(" * italic(p) == .(p_parts_cd4$m) * " × 10"^.(p_parts_cd4$e) * ")")),
        x = snp_cd4_tcm,
        y = peak_name
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 14)
    )

p2 <- ggplot(plot_df_cd14_mono, aes(x = .data[[snp_cd14_mono]], y = .data[[peak_name]])) +
    geom_violin(aes(group = .data[[snp_cd14_mono]]), fill = "#ffcd9a", colour = "#f7b874") +
    geom_boxplot(aes(group = .data[[snp_cd14_mono]]), width = 0.2, outlier.shape = NA, fill = "#f7b874", colour = "#f7b874") +
    geom_smooth(method = "lm", colour = "#f7b874", se = FALSE) +
    labs(
        title = bquote(atop("CD14"["Mono"], "(" * italic(p) == .(p_parts_cd14$m) * " × 10"^.(p_parts_cd14$e) * ")")),
        x = snp_cd14_mono,
        y = ""
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 14)
    )

# Locus and gene plot
library(AnnotationHub)
library(GenomeInfoDb)
library(Signac)

source("/g/data/ei56/od8037/Plotting/QTL_Boxplots/custom_gene_plot.R")

annotations <- readRDS("/g/data/ei56/od8037/NewGenotypes/Coloc/caQTL2eQTL/Plot_Examples/annotations.rds")
pbmc <- readRDS("/g/data/ei56/ax3061/proj/tenk10k/caQTL/data/QCed/tob_atac_S0220_1_QCed.rds")

peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]
Annotation(pbmc) <- annotations

plot_df <- rbind(
    cd14_mono_caqtl_raw[, .(position, log_p, celltype)],
    cd4_tcm_caqtl_raw[, .(position, log_p, celltype)]
) %>%
    dplyr::mutate(
        outline = case_when(
            celltype == "CD14 Mono" & position == as.integer(sub("^[^:]+:([^:]+):.*", "\\1", snp_cd14_mono)) ~ "black",
            celltype == "CD4 TCM" & position == as.integer(sub("^[^:]+:([^:]+):.*", "\\1", snp_cd4_tcm)) ~ "black",
            TRUE ~ "transparent"
        )
    )

start_bp <- min(plot_df$position) + (max(plot_df$position) - min(plot_df$position)) * 0.15
end_bp <- max(plot_df$position) - (max(plot_df$position) - min(plot_df$position)) * 0.15

match <- str_match(peak_name, "^chr(\\d+):(\\d+)-(\\d+)$")

p3 <- plot_df %>%
    ggplot(aes(x = position, y = log_p)) +
    geom_point(
        data = filter(plot_df, outline != "black"),
        aes(fill = celltype),
        shape = 21,
        size = 1.2,
        stroke = 0,
        colour = "transparent"
    ) +
    geom_point(
        data = filter(plot_df, outline == "black"),
        aes(fill = celltype),
        shape = 21,
        size = 1.8,
        stroke = 0.5,
        colour = "black"
    ) +
    scale_fill_manual(values = c("CD14 Mono" = "#f7b874", "CD4 TCM" = "#c5c5c5")) +
    annotate(
        "rect",
        xmin = as.integer(match[3]) + 3500,
        xmax = as.integer(match[4]) - 3500,
        ymin = -Inf,
        ymax = Inf,
        fill = "dodgerblue",
        alpha = 0.35
    ) +
    labs(
        x = glue("Position on Chromosome {chr_num}"),
        y = expression(-log["10"](italic(P)["caQTL"]))
    ) +
    xlim(start_bp, end_bp) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 14)
    )

region_string <- paste0("chr", chr_num, "-", start_bp, "-", end_bp)

gene_plot <- AnnotationPlot(object = pbmc, gapwidth = 70000, region = region_string) +
    scale_x_continuous(
        name = glue("Position on Chromosome {chr_num} (Mbp)"),
        breaks = seq(round(start_bp / 1e5) * 1e5, end_bp, by = 500000),
        labels = number_format(scale = 1e-6, accuracy = 0.01),
        limits = c(start_bp, end_bp)
    ) +
    labs(x = glue("Position on Chromosome {chr_num} (Mbp)")) +
    theme(
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 13)
    )

p4 <- p3 / gene_plot +
    plot_layout(heights = c(1, 0.75))

# Combine plots
p <- p1 + p2 + p4 +
    plot_layout(ncol = 3, widths = c(1, 1, 2.7))

# Save plot
out_file <- glue("/g/data/ei56/od8037/Plotting/QTL_Boxplots/Plots/{img_name}.png")
ggsave(
    filename = out_file,
    plot = p,
    width = 10,
    height = 4
)
