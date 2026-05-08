
library(ggplot2)
library(scales)
library(gridExtra)
library(ggrastr)
library(MASS)
library(data.table)

peak <- 'chr16:50326863-50327752'
snp <- '16:50336792:G:C'
processed_df <- fread("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/plots/chr16_50326863_50327752_16_50336792_G_C_df_20bin.csv")
pseudotime_step <- 0.05
processed_df[, pseudotime := (as.numeric(sub("Bin", "", pseudotime_bin)) - 1) * pseudotime_step + runif(.N, 0, pseudotime_step)]

df_plot <- data.frame(expression = processed_df$chr16_50326863_50327752, cellfunction = processed_df$pseudotime, genotype = processed_df$genotype, genotype_label = processed_df$genotype_label)

# Create the plot
p <- ggplot(df_plot, aes(x = cellfunction, y = expression, color = as.factor(genotype))) +
    geom_point_rast(alpha = 0.01, size = 0.5, raster.dpi = 300) +
    geom_density_2d(bins = 4, alpha = 0.6) +  
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, fullrange = TRUE, linetype = "dashed") +
    scale_color_manual(
    values = c("0" = "#5cc5ef", "1" = "#ffb745", "2" = "#e7552c"),
    labels = setNames(df_plot$genotype_label, df_plot$genotype)  # maps 0/1/2 -> label
    ) +
    scale_x_continuous(limits = c(-0.5, 2)) +
    theme_classic() +
    labs(x = "Cell State", y = "ATAC Peak Accessibility", color = "Genotype",
         title = paste("Peak:", peak, "\nSNP:", snp))

ggsave("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/plots/plot.png", plot = p, width = 6, height = 4.5, dpi = 300)
