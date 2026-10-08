library(data.table)
library(dplyr)
library(ggplot2)
library(ggdist)

xx <- fread("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/plots/chr16_50326863_50327752_16_50336792_G_C_df_3bin.csv")

donor_means <- copy(xx)
setnames(donor_means, "pseudotime_bin",          "pseudotime_group")
setnames(donor_means, "chr16_50326863_50327752", "value")

ns <- donor_means[, .(N = uniqueN(donor_id)), by = .(genotype_label, pseudotime_group)]

slopes <- donor_means[, {
  m <- lm(value ~ genotype)
  list(slope = sprintf("%.3f", coef(m)[2]))
}, by = pseudotime_group]

slope_labels <- slopes[, setNames(slope, pseudotime_group)]

# ── 3. Factor ordering ────────────────────────────────────────────────────────
donor_means[, genotype_label  := factor(genotype_label,  levels = c("G/G","G/C","C/C"))]
donor_means[, pseudotime_group := factor(pseudotime_group, levels = c("Early","Mid","Late"))]

# Regression lines summary (group centroid x = genotype int)
reg_lines <- donor_means[, .(
  x    = mean(genotype),
  yhat = mean(value)
), by = .(pseudotime_group, genotype_label, genotype)]

# ── 4. Palette & theme ────────────────────────────────────────────────────────
pal <- c(
  "Early" = "#5B9BD5",   # cool blue
  "Mid"   = "#70B77E",   # sage green
  "Late"  = "#E8A838"    # amber
)

theme_caQTL <- function() {
  theme_classic(base_family = "sans") +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),

      panel.grid.major.y = element_line(color = "#EBEBEB", linewidth = 0.4),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),

      axis.text  = element_text(color = "#444444", size = 24, family = "sans"),
      axis.title = element_text(color = "#222222", size = 26, family = "sans"),
      axis.line  = element_line(color = "#CCCCCC", linewidth = 0.4),
      axis.ticks = element_line(color = "#CCCCCC", linewidth = 0.4),

      plot.title    = element_text(
        family = "sans", size = 26, color = "#111111",
        face = "bold", lineheight = 1.3, margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        family = "sans", size = 24, color = "#555555", margin = margin(b = 6)
      ),
      plot.caption = element_text(color = "#888888", size = 20, family = "sans",
                                  hjust = 0, margin = margin(t = 10)),
      plot.margin  = margin(10, 24, 16, 20),

      legend.background = element_rect(fill = "white", color = "#DDDDDD", linewidth = 0.4),
      legend.key        = element_rect(fill = NA),
      legend.title      = element_text(color = "#222222", size = 24, family = "sans", face = "bold"),
      legend.text       = element_text(color = "#444444", size = 22, family = "sans"),
      legend.margin     = margin(6, 10, 6, 10),

      strip.text = element_text(color = "#222222", family = "sans", size = 18)
    )
}

legend_labels <- c(
  "Early" = paste0("Early (slope=", 1.107, ")"),
  "Mid"   = paste0("Mid (slope=",   0.630,   ")"),
  "Late"  = paste0("Late (slope=",  0.357,  ")")
)

dodge_offsets <- c("Early" = -0.85/3, "Mid" = 0, "Late" = 0.85/3)

reg_segments <- donor_means[, {
  m <- lm(value ~ genotype)
  list(
    y_start = predict(m, newdata = data.frame(genotype = 0)),
    y_end   = predict(m, newdata = data.frame(genotype = 2))
  )
}, by = pseudotime_group]

reg_segments[, x_start := 1 + dodge_offsets[pseudotime_group]]
reg_segments[, x_end   := 3 + dodge_offsets[pseudotime_group]]
reg_segments[, pseudotime_group := factor(pseudotime_group, levels = c("Early","Mid","Late"))]

p <- ggplot(donor_means,
            aes(x      = genotype_label,
                y      = value,
                fill   = pseudotime_group,
                color  = pseudotime_group)) +

  geom_violin(
    aes(group = interaction(genotype_label, pseudotime_group)),
    position     = position_dodge(width = 0.85),
    width        = 0.8,
    alpha        = 0.25,
    linewidth    = 0.3,
    trim         = FALSE,
    scale        = "width"
  ) +

  geom_boxplot(
    aes(group = interaction(genotype_label, pseudotime_group)),
    position      = position_dodge(width = 0.85),
    width         = 0.12,
    alpha         = 0.70,
    linewidth     = 0.55,
    outlier.shape = NA,
    coef          = 0            # whiskers = IQR only
  ) +

  stat_summary(
    aes(group = interaction(genotype_label, pseudotime_group)),
    fun         = median,
    geom        = "point",
    shape       = 21,
    size        = 0.2,
    fill        = "white",
    color       = "white",
    stroke      = 0.8,
    position    = position_dodge(width = 0.85)
  ) +

  geom_segment(
    data = reg_segments,
    aes(x     = x_start,
        xend  = x_end,
        y     = y_start,
        yend  = y_end,
        color = pseudotime_group),
    linewidth   = 0.45,
    alpha       = 0.9,
    inherit.aes = FALSE
  ) +

  ## ── Sample-size labels ────────────────────────────────────────────────────
  geom_text(
    data = ns,
    aes(x     = genotype_label,
        y     = 3.9,
        label = N,
        group = pseudotime_group),
    position  = position_dodge(width = 0.85),
    size      = 8,
    family    = "sans",
    color     = "#888888",
    fontface  = "plain",
    inherit.aes = FALSE
  ) +

  ## ── Scales ────────────────────────────────────────────────────────────────
  scale_fill_manual(values = pal, name = "Pseudotime") +
  scale_color_manual(values = pal, name = "Pseudotime") +
  scale_y_continuous(
    breaks = seq(-3, 3, 1),
    limits = c(-3, 4.1),
    expand = c(0, 0)
  ) +

  ## ── Labels ────────────────────────────────────────────────────────────────
  labs(
      title    = "Peak: chr16:50326863-50327752 | SNP: 16:50336792:G:C",
      subtitle = "p_interaction = 6.76e-28",
      x        = "Genotype",
      y        = "Donor mean log-CPM accessibility"
    ) +

  ## ── Slope annotation: colored text per group, no legend boxes ────────────
  annotate("text", x = Inf, y = Inf,
           label = "Pseudotime",
           hjust = 1.05, vjust = 40.5,
           size = 6.5, family = "sans", fontface = "bold", color = "#333333") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Early  slope = ", slope_labels["Early"]),
           hjust = 1.05, vjust = 45.5,
           size = 6.2, family = "sans", color = pal["Early"]) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Mid    slope = ", slope_labels["Mid"]),
           hjust = 1.05, vjust = 47.5,
           size = 6.2, family = "sans", color = pal["Mid"]) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Late   slope = ", slope_labels["Late"]),
           hjust = 1.05, vjust = 49.5,
           size = 6.2, family = "sans", color = pal["Late"]) +

  ## ── Theme ─────────────────────────────────────────────────────────────────
  theme_caQTL() +
  theme(legend.position = "none")

ggsave(
  "/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/plots/violin_plot3.pdf",
  plot   = p,
  width  = 16,
  height = 11,
  units  = "in",
  device = cairo_pdf,
  bg     = "white"
)

message("✓ Saved: violin_caQTL_beautiful.pdf / .png")




