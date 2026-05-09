library(tidyverse)
library(data.table)

out_dir <- "/g/data/fy54/od8037/TenK10K/RarecaQTL/Plots"

plot_df <- fread("/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/Conditional/summary.csv") %>%
    mutate(
        rare_before = Sig_Before_And_After + Sig_Before_Not_Sig_After,
        rare_after = Sig_Before_And_After + Not_Sig_Before_Sig_After,
        rare_in_common_caPeaks = Sig_Before_And_After
    )

colour_df <- fread("/g/data/ei56/od8037/NewGenotypes/caQTL/colour_palette_table.tsv") 
celltype_order <- colour_df$wg2_scpred_prediction[
    colour_df$wg2_scpred_prediction %in% celltype_list
]

stack_df <- plot_df %>%
    bind_rows(
        tibble(
            Celltype = c("ILC", "CD8_Proliferating"),
            Neither_Sig = 0,
            Sig_Before_And_After = 0,
            Not_Sig_Before_Sig_After = 0,
            Sig_Before_Not_Sig_After = 0,
            Total = 0
        )
    ) %>%
    select(
        Celltype,
        Sig_Before_And_After,
        Not_Sig_Before_Sig_After,
        Sig_Before_Not_Sig_After
    ) %>%
    pivot_longer(
        cols = -Celltype,
        names_to = "category",
        values_to = "count"
    ) %>%
    mutate(
        category = recode(
            category,
            Sig_Before_And_After = "Retained after conditioning",
            Not_Sig_Before_Sig_After = "Significant after conditioning only",
            Sig_Before_Not_Sig_After = "Lost after conditioning"
        ),
        category = factor(
            category,
            levels = c(
                "Lost after conditioning",
                "Retained after conditioning",
                "Significant after conditioning only"
            )
        ),
        Celltype = factor(Celltype, levels = celltype_order)
    )

total_df <- stack_df %>%
    group_by(Celltype) %>%
    summarise(
        total_count = sum(count),
        .groups = "drop"
    )

p <- ggplot(stack_df, aes(x = Celltype, y = count, fill = category)) +
    geom_col(width = 0.8) +

    geom_text(
        data = total_df,
        aes(
            x = Celltype,
            y = total_count,
            label = total_count
        ),
        inherit.aes = FALSE,
        angle = 90,
        hjust = -0.14,
        vjust = 0.5,
        size = 5.5
    ) +

    scale_fill_manual(
        values = c(
            "Lost after conditioning" = "grey75",
            "Retained after conditioning" = "#2166AC",
            "Significant after conditioning only" = "#67A9CF"
        )
    ) +

    scale_y_continuous(
        labels = scales::comma,
        expand = expansion(mult = c(0, 0.22))
    ) +

    labs(
        x = "Cell Type",
        y = "Number of Rare caQTLs",
        fill = NULL
    ) +

    theme_classic() +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.8, "cm"),

        axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1,
            size = 20
        ),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_blank()
    ) +

    scale_x_discrete(
        labels = function(x) {
            sapply(x, function(label) {
                if (grepl("_", label)) {
                    parse(text = gsub("_", "[", label, fixed = TRUE) %>% paste0("]"))
                } else if (label == "Treg") {
                    parse(text = "T[reg]")
                } else {
                    label
                }
            })
        }
    )

ggsave(
    filename = file.path(out_dir, "bar_plot.png"),
    plot = p,
    width = 15.5,
    height = 8,
    dpi = 300
)
