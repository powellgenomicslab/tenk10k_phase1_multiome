# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(patchwork)

# Load metadata
result_df <- data.frame(
    celltype = character(),
    count = integer(),
    library_name = character()
)
for (cur_library in readLines("/g/data/ei56/od8037/TenK10K/PeakCalling/library_names.txt")) {
    meta_df <- fread(glue(
        "/g/data/ei56/jf1058/TenK10K/Multiome/data/annotated/annotated/tob_atac_{cur_library}_annotated_meta_new_ref_nonPBMC_excluded.csv"
    ))[!is.na(donor_id)] %>% 
        group_by(predicted.id) %>%
        summarise(count = n(), .groups = "drop") %>%
        as.data.frame() %>%
        rename(celltype = predicted.id) %>%
        mutate(
            celltype = gsub("_", " ", celltype),
            library_name = cur_library
        )
    result_df <- rbind(result_df, meta_df)
}

# Load celltype colours
colour_df <- fread("/g/data/ei56/od8037/NewGenotypes/caQTL/colour_palette_table.tsv") 
celltype_order <- colour_df$wg2_scpred_prediction[
    colour_df$wg2_scpred_prediction %in% full_results$celltype
]

# Plot results
p1 <- result_df %>%
    group_by(celltype) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(celltype = gsub(" ", "_", celltype)) %>%
    left_join(colour_df, by = c("celltype" = "wg2_scpred_prediction")) %>%
    mutate(celltype = factor(celltype, levels = celltype_order)) %>%
    ggplot(aes(x = celltype, y = count, fill = color)) +
        geom_bar(stat = "identity") +
        scale_fill_identity() +  
        theme_classic(base_size = 20) +
        labs(
            x = "Cell Type",
            y = "Number of Cells"
        ) +
        theme(
            legend.position = "none",
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_text(size = 25)
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
        ) +
        scale_y_log10(
            breaks = 10^seq(0, 6, by = 2),
            labels = scales::comma_format(),
            expand = expansion(mult = c(0, 0.1))
        )

p2 <- result_df %>%
    mutate(celltype = gsub(" ", "_", celltype)) %>%
    left_join(colour_df, by = c("celltype" = "wg2_scpred_prediction")) %>%
    mutate(celltype = factor(celltype, levels = celltype_order)) %>%
    ggplot(aes(x = celltype, y = count, fill = color, colour = color)) +
        scale_fill_identity() +  
        scale_colour_identity() +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +  # violin as density
        geom_jitter(width = 0) + 
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        theme_classic(base_size = 20) +
        labs(x = "Cell Type", y = "Number of Cells") +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
            axis.title.x = element_text(size = 25),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_text(size = 25),
            legend.position = "none"
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

p <- p1 / p2 + plot_layout(heights = c(1, 1))
ggsave(
    p,
    filename = glue("/g/data/ei56/od8037/Plotting/Plots/celltype_counts_combined.png"),
    width = 15.5, height = 9.9
)
