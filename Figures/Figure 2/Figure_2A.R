# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(patchwork)

## caQTL Summary
# Collect caQTL results
result_df <- data.frame(
    celltype = character(),
    chr_num = integer(),
    num_caqtl = integer()
)
for (celltype in readLines("/g/data/ei56/od8037/Final_caQTL/celltype_names.txt")) {
    for (chr_num in 1:22) {
        num_caqtl <- fread(glue(
            "/g/data/ei56/od8037/NewGenotypes/caQTL/Runs/{celltype}/Results/TenK10K.sig_cis_qtl_pairs.chr{chr_num}.csv"
        )) %>% filter(qval < 0.05) %>% nrow()
        result_df <- rbind(result_df, data.frame(
            celltype = celltype,
            chr_num = chr_num,
            num_caqtl = num_caqtl
        ))
    }
}
full_results <- result_df %>% 
    group_by(celltype) %>% 
    summarise(num_caqtl = sum(num_caqtl), .groups = "drop") %>%
    as.data.frame() %>%
    rbind(data.frame(
        celltype = c("CD8_Proliferating", "ILC"),
        num_caqtl = c(0, 0)
    ))

# Plot results
colour_df <- fread("/g/data/ei56/od8037/NewGenotypes/caQTL/colour_palette_table.tsv") 
celltype_order <- colour_df$wg2_scpred_prediction[
    colour_df$wg2_scpred_prediction %in% full_results$celltype
]
plot_df <- full_results %>%
    left_join(colour_df, by = c("celltype" = "wg2_scpred_prediction")) %>%
    mutate(celltype = factor(celltype, levels = celltype_order)) 

p <- plot_df %>%
    ggplot(aes(x = celltype, y = num_caqtl, fill = color)) +
        geom_bar(stat = "identity") + 
        geom_text(
            aes(label = num_caqtl),
            angle = 90,
            hjust = -0.14,
            vjust = 0.5,
            size = 5.5
        ) +
        scale_fill_identity() +        
        theme_classic(base_size = 18) +              
        labs(
            x = "Cell Type",
            y = "Number of caQTLs"
        ) +
        theme(
            legend.position = "none",
            axis.text.x = element_text(
                angle = 90, vjust = 0.5, hjust = 1, size = 20
            ),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 25),
            axis.title.y = element_text(size = 25),
            plot.title = element_blank()
        ) +
        scale_y_continuous(
            expand = expansion(mult = c(0, 0.22))
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
ggplot2::ggsave(p, filename = glue("/g/data/ei56/od8037/Plotting/Plots/caqtl_by_celltype.png"), width = 15.5, height = 7)
