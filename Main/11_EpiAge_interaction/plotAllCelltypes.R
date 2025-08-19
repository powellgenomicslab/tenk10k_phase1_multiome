library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)

# Load the data
csv_files <- list.files(path = "/directflow/SCCGGroupShare/projects/petall/ipfLung-atac/tenk10k_atac/epitrace_output",
                        pattern = "epitrace_meta_tob_atac_S.*\\.csv", full.names = TRUE)

epitrace_list <- list()

for(i in seq(length(csv_files))) {
    df <- read.csv(csv_files[i])
    name <- sub(".*epitrace_meta_tob_atac_(S.*)\\.csv", "\\1", basename(csv_files[i]))
    df$cellName <- paste0(name, "_", df$X)
    rownames(df) <- df$cellName
    epitrace_list[[i]] <- df
    names(epitrace_list)[i] <- name
}

epitrace_full <- rbindlist(epitrace_list)

# plot epigenetic age by cell type

color_palatte <- read_delim("/directflow/SCCGGroupShare/projects/petall/tenk10k_epiage_tensorQTL/data/colour_palette_table.tsv", trim_ws = TRUE)

epitrace_full <- epitrace_full %>%
    dplyr::filter(!(predicted.id %in% c("Eryth", "Platelet", "CD4 Proliferating")))

epitrace_full$predicted.id <- gsub("NK_CD56bright", "NK CD56bright", epitrace_full$predicted.id)

epitrace_full <- merge(epitrace_full, color_palatte, by.x = "predicted.id", by.y = "cell_type", all.x = TRUE)

# Calculate the median EpiTraceAge_iterative for each cell type
celltype_medians <- epitrace_full %>%
    dplyr::group_by(predicted.id) %>%
    dplyr::summarise(median_age = median(EpiTraceAge_iterative))

# Order the cell types by median age
celltype_order <- celltype_medians %>%
    dplyr::arrange(median_age) %>%
    dplyr::pull(predicted.id)

epitrace_full$predicted.id <- factor(epitrace_full$predicted.id, levels = celltype_order, ordered = TRUE)

# Define custom labels with subscripts
celltype_labels <- c(
    "CD14 Mono" = expression(CD14[Mono]),
    "CD16 Mono" = expression(CD16[Mono]),
    "CD4 Naive" = expression(CD4[Naive]),
    "CD8 Naive" = expression(CD8[Naive]),
    "B naive" = expression(B[naive]),
    "B intermediate" = expression(B[intermediate]),
    "B memory" = expression(B[memory]),
    "CD4 TCM" = expression(CD4[TCM]),
    "CD8 TCM" = expression(CD8[TCM]),
    "CD8 Proliferating" = expression(CD8[Proliferating]),
    "NK CD56bright" = expression(NK[CD56bright]),
    "CD4 TEM" = expression(CD4[TEM]),
    "CD8 TEM" = expression(CD8[TEM]),
    "CD4 CTL" = expression(CD4[CTL])
)

# In your ggplot, add scale_x_discrete(labels = ...) after labs()
p1 <- ggplot(epitrace_full, aes(x = predicted.id, y = EpiTraceAge_iterative, fill = predicted.id)) +
        geom_violin(scale = "width", alpha = 0.8) +
        geom_boxplot(width = 0.1, alpha = 0.5) +
        scale_fill_manual(values = setNames(epitrace_full$color, epitrace_full$predicted.id)) +
        theme_classic(base_size = 20) +
        labs(title = paste(""),
                 x = "",
                 y = "Relative Epigenetic Age",
                 fill = "Cell Type") +
        scale_x_discrete(labels = celltype_labels) +
        theme(
                legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
                axis.text.y = element_text(size = 16),
                axis.title.x = element_text(size = 20),
                axis.title.y = element_text(size = 18)
        )

ggsave(filename = paste0("figures/fig5a-tenk10k-EpiTraceAge-Celltype.png"), plot = p1, width = 16, height = 5, units = "in")