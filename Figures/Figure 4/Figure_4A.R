library(tidyverse)
library(data.table)
library(glue)
library(ComplexUpset)

celltype_list <- readLines("/g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/celltype_names.txt")
colour_df <- fread("/g/data/ei56/od8037/NewGenotypes/caQTL/colour_palette_table.tsv") 
celltype_order <- colour_df$wg2_scpred_prediction[
    colour_df$wg2_scpred_prediction %in% celltype_list
]

pg_list <- list()
for (celltype in celltype_list) {
    print(glue("Processing {celltype}..."))

    df_list <- vector("list", 22)
    for (chr_num in 1:22) {
        df_list[[chr_num]] <- fread(glue(
            "/g/data/ei56/od8037/NewGenotypes/Coloc/caQTL2eQTL/Coloc_Results/{celltype}/chr{chr_num}.csv"
        ))[PP.H4.abf > 0.8]
    }

    df <- rbindlist(df_list, fill = TRUE)

    if (nrow(df) > 0) {
        pg_list[[celltype]] <- unique(paste0(df$peak, "_", df$gene))
    } else {
        pg_list[[celltype]] <- character(0)
    }
}

all_pairs <- unique(unlist(pg_list))
upset_df <- data.table(peak_gene = all_pairs)
for (celltype in names(pg_list)) {
    upset_df[, (celltype) := peak_gene %in% pg_list[[celltype]]]
}

celltype_cols <- colour_df$color
names(celltype_cols) <- colour_df$wg2_scpred_prediction

p <- upset(
    upset_df,
    intersect = celltype_order,
    n_intersections = 25,
    sort_sets = FALSE,
    width_ratio = 0.2,
    height_ratio = 0.65,
    name = "Unique or Shared Peak-Gene Pairs",
    queries = list(
        upset_query(set = "pDC", fill = "#795447FF"),
        upset_query(set = "cDC2", fill = "#5D3F37FF"),
        upset_query(set = "CD16_Mono", fill = "#FFAB91FF"),
        upset_query(set = "CD14_Mono", fill = "#D84314FF"),
        upset_query(set = "B_memory", fill = "#FFEB3AFF"),
        upset_query(set = "B_intermediate", fill = "#FABF2CFF"),
        upset_query(set = "B_naive", fill = "#F47F17FF"),
        upset_query(set = "NK", fill = "#2D7D32FF"),
        upset_query(set = "CD8_Naive", fill = "#1E87E5FF"),
        upset_query(set = "CD8_TEM", fill = "#0C46A0FF"),
        upset_query(set = "Treg", fill = "#B29DDAFF"),
        upset_query(set = "CD4_CTL", fill = "#9474CCFF"),
        upset_query(set = "CD4_Naive", fill = "#512CA7FF"),
        upset_query(set = "CD4_TCM", fill = "#311A92FF")
    ),
    set_sizes = (
        upset_set_size() +
        labs(y = "Peak-Gene Pairs per Cell Type")
    ),
    base_annotations = list(
        "Intersection size" = (
            intersection_size(
                text = list(
                    angle = 90, vjust = 0.5, hjust = -0.08
                )
            ) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                labs(y = "Peak-Gene Pairs")
        )
    )
) 

ggsave(
    p,
    filename = "/g/data/fy54/od8037/TenK10K/caQTL2eQTLPlot/plot.png",
    width = 14,
    height = 8,
    dpi = 600
)
