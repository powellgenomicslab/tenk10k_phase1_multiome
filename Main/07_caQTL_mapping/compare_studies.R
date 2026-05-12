# Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(readxl)
library(rtracklayer)
library(stringr)
library(GenomicRanges)

## Replication study for Mu et al ----------------------------------------
# 1. Load and filter published results (hg19)
mu_results <- read_excel(
    "/g/data/ei56/od8037/Plotting/Compare_Studies/mu_et_al.xlsx",
    skip = 2,
    sheet = 10
) %>% 
    as.data.frame() %>%
    janitor::clean_names() 

# 2. Convert hg19 peaks to hg38 using liftOver
# Extract chr, start and end 
peak_parts <- str_match(mu_results$peak, "(chr[0-9XY]+)_(\\d+)-(\\d+)")
gr_hg19 <- GRanges(
    seqnames = peak_parts[, 2],
    ranges = IRanges(start = as.integer(peak_parts[, 3]), end = as.integer(peak_parts[, 4]))
)

# Import chain file and liftOver
chain <- import.chain("/g/data/ei56/od8037/Plotting/Compare_Studies/hg19ToHg38.over.chain")
lifted <- liftOver(gr_hg19, chain)

# Get indices where liftOver succeeded
valid_idx <- which(lengths(lifted) > 0)

# Extract components only for successful liftOvers
seqs <- vapply(lifted[valid_idx], function(x) as.character(seqnames(x)[1]), character(1))
starts <- vapply(lifted[valid_idx], function(x) start(x)[1], integer(1))
ends <- vapply(lifted[valid_idx], function(x) end(x)[1], integer(1))

# Construct GRanges without NAs
gr_hg38 <- GRanges(seqnames = seqs, ranges = IRanges(start = starts, end = ends))

# Add to dataframe
mu_results$peak_hg38 <- NA_character_
mu_results$peak_hg38[valid_idx] <- paste0(seqs, ":", starts, "-", ends)
mu_results <- mu_results[valid_idx, ]

# 3. Compare with our results
# Cell type mapping
recode_map <- list(
    B = c("B_naive", "B_memory", "B_intermediate", "Plasmablast"),
    CD4.T = c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD4_CTL", "CD4_Proliferating", "Treg"),
    CD8.T = c("CD8_Naive", "CD8_TCM", "CD8_TEM"),
    DC = c("cDC1", "cDC2", "pDC", "ASDC"),
    Mono = c("CD14_Mono", "CD16_Mono"),
    NK = c("NK", "NK_CD56bright", "NK_Proliferating"),
    other.T = c("gdT", "MAIT", "dnT")
)
mu_results$overlap <- "N"

# Get all unique mapped cell types from mu_results
mapped_cells <- unique(mu_results$ca_qtl_cell)
mapped_cells <- mapped_cells[mapped_cells %in% names(recode_map)]

# Loop through each unique celltype once
for (cell in mapped_cells) {
    print(cell)
    our_celltypes <- recode_map[[cell]]
    if (is.na(our_celltypes[1])) next  # skip unmapped
    
    # Subset rows in mu_results for this cell type
    idx <- which(mu_results$ca_qtl_cell == cell)
    mu_subset <- mu_results[idx, ]

    # Build GRanges of all peaks from mu_subset
    peak_parts <- str_match(mu_subset$peak_hg38, "(chr[0-9XY]+):(\\d+)-(\\d+)")
    mu_gr <- GRanges(
        seqnames = peak_parts[, 2],
        ranges = IRanges(start = as.integer(peak_parts[, 3]), end = as.integer(peak_parts[, 4]))
    )

    # Collect all peaks from mapped celltypes
    our_peaks <- list()
    for (our_ct in our_celltypes) {
        dt <- fread(glue(
            "/g/data/ei56/od8037/NewGenotypes/caQTL/ProcessedResults/caQTLSummary/{our_ct}.csv"
        ))[qval < 0.05]
        
        # Parse phenotype_id into GRanges
        our_peaks[[our_ct]] <- GRanges(
            seqnames = tstrsplit(dt$phenotype_id, ":", fixed = TRUE)[[1]],
            ranges = IRanges(
                start = as.integer(tstrsplit(tstrsplit(dt$phenotype_id, ":", fixed = TRUE)[[2]], "-", fixed = TRUE)[[1]]),
                end = as.integer(tstrsplit(tstrsplit(dt$phenotype_id, ":", fixed = TRUE)[[2]], "-", fixed = TRUE)[[2]])
            )
        )
    }

    # Combine GRanges from all mapped celltypes
    if (length(our_peaks) == 0) next
    our_combined_gr <- do.call(c, unname(our_peaks))

    # Check overlaps and update mu_results
    overlaps <- countOverlaps(mu_gr, our_combined_gr) > 0
    mu_results$overlap[idx[overlaps]] <- "Y"
}

# 4. Plot results
# Count number of replicated vs non-replicated caQTLs per cell type
mu_plot_df <- mu_results %>%
    mutate(status = ifelse(overlap == "Y", "Replicated", "Not Replicated")) %>%
    group_by(ca_qtl_cell, status) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(ca_qtl_cell) %>%
    mutate(total = sum(n)) %>%
    ungroup()

## Replication study for Benaglio et al ----------------------------------------
# 1. Load and filter published results (hg19)
benaglio_results <- read_excel(
    "/g/data/ei56/od8037/Plotting/Compare_Studies/benaglio_et_al.xlsx",
    skip = 1
) %>% 
    as.data.frame() %>%
    janitor::clean_names() 

# 2. Convert hg19 peaks to hg38 using liftOver
# Extract chr, start and end 
peak_parts <- str_match(benaglio_results$feature_id, "(chr[0-9XY]+):(\\d+)-(\\d+)")
gr_hg19 <- GRanges(
    seqnames = peak_parts[, 2],
    ranges = IRanges(start = as.integer(peak_parts[, 3]), end = as.integer(peak_parts[, 4]))
)

# Import chain file and liftOver
chain <- import.chain("/g/data/ei56/od8037/Plotting/Compare_Studies/hg19ToHg38.over.chain")
lifted <- liftOver(gr_hg19, chain)

# Get indices where liftOver succeeded
valid_idx <- which(lengths(lifted) > 0)

# Extract components only for successful liftOvers
seqs <- vapply(lifted[valid_idx], function(x) as.character(seqnames(x)[1]), character(1))
starts <- vapply(lifted[valid_idx], function(x) start(x)[1], integer(1))
ends <- vapply(lifted[valid_idx], function(x) end(x)[1], integer(1))

# Construct GRanges without NAs
gr_hg38 <- GRanges(seqnames = seqs, ranges = IRanges(start = starts, end = ends))

# Add to dataframe
benaglio_results$peak_hg38 <- NA_character_
benaglio_results$peak_hg38[valid_idx] <- paste0(seqs, ":", starts, "-", ends)
benaglio_results <- benaglio_results[valid_idx, ]

# 3. Compare with our results
# Cell type mapping
recode_map <- list(
    act_cd4_t     = NA,
    adaptive_NK   = "NK",
    b             = NA,
    bulk          = NA,
    cDC           = c("cDC1", "cDC2"),
    cMono         = "CD14_Mono",
    cyto_cd8_t    = NA,
    cyto_nk       = "NK",
    iMono         = "CD14_Mono",
    mem_b         = "B_memory",
    mem_cd8_t     = c("CD8_TCM", "CD8_TEM"),
    mkc           = NA,
    mono          = NA,
    naive_b       = "B_naive",
    naive_cd4_t   = "CD4_Naive",
    naive_cd8_t   = "CD8_Naive",
    ncMono        = "CD16_Mono",
    nk            = NA,
    t             = NA,
    tReg          = "Treg"
)
benaglio_results$overlap <- "N"

# Get all unique mapped cell types from benaglio_results
mapped_cells <- unique(benaglio_results$cell_type)
mapped_cells <- mapped_cells[mapped_cells %in% names(recode_map)]

# Loop through each unique benaglio celltype once
for (cell in mapped_cells) {
    print(cell)
    our_celltypes <- recode_map[[cell]]
    if (is.na(our_celltypes[1])) next

    # Subset rows in benaglio_results for this cell type
    idx <- which(benaglio_results$cell_type == cell)
    benaglio_subset <- benaglio_results[idx, ]

    # Build GRanges of all peaks from benaglio_subset
    peak_parts <- str_match(benaglio_subset$peak_hg38, "(chr[0-9XY]+):(\\d+)-(\\d+)")
    benaglio_gr <- GRanges(
        seqnames = peak_parts[, 2],
        ranges = IRanges(start = as.integer(peak_parts[, 3]), end = as.integer(peak_parts[, 4]))
    )

    # Collect all peaks from mapped celltypes
    our_peaks <- list()
    for (our_ct in our_celltypes) {
        dt <- fread(glue(
            "/g/data/ei56/od8037/NewGenotypes/caQTL/ProcessedResults/caQTLSummary/{our_ct}.csv"
        ))[qval < 0.05]
        
        # Parse phenotype_id into GRanges
        our_peaks[[our_ct]] <- GRanges(
            seqnames = tstrsplit(dt$phenotype_id, ":", fixed = TRUE)[[1]],
            ranges = IRanges(
                start = as.integer(tstrsplit(tstrsplit(dt$phenotype_id, ":", fixed = TRUE)[[2]], "-", fixed = TRUE)[[1]]),
                end = as.integer(tstrsplit(tstrsplit(dt$phenotype_id, ":", fixed = TRUE)[[2]], "-", fixed = TRUE)[[2]])
            )
        )
    }

    # Combine GRanges from all mapped celltypes
    if (length(our_peaks) == 0) next
    our_combined_gr <- do.call(c, unname(our_peaks))

    # Check overlaps and update benaglio_results
    overlaps <- countOverlaps(benaglio_gr, our_combined_gr) > 0
    benaglio_results$overlap[idx[overlaps]] <- "Y"
}

# 4. Plot results
# Count number of replicated vs non-replicated caQTLs per cell type
to_exclude <- c("act_cd4_t", "cyto_cd8_t", "mkc", "nk", "mono", "t", "b", "bulk")
benaglio_plot_df <- benaglio_results %>%
    mutate(status = ifelse(overlap == "Y", "Replicated", "Not Replicated")) %>%
    filter(!(cell_type %in% to_exclude)) %>%
    group_by(cell_type, status) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cell_type) %>%
    mutate(total = sum(n)) %>%
    ungroup()

## Plot combined results ----------------------------------------
# Combine datasets with a source column
benaglio_df <- benaglio_plot_df %>%
    mutate(source = "Benaglio et al.")
mu_df <- mu_plot_df %>%
    dplyr::rename(cell_type = ca_qtl_cell) %>%
    mutate(source = "Mu et al.")
combined_df <- bind_rows(benaglio_df, mu_df) %>%
    mutate(
        cell_type = gsub("_", " ", cell_type),
        cell_type = gsub("\\.", " ", cell_type)
    )

# Calculate percent replicated for line plot
replication_pct <- combined_df %>%
    group_by(source, cell_type) %>%
    summarise(
        total = sum(n),
        replicated = sum(n[status == "Replicated"]),
        .groups = "drop"
    ) %>%
    mutate(percent = 100 * replicated / total)

# Plot combined results
max_n <- 3000
p <- ggplot(combined_df, aes(x = fct_reorder(cell_type, -total), y = n, fill = status)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    geom_line(
        data = replication_pct,
        aes(x = cell_type, y = percent * max_n / 100, group = 1),
        inherit.aes = FALSE,
        color = "#558dbb"
    ) +
    geom_point(
        data = replication_pct,
        aes(x = cell_type, y = percent * max_n / 100),
        inherit.aes = FALSE,
        color = "#558dbb"
    ) +
    facet_wrap(~source, scales = "free_x") +
    scale_y_continuous(
        name = "Number of caQTLs",
        limits = c(0, max_n),
        sec.axis = sec_axis(~ . / max_n * 100, name = "Percent of Replication", breaks = seq(0, 100, 25)),
        expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(values = c("Replicated" = "#d22818", "Not Replicated" = "#ffc9d1")) +
    theme_classic() +
    labs(x = "Cell Type", y = NULL, fill = NULL) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        axis.title.y.right = element_text(color = "#558dbb"),
        axis.text.y.right = element_text(color = "#558dbb"),
        axis.ticks.y.right = element_line(color = "#558dbb"),
        axis.line.y.right = element_line(color = "#558dbb")
    )

ggsave(p, filename = "/g/data/ei56/od8037/Plotting/Compare_Studies/combined_plot.png", width = 8, height = 4)
