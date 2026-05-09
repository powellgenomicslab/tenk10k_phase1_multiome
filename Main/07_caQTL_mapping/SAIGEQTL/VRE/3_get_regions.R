library(tidyverse)
library(data.table)
library(glue)

celltype_list <- readLines("/g/data/fy54/od8037/Brenner/SAIGEQTL/Phenotypes/celltype_names.txt")
for (celltype in celltype_list) {
    print(glue("Processing {celltype}"))
    all_peaks <- readLines(glue(
        "/g/data/fy54/od8037/Brenner/SAIGEQTL/Phenotypes/ValidPeaks/{celltype}.txt"
    ))

    for (chr_num in 1:22) {
        print(glue("Processing chromosome {chr_num}"))
        peaks <- all_peaks[grepl(glue("^chr{chr_num}:"), all_peaks)]

        peak_df <- data.table(peak_id = peaks) %>%
            separate(peak_id, into = c("chr", "start", "end"), sep = "[:-]", remove = FALSE) %>%
            dplyr::mutate(
                start = as.numeric(start),
                end = as.numeric(end),
                centre_start = start + as.integer((end - start)/2),
                centre_end = centre_start + 1,
                chr_num = chr_num,
                cis_start = centre_start - 1000000,
                cis_end = centre_end + 1000000
            ) %>%
            dplyr::mutate(
                cis_start = ifelse(cis_start < 0, 0, cis_start),
                peak_id = gsub(":", "_", peak_id),
                peak_id = gsub("-", "_", peak_id)
            ) %>%
            dplyr::select(peak_id, chr_num, cis_start, cis_end)

        dir.create(glue("/g/data/fy54/od8037/Brenner/SAIGEQTL/Common/Regions/{celltype}"), showWarnings = FALSE)
        write.table(
            peak_df,
            file = glue("/g/data/fy54/od8037/Brenner/SAIGEQTL/Common/Regions/{celltype}/chr{chr_num}.csv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
    }
}
