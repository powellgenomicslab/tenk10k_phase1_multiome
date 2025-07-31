# Load libraries
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

# Define function to read CSV files and convert to GRanges
read_csv_to_gr <- function(file) {
    print(paste("Reading CSV:", file))
    peaks_df <- read.csv(file) %>%
        dplyr::select(chrom, start, end) %>%
        dplyr::rename(chr = chrom) %>%
        dplyr::filter(chr %in% paste0("chr", 1:22)) # Keep only autosomes
    makeGRangesFromDataFrame(peaks_df)
}

# Combine peak halves for CD14_Mono
cd14_gr_half1 <- read_csv_to_gr("/g/data/ei56/od8037/TenK10K/PeakCalling/BatchedPeakCalling/CelltypePeaks/CD14_Mono_half1.csv")
cd14_gr_half2 <- read_csv_to_gr("/g/data/ei56/od8037/TenK10K/PeakCalling/BatchedPeakCalling/CelltypePeaks/CD14_Mono_half2.csv")
cd14_combined_peaks <- reduce(c(cd14_gr_half1, cd14_gr_half2))

# Combine peak halves for CD4_TCM
cd4_gr_half1 <- read_csv_to_gr("/g/data/ei56/od8037/TenK10K/PeakCalling/BatchedPeakCalling/CelltypePeaks/CD4_TCM_half1.csv")
cd4_gr_half2 <- read_csv_to_gr("/g/data/ei56/od8037/TenK10K/PeakCalling/BatchedPeakCalling/CelltypePeaks/CD4_TCM_half2.csv")
cd4_combined_peaks <- reduce(c(cd4_gr_half1, cd4_gr_half2))

# Load full CSV-based peaks 
csv_files <- list.files(
    "/g/data/ei56/od8037/TenK10K/PeakCalling/CelltypePeaks", 
    pattern = "\\.csv$", 
    full.names = TRUE
)
gr_list <- lapply(csv_files, read_csv_to_gr)

# Add combined peaks to the list
full_gr_list <- c(gr_list, list(cd14_combined_peaks, cd4_combined_peaks))

# Merge peaks
combined_peaks <- reduce(x = do.call(c, full_gr_list))

# Filter out bad peaks based on length
peakwidths <- width(combined_peaks)
combined_peaks <- combined_peaks[peakwidths < 10000 & peakwidths > 20]

# Save combined peaks as BED file
export.bed(combined_peaks, con = "/g/data/ei56/od8037/TenK10K/PeakCalling/MACS3_Combined_Peaks.bed")