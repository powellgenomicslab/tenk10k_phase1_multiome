# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Load libraries
library(tidyverse)
library(data.table)
library(R.utils)
library(BSgenome)
library(GenomicRanges)
library(EDASeq)
library(glue)

# Get celltype for current run
celltype_list <- readLines("celltype_names.txt")
celltype <- celltype_list[as.numeric(args[1])]
print(paste("Processing celltype:", celltype))

# Load the genome reference
print("Loading genome reference")
genome <- getBSgenome("hg38")

# Read pseudobulk matrix and convert to a data frame
file_path <- glue("PseudobulkMatrices/{celltype}.csv")
if (file.exists(file_path)) {
    print("Reading pseudobulk matrix")
    mat <- fread(file_path)
} else {
    stop(glue("File {file_path} does not exist. Stopping script."))
}

mat <- as.data.frame(mat) %>%
    drop_na()
print(glue("Peaks remaining after filtering: {nrow(mat)}"))

# Create directory for runs of this celltype
run_dir <- glue("Runs/{celltype}")
dir.create(run_dir, showWarnings = FALSE)

# Extract peak identifiers and set them as rownames
peaks <- mat$V1
rownames(mat) <- peaks
mat <- dplyr::select(mat, -V1) %>% as.matrix()

# Create a data frame with peak location details
print("Obtaining peak location details")
cnt_table <- data.frame(loc = peaks) %>%
    tidyr::separate(loc, into = c("chr", "range"), sep = ":", remove = FALSE) %>% 
    tidyr::separate(range, into = c("start", "end"), sep = "-", remove = TRUE) %>% 
    dplyr::mutate(start = as.numeric(start), end = as.numeric(end))  

# Create a GRanges object for the peaks
print("Creating GRanges object")
gr <- GRanges(
    seqnames = cnt_table$chr, 
    ranges = IRanges(cnt_table$start, cnt_table$end), 
    strand = "*", 
    mcols = data.frame(peakID = cnt_table$loc)
)

# Calculate GC content for each peak
print("Calculating GC content")
peakSeqs <- getSeq(x = genome, gr)
gcContentPeaks <- letterFrequency(peakSeqs, "GC", as.prob = TRUE)[,1]

# Normalise read counts based on GC content
print("Normalising read counts")
data <- withinLaneNormalization(x = mat, y = gcContentPeaks, num.bins = 20, which = "full")

# Convert normalised data to a data frame
mx <- as.data.frame(data)

# Save the normalised matrix to a file
print("Saving normalised matrix")
mx_save <- data.frame(mx)
mx_save$loc <- rownames(mx_save)
fwrite(mx_save, glue("{run_dir}/PseudobulkMatNorm.csv"))

print("Finished!")
