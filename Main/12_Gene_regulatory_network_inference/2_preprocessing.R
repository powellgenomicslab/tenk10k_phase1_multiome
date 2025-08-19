# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Load libraries
library(glue)
library(tidyverse)
library(data.table)
library(R.utils)

# Get celltype for current run
celltype <- as.character(args[1])
print(glue("Processing celltype: {celltype}"))


## Load normalised pseudobulk matrix for repeat 1
print("Loading normalised pseudobulk matrix for repeat 1")
repeat_1 <- fread(glue("/g/data/ei56/od8037/TenK10K/caQTL_Full/Repeat_1/Runs/{celltype}/PseudobulkMatNorm.csv"))
colnames(repeat_1) <- sub("^X", "", colnames(repeat_1))
repeat_1 <- as.data.frame(repeat_1) %>%
    column_to_rownames("loc")
# Store original peak order from repeat 1 before filtering
original_order <- rownames(repeat_1)
# Obtain peaks with < 95% zeros
repeat_1_peaks <- rownames(repeat_1)[rowSums(repeat_1 == 0) / ncol(repeat_1) < 0.95]
print(glue("Number of peaks in repeat 1: {length(repeat_1_peaks)}"))

## Load normalised pseudobulk matrix for repeat 2
print("Loading normalised pseudobulk matrix for repeat 2")
repeat_2 <- fread(glue("/g/data/ei56/od8037/TenK10K/caQTL_Full/Repeat_2/Runs/{celltype}/PseudobulkMatNorm.csv"))
colnames(repeat_2) <- sub("^X", "", colnames(repeat_2))
repeat_2 <- as.data.frame(repeat_2) %>%
    column_to_rownames("loc")
# Obtain peaks with < 95% zeros
repeat_2_peaks <- rownames(repeat_2)[rowSums(repeat_2 == 0) / ncol(repeat_2) < 0.95]
print(glue("Number of peaks in repeat 2: {length(repeat_2_peaks)}"))


# Obtain union of peaks
print("Obtaining union of peaks")
peaks_union_unsorted <- union(repeat_1_peaks, repeat_2_peaks)
peaks_union <- original_order[original_order %in% peaks_union_unsorted]
print(glue("Number of peaks in union: {length(peaks_union)}"))


## Process repeat 1
print("Processing repeat 1")
repeat_1 <- repeat_1[peaks_union, ]
print("Calculating CPM")
lib_size <- colSums(repeat_1)
repeat_1 <- sweep(repeat_1, 2, lib_size, FUN = "/")
repeat_1 <- repeat_1 * 1e6
print("Standardising counts")
repeat_1 <- t(scale(t(repeat_1)))

## Process repeat 2
print("Processing repeat 2")
repeat_2 <- repeat_2[peaks_union, ]
print("Calculating CPM")
lib_size <- colSums(repeat_2)
repeat_2 <- sweep(repeat_2, 2, lib_size, FUN = "/")
repeat_2 <- repeat_2 * 1e6
print("Standardising counts")
repeat_2 <- t(scale(t(repeat_2)))


# Create directory for runs of this celltype
run_dir <- glue("/g/data/ei56/od8037/Final_caQTL/Runs/{celltype}")
dir.create(run_dir, showWarnings = FALSE)

# Create folder to store TensorQTL logs
dir.create(glue("{run_dir}/Logs"), showWarnings = FALSE)

# Find intersecting and unique donors
print("Finding intersecting and unique donors")
donors_intersect <- intersect(colnames(repeat_1), colnames(repeat_2))
donors_only_1 <- setdiff(colnames(repeat_1), donors_intersect)
donors_only_2 <- setdiff(colnames(repeat_2), donors_intersect)

# Average intersecting donors
print("Averaging intersecting donors")
mx_intersect <- (repeat_1[, donors_intersect] + repeat_2[, donors_intersect]) / 2

# Combine all columns: averaged + unique from each repeat
print("Combining all donor matrices")
mx <- cbind(
    mx_intersect,
    repeat_1[, donors_only_1, drop = FALSE],
    repeat_2[, donors_only_2, drop = FALSE]
)

# Re-standardise the averaged matrix
print("Re-standardising counts")
mx <- as.data.frame(t(scale(t(mx)))) %>%
    tidyr::drop_na() # Drop rows with NAs (ie. rows with 0 variance)

# Save the standardised matrix to a file
print("Saving standardised matrix")
mx_loc <- data.frame(mx)
colnames(mx_loc) <- sub("^X", "", colnames(mx_loc))
mx_loc$loc <- rownames(mx_loc)
fwrite(mx_loc, glue("{run_dir}/PseudobulkMatStandardised.csv"))

# Add back location information
print("Generating location information")
mx_loc <- mx_loc %>%
    tidyr::separate(loc, into = c("#chr", "range"), sep = ":", remove = TRUE) %>%  
    tidyr::separate(range, into = c("start", "end"), sep = "-", remove = TRUE) %>%
    dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) 

# Process each chromosome
print("Processing chromosomes")
for (chr_num in 1:22) {
    print(paste("Chromosome number:", chr_num))
    mx_subset <- mx_loc %>%
        dplyr::filter(`#chr` == paste0("chr", chr_num))

    # Add the TSS position
    mx_subset$start <- mx_subset$start + as.integer((mx_subset$end - mx_subset$start)/2)
    mx_subset$end <- mx_subset$start + 1
    mx_subset$phenotype_id <- rownames(mx_subset)

    # Rearrange columns
    data <- mx_subset %>%
        dplyr::select(`#chr`, start, end, phenotype_id, everything())

    # Save expression matrix
    dir.create(glue("{run_dir}/ExpressionBeds/"), showWarnings = FALSE)
    write.table(
        data, 
        glue("{run_dir}/ExpressionBeds/chr{chr_num}.bed"),
        row.names = F, col.names = T, quote = F, sep = "\t"
    )
    gzip(
        glue("{run_dir}/ExpressionBeds/chr{chr_num}.bed"), 
        destname = glue("{run_dir}/ExpressionBeds/chr{chr_num}.bed.gz"),
        overwrite = T
    )
}

# Perform PCA
print("Performing PCA")
expr <- t(as.matrix(mx))
prcompResult <- prcomp(expr, center = TRUE, scale. = TRUE) 
PCs <- prcompResult$x 

# Elbow plot
resultRunElbow <- PCAForQTL::runElbow(prcompResult = prcompResult)
print(paste("Recommended PCs:", resultRunElbow))

# Visualise result
print("Saving scree plot")
p <- PCAForQTL::makeScreePlot(prcompResult, labels = c("Elbow") ,values = c(resultRunElbow), titleText = "Scree Plot")
ggsave(p, filename = glue("{run_dir}/ScreePlot.png"))

# Generate data frame mapping each donor to their covariates
print("Mapping donor covariate data")
cov_dt <- fread("/g/data/ei56/od8037/TenK10K/caQTL_Combined/CD4_NC_covar_peer_factors_PF0.txt")
cov_dt <- t(cov_dt[1:8, ]) 
colnames(cov_dt) <- cov_dt["id", ]
cov_dt <- cov_dt[-1, ] %>% as.data.frame()
cov_dt[] <- lapply(cov_dt, as.numeric)
knownCovariates <- cov_dt[rownames(expr), ]
colnames(knownCovariates)[2:7] <- paste0("genotype_pc", 1:6) # This is to avoid potential confusion

# Select top PCs
print("Selecting top PCs")
PCsTop <- PCs[,1:resultRunElbow]

# Combine PCs and covariates
print("Combining PCs and covariates")
PCsTop <- scale(PCsTop) 
covariatesToUse <- cbind(knownCovariates, PCsTop) %>% 
    t() %>%
    as.data.frame()
covariatesToUse$id <- rownames(covariatesToUse)
covariatesToUse <- covariatesToUse[, c("id", setdiff(names(covariatesToUse), "id"))]

# Save covariates
print("Saving covariates")
fwrite(covariatesToUse, glue("{run_dir}/covariates.txt"))

print("Finished!")
