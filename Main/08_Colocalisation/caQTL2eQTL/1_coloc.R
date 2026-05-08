# Load libraries
library(coloc)
library(tidyverse)
library(glue)
library(data.table)

## 1. Command line arguments -------------------------------------------------------------
# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
celltype <- args[1]
print(paste("Processing celltype:", celltype))
chr_num <- args[2]
print(paste("Chromosome:", chr_num))

## 2.1. Preprocessing caQTL -------------------------------------------------------------
# Load raw caQTL data
print("Loading raw caQTL data ...")
caqtl_raw <- fread(glue(
    "/g/data/ei56/od8037/TenK10K/caQTLNew/Results/{celltype}/TenK10K.cis_qtl_pairs.chr{chr_num}.csv"
), select = c("phenotype_id", "variant_id", "af", "slope", "slope_se", "pval_nominal"))

# Remove INDELs
print("Removing INDELs ...")
caqtl_raw[, c("chr_num", "position", "ref", "effect") := tstrsplit(variant_id, ":", fixed = TRUE)]
caqtl_raw <- caqtl_raw[nchar(ref) == 1 & nchar(effect) == 1]

# Extract significant caQTL peaks
caqtl_threshold <- 5e-8
print(glue("Extracting significant caQTL peaks with p-value < {caqtl_threshold} ..."))
sig_result_df <- caqtl_raw[pval_nominal < caqtl_threshold]
sig_peaks <- unique(sig_result_df$phenotype_id)

print(paste0("Number of significant peaks for Chromosome ", chr_num, ": ", length(sig_peaks)))
if (length(sig_peaks) == 0) {
    print(paste("No significant peaks found for Chromosome", chr_num, "- Terminating script."))
    stop(paste("No significant peaks found for Chromosome", chr_num, "- Terminating script."))
}

# Load number of donors for each celltype
donor_df <- read.csv("/g/data/ei56/od8037/Final_Coloc/caQTL2GWAS/celltype_donors.csv") 
num_donors <- donor_df[donor_df$celltype == celltype, "num_donors"]


## Process caQTL data
# Apply efficient vectorised transformations in-place
print("Processing caQTL data ...")
caqtl_raw[, `:=`(
    beta = slope,  
    MAF = pmin(af, 1 - af),  
    varbeta = slope_se^2,  
    snp = paste0("chr", chr_num, "_", position, "_", ref, "_", effect), 
    N = num_donors
)]
caqtl_raw <- caqtl_raw[, .(phenotype_id, beta, varbeta, position, snp, N, MAF)]
##


## 2.2. Preprocessing eQTL -------------------------------------------------------------
# Load eQTL data
eqtl_df <- fread(glue(
    "/g/data/ei56/ax3061/proj/tenk10k/eQTL/{celltype}_common_all_cis_raw_pvalues.tsv"
))[CHR == chr_num]

# Process eQTL data
eqtl_df[, `:=`(
    beta = BETA,
    varbeta = SE^2,
    position = POS,
    snp = paste0("chr", gsub(":", "_", MarkerID)),
    MAF = pmin(AF_Allele2, 1 - AF_Allele2)
)]
eqtl_df <- eqtl_df[
    !is.na(beta),
    .(beta, varbeta, position, snp, N, MAF, gene)
]

# Add gene positions
gene_pos_df <- fread("/g/data/ei56/od8037/Final_Coloc/caQTL2eQTL/gene_pos.csv")[, .(gene_id, centre)]
eqtl_df <- merge(
    eqtl_df, 
    gene_pos_df, 
    by.x = "gene", 
    by.y = "gene_id", 
    all.x = TRUE
)
##


## 3. Colocalisation analysis -------------------------------------------------------------
# Set results directory
results_dir <- glue("/g/data/ei56/od8037/NewGenotypes/Coloc/caQTL2eQTL/Coloc_Results/{celltype}")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Set cis window size
cis_window_size <- 1000000
print(glue("cis window size: {cis_window_size}"))

# Loop through caQTL peaks and perform colocalisation analysis
result_df <- data.frame(
    peak = character(),
    gene = character(),
    nsnps_coloc_tested = integer(),
    PP.H0.abf = numeric(),
    PP.H1.abf = numeric(),
    PP.H2.abf = numeric(),
    PP.H3.abf = numeric(),
    PP.H4.abf = numeric(),
    celltype = character(),
    chrom = character(),
    top_snp = character()
)
for (peak in sig_peaks) {

    # Extract cis window coordinates
    peak_start <- as.numeric(strsplit(strsplit(peak, ":")[[1]][2], "-")[[1]][1]) 
    peak_end <- as.numeric(strsplit(strsplit(peak, ":")[[1]][2], "-")[[1]][2]) 

    cis_start <- peak_start + as.integer((peak_end - peak_start)/2) - cis_window_size
    cis_end <- peak_start + as.integer((peak_end - peak_start)/2) + 1 + cis_window_size

    # Filter eQTL data based on cis window
    eqtl_df_subset <- eqtl_df[position >= cis_start & position <= cis_end & centre >= cis_start & centre <= cis_end]

    if (nrow(eqtl_df_subset) == 0) {
        print(glue("No SNP eQTL data for {peak} in the cis-window: skipping ..."))
        next
    }
    print(glue("Extracted SNP eQTL data for {peak}!"))

    # Filter caQTL data based on chosen peak
    caqtl_data <- caqtl_raw[
        phenotype_id == peak & (beta != 0 | varbeta != 0) & !is.na(beta) & !is.na(varbeta) & !is.na(position)
    ]

    # Process caQTL data
    caqtl_data <- as.list(caqtl_data[, !"phenotype_id"])
    caqtl_data$type <- "quant"

    # Loop through each gene
    for (gene_name in unique(eqtl_df_subset$gene)) {
        print(glue("Processing gene: {gene_name} ..."))
        eqtl_data <- as.list(eqtl_df_subset[gene == gene_name, !"gene"])
        eqtl_data$type <- "quant"

        shared_snps <- length(intersect(eqtl_data$snp, caqtl_data$snp))
        if (shared_snps == 0) {
            print(glue("No common SNPs between eQTL and caQTL data for {peak}: skipping ..."))
            next
        }
        print(glue("Number of SNPs shared between eQTL and caQTL data for {peak}: {shared_snps}"))

        # Perform colocalisation analysis
        my.res <- coloc.abf(
            dataset1 = caqtl_data,
            dataset2 = eqtl_data
        )
        
        # Save results to file
        p_df <- data.frame(
            peak = peak,
            gene = gene_name,
            nsnps_coloc_tested = my.res$summary[1],
            PP.H0.abf = my.res$summary[2],
            PP.H1.abf = my.res$summary[3],
            PP.H2.abf = my.res$summary[4],
            PP.H3.abf = my.res$summary[5],
            PP.H4.abf = my.res$summary[6],
            celltype = celltype,
            chrom = paste0("chr", chr_num),
            top_snp = arrange(my.res$results, desc(SNP.PP.H4))[1,1],
            stringsAsFactors = FALSE
        )
        result_df <- rbind(result_df, p_df)
    }
}

# Print number of significant signals
print(glue("Number of significant signals: {sum(result_df$PP.H4.abf > 0.8)}"))

# Save results to CSV
write.csv(result_df, glue("{results_dir}/chr{chr_num}.csv"), row.names=FALSE)
print("Finished!")
