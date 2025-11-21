# Load libraries
library(coloc)
library(tidyverse)
library(glue)
library(data.table)
library(fst)

## 1. Command line arguments -------------------------------------------------------------
# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
celltype <- args[1]
print(paste("Processing celltype:", celltype))
condition <- args[2]
print(paste("Condition:", condition))

window_num <- 1000000
print(glue("Window size in base pairs: {format(window_num, scientific = FALSE)}"))

# Set results directory
results_dir <- glue("/g/data/ei56/od8037/NewGenotypes/Coloc/BloodTraits/Coloc_Results/{condition}/{celltype}")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Load number of donors for each celltype
donor_df <- read.csv("/g/data/ei56/od8037/Final_Coloc/caQTL2GWAS/celltype_donors.csv") 
num_donors <- donor_df[donor_df$celltype == celltype, "num_donors"]

for (chr_num in 1:22) {
    ## 2. Preprocessing -------------------------------------------------------------
    # Load caQTL data
    sig_result_df <- fread(glue(
        "/g/data/ei56/od8037/NewGenotypes/caQTL/ProcessedResults/caQTLSummary/{celltype}.csv"
    ))[grepl(glue("^chr{chr_num}:"), phenotype_id), .(phenotype_id)]

    print(paste0("Number of caQTLs for Chromosome ", chr_num, ": ", nrow(sig_result_df)))
    if (nrow(sig_result_df) == 0) {
        print(paste("No caQTLs found for Chromosome", chr_num, "- Skipping"))
        next
    }

    # Define cis coordinates
    peak_coords <- tstrsplit(sig_result_df$phenotype_id, ":", fixed = TRUE)
    peak_range <- tstrsplit(peak_coords[[2]], "-", fixed = TRUE)
    sig_result_df[, peak_start := as.integer(peak_range[[1]])]
    sig_result_df[, peak_end := as.integer(peak_range[[2]])]
    sig_result_df[, peak_mid := peak_start + as.integer((peak_end - peak_start)/2)]
    sig_result_df[, `:=`(
        cis_start = peak_mid - window_num,
        cis_end = peak_mid + window_num + 1
    )]

    # Load processed caQTL data
    print("Loading caQTL data ...")
    caqtl_result_dir <- glue(
        "/g/data/ei56/od8037/NewGenotypes/caQTL/ProcessedResults/{celltype}/chr{chr_num}.fst"
    )
    caqtl_df <- read_fst(caqtl_result_dir, as.data.table = TRUE)


    ## 3. Colocalisation analysis -------------------------------------------------------------
    # Load full GWAS data
    print("Loading GWAS data ...")
    gwas_path <- glue("/g/data/ei56/od8037/NewGenotypes/Coloc/BloodTraits/GWAS_Results/{condition}/chr{chr_num}.fst")
    gwas_df <- read_fst(gwas_path, as.data.table = TRUE)

    # Set keys for faster subsetting
    setkey(gwas_df, position)
    setkey(caqtl_df, phenotype_id)

    # Loop through caQTL peaks and perform colocalisation analysis
    result_list <- vector("list", length(sig_result_df$phenotype_id))
    for (i in seq_along(sig_result_df$phenotype_id)) {
        peak <- sig_result_df$phenotype_id[i]

        # Extract cis window coordinates
        cis_start <- sig_result_df$cis_start[i]
        cis_end <- sig_result_df$cis_end[i]

        # Filter GWAS data based on cis window
        gwas_df_subset <- gwas_df[position %between% c(cis_start, cis_end)]

        if (nrow(gwas_df_subset) == 0) {
            print(glue("No SNP GWAS data for {peak} in the cis-window: skipping ..."))
            next
        }
        print(glue("Extracted SNP GWAS data for {peak}!"))

        # Process GWAS data
        gwas_data <- as.list(gwas_df_subset)
        gwas_data$type <- "cc"

        # Filter caQTL data based on chosen peak
        caqtl_data_subset <- caqtl_df[J(peak)]

        # Process caQTL data
        caqtl_data <- as.list(caqtl_data_subset[, !"phenotype_id"])
        caqtl_data$type <- "quant"

        # Find number of shared SNPs
        shared_snps <- length(intersect(gwas_data$snp, caqtl_data$snp))
        if (shared_snps == 0) {
            print(glue("No common SNPs between GWAS and caQTL data for {peak}: skipping ..."))
            next
        }
        print(glue("Number of SNPs shared between GWAS and caQTL data for {peak}: {shared_snps}"))

        # Perform colocalisation analysis
        my.res <- coloc.abf(
            dataset1 = caqtl_data,
            dataset2 = gwas_data
        )
        
        # Save results to file
        p_df <- data.frame(
            peak = peak,
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
        result_list[[i]] <- p_df
    }

    result_df <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
    fwrite(result_df, glue("{results_dir}/chr{chr_num}.csv"), row.names=FALSE)
}

print("Finished!")
