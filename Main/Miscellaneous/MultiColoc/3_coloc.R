# Load libraries
library(coloc)
library(tidyverse)
library(glue)
library(data.table)
library(fst)
library(qvalue)
source("/g/data/fy54/od8037/TenK10K/MultiColoc/analyse_peak.R")

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
celltype <- args[1]
print(paste("Processing celltype:", celltype))
condition <- args[2]
print(paste("Condition:", condition))
tmpdir <- args[3]
print(paste("Temporary directory:", tmpdir))
chr_num <- args[4]
print(paste("Chromosome number:", chr_num))

window_num <- 1000000
print(glue("Window size in base pairs: {format(window_num, scientific = FALSE)}"))
min_p_gwas <- 1e-4
print(glue("Minimum GWAS p-value: {min_p_gwas}"))

# Set results directory and output file
results_dir <- glue("/g/data/fy54/od8037/TenK10K/MultiColoc/caQTL2GWAS/Results/{condition}/{celltype}")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- glue("{results_dir}/chr{chr_num}.csv")
if (file.exists(output_file)) {
    file.remove(output_file)
}

# Load number of donors for each celltype
donor_df <- fread("/g/data/fy54/od8037/TenK10K/MultiColoc/eQTL2GWAS/num_donors.csv")[celltype_name == celltype] 
num_donors <- donor_df$num_donors[1]

# Load caQTL summary data
sig_result_df <- fread(glue(
    "/g/data/ei56/od8037/NewGenotypes/caQTL/ProcessedResults/caQTLSummary/{celltype}.csv"
))[, .(phenotype_id)]

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

# Load caPeak locations for the current chromosome
capeak_loc_df <- sig_result_df[grepl(glue("^chr{chr_num}:"), phenotype_id)]
print(glue("Number of caQTLs for chromosome {chr_num}: {nrow(capeak_loc_df)}"))
if (nrow(capeak_loc_df) == 0) {
    message(glue("No caQTLs found for chromosome {chr_num} - Terminating script for this task."))
    quit(save = "no", status = 0)
}

# Load full GWAS and caQTL data
print("Loading GWAS and caQTL data ...")
gwas_df <- read_fst(glue(
    "/g/data/fy54/od8037/TenK10K/MultiColoc/caQTL2GWAS/GWAS/{condition}/chr{chr_num}.fst"
), as.data.table = TRUE)
setkey(gwas_df, position)

caqtl_df <- read_fst(glue(
    "/g/data/ei56/od8037/TenK10K/caQTLNew/ProcessedResults/{celltype}/chr{chr_num}.fst"
), as.data.table = TRUE)
setkey(caqtl_df, phenotype_id)

# Set genotype path
plink_bfile <- glue(
    "/g/data/ei56/ax3061/proj/tenk10k/caQTL/data/genotype/TenK10K_TOB_ATAC_renamed_chr{chr_num}_common_variants_qced"
)

for (i in seq_along(capeak_loc_df$phenotype_id)) {
    # Extract gene name
    peak_name <- capeak_loc_df$phenotype_id[i]

    # Extract cis window coordinates
    cis_start <- capeak_loc_df$cis_start[i]
    if (cis_start < 0) cis_start <- 0
    cis_end <- capeak_loc_df$cis_end[i]

    # Filter GWAS data based on cis window
    gwas_df_subset <- gwas_df[position %between% c(cis_start, cis_end)]
    if (nrow(gwas_df_subset) == 0) {
        print(glue("No SNP GWAS data for {peak_name}: skipping ..."))
        next
    } else if (min(gwas_df_subset$P, na.rm = TRUE) >= min_p_gwas) {
        print(glue("No significant GWAS hits for {peak_name}: skipping ..."))
        next
    }
    gwas_df_subset[, A1:= sub(".*:", "", snp)]

    # Filter caQTL data based on chosen peak
    caqtl_df_subset <- caqtl_df[J(peak_name)]

    # Find number of shared SNPs
    shared_snps <- intersect(unique(gwas_df_subset$snp), unique(caqtl_df_subset$snp))
    if (length(shared_snps) == 0) {
        print(glue("No common SNPs for {peak_name}: skipping ..."))
        next
    }

    # Generate LD matrix using PLINK
    shared_snps_df <- gwas_df_subset[snp %in% shared_snps] %>%
        dplyr::arrange(position)
    tmp_allele <- glue("{tmpdir}/allele_{i}.txt")
    tmp_extract <- glue("{tmpdir}/extract_{i}.txt")
    tmp_ld <- glue("{tmpdir}/ld_{i}")
    writeLines(paste(shared_snps_df$A1, shared_snps_df$snp, sep = "\t"), tmp_allele)
    writeLines(shared_snps_df$snp, tmp_extract)

    exit_code <- tryCatch(
        sys::exec_wait(
            "/g/data/ei56/od8037/software/PLINK/plink", 
            c(
                "--bfile",     plink_bfile,
                "--r",         "square",
                "--make-just-bim",
                "--chr",       chr_num,
                "--from-bp",   cis_start,
                "--to-bp",     cis_end,
                "--extract",   tmp_extract,
                "--silent",
                "--a1-allele", tmp_allele, "1", "2",
                "--threads",   1,
                "--out",       tmp_ld
            )
        ),
        error = function(e) 1L
    )
    unlink(c(tmp_extract, tmp_allele))

    if (exit_code != 0 || !file.exists(glue("{tmp_ld}.ld"))) {
        unlink(glue("{tmp_ld}*"))
        print(glue("LD matrix generation failed for {peak_name}: skipping ..."))
        next
    }

    bim_file <- fread(glue("{tmp_ld}.bim"))
    mat_ld <- as.matrix(fread(glue("{tmp_ld}.ld")))
    dimnames(mat_ld) <- list(bim_file$V2, bim_file$V2)
    unlink(glue("{tmp_ld}*"))

    mat_snps <- intersect(bim_file$V2, shared_snps)
    snps_na <- unique(rownames(which(is.na(mat_ld), arr.ind = TRUE)))
    final_snps <- setdiff(mat_snps, snps_na)
    mat_ld <- mat_ld[final_snps, final_snps]

    # Process GWAS and caQTL data
    gwas_data <- as.list(gwas_df_subset[snp %in% final_snps][order(match(snp, final_snps))])
    gwas_data$type <- "cc"
    gwas_data$N <- 398668
    gwas_data$LD <- mat_ld

    caqtl_data <- as.list(caqtl_df_subset[snp %in% final_snps][order(match(snp, final_snps))])
    caqtl_data$type <- "quant"
    caqtl_data$N <- num_donors
    caqtl_data$LD <- mat_ld

    # Perform SuSiE and colocalisation analysis
    result <- tryCatch({
        analyse_peak(peak_name, celltype, chr_num, gwas_data, caqtl_data)
    }, error = function(e) {
        print(glue("Error during SuSiE/coloc for {peak_name}: {e$message}"))
        return(NULL)
    })

    if (is.null(result)) next

    write_header <- !file.exists(output_file)
    fwrite(result, file = output_file, append = TRUE, col.names = write_header)

    print(glue("Coloc completed for {peak_name}!"))
}

if (!file.exists(output_file)) {
    message(glue("No coloc results were generated for chromosome {chr_num}."))
}

print("Finished!")
