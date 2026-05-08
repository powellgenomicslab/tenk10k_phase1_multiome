# Load libraries
library(data.table)
library(susieR)
library(glue)
library(fst)
library(pbmcapply)
setDTthreads(1)

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
celltype_name <- args[1]
chr_num <- as.integer(args[2])
tmpdir <- args[3]
print(glue("Temp directory: {tmpdir}"))

# Load number of donors
donor_df <- fread("/g/data/ei56/od8037/Final_Coloc/caQTL2GWAS/celltype_donors.csv") 
num_donors <- as.integer(donor_df[donor_df$celltype == celltype_name, "num_donors"])

# Load caQTL results
print("Loading caQTL results...")
caqtl_df <- read_fst(glue(
    "/g/data/fy54/od8037/TenK10K/FineMapping/caQTLResults/1Mb/{celltype_name}/chr{chr_num}.fst"
), as.data.table = TRUE)
setkey(caqtl_df, phenotype_id)

# Set genotype path
plink_bfile <- glue(
    "/g/data/ei56/ax3061/proj/tenk10k/caQTL/data/genotype/TenK10K_TOB_ATAC_renamed_chr{chr_num}_common_variants_qced"
)

# Define function to process each peak
process_peak <- function(peak_name) {
    print(glue("{Sys.time()}: {peak_name}"))
    tryCatch({
        setDTthreads(1) 
        
        # Extract SNPs for this peak
        peak_df <- caqtl_df[J(peak_name)]
        snps <- peak_df$variant_id

        # Generate LD matrix using PLINK
        peak_save_name <- gsub("-", "_", gsub(":", "_", peak_name))
        snp_file <- glue("{tmpdir}/snps_{peak_save_name}.txt")
        out_prefix <- glue("{tmpdir}/plink_{peak_save_name}")
        fwrite(data.table(snps), snp_file, col.names = FALSE)

        exit_code <- system2(
            command = "/g/data/ei56/od8037/software/PLINK/plink",
            args = c(
                "--bfile", plink_bfile,
                "--extract", snp_file,
                "--r", "square",
                "--out", out_prefix,
                "--threads", "1",
                "--keep-allele-order",
                "--silent",
                "--write-snplist" 
            ),
            stdout = NULL,
            stderr = NULL
        )

        # Check for PLINK failure
        if (exit_code != 0 || !file.exists(glue("{out_prefix}.ld"))) {
            print(glue("PLINK failed for {peak_name} with exit code {exit_code}. Skipping."))
            unlink(c(snp_file, glue("{out_prefix}*")))
            return(NULL) 
        }

        ld_mat <- as.matrix(fread(glue("{out_prefix}.ld"), header = FALSE))
        snp_id_list <- readLines(glue("{out_prefix}.snplist"))
        dimnames(ld_mat) <- list(snp_id_list, snp_id_list)
        unlink(c(snp_file, glue("{out_prefix}*")))

        # Remove SNPs with NA in LD matrix
        keep_snps <- rownames(ld_mat)[complete.cases(ld_mat)]
        
        # Guard against empty LD matrix
        if(length(keep_snps) == 0) return(NULL)
        ld_mat <- ld_mat[keep_snps, keep_snps]

        # Extract Z-scores
        snp_order <- match(keep_snps, peak_df$variant_id)
        z <- peak_df$Z[snp_order]

        # Run SuSiE
        susie_fit <- susie_rss(z = z, R = ld_mat, L = 10, n = num_donors)

        # Create full PIP table
        pip_df <- data.table(
            celltype = celltype_name,
            peak = peak_name,
            SNP = keep_snps,
            PIP = susie_fit$pip
        )

        # Add credible set membership
        cs_info <- susie_fit$sets
        if (!is.null(cs_info)) {
            cs_membership <- rep(NA_integer_, length(keep_snps))
            cs_membership[unlist(cs_info$cs)] <- rep(seq_along(cs_info$cs), lengths(cs_info$cs))
            pip_df$Credible_Set <- cs_membership
        } else {
            pip_df$Credible_Set <- NA_integer_
        }
        rm(ld_mat, z, susie_fit)
        
        return(pip_df[!is.na(Credible_Set)])
        
    }, error = function(e) {
        message(glue("Error processing {peak_name}: {e$message}"))
        return(NULL)
    })
}

# Run with progress bar
peak_list <- unique(caqtl_df$phenotype_id)
print(glue("Processing {length(peak_list)} peaks..."))
res_list <- pbmclapply(
    peak_list,
    process_peak,
    mc.cores = 16
)

valid_results <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, res_list)
print(glue("Successfully processed {length(valid_results)} out of {length(peak_list)} peaks."))

# Combine results into dataframes
print("Saving credible sets...")
credible_save_dir <- glue("/g/data/fy54/od8037/TenK10K/FineMapping/SuSiEResults/{celltype_name}")

if (!dir.exists(credible_save_dir)) {
    dir.create(credible_save_dir, recursive = TRUE)
}

if (length(valid_results) > 0) {
    credible_table <- rbindlist(valid_results, use.names = TRUE, fill = TRUE)
    fwrite(credible_table, glue("{credible_save_dir}/chr{chr_num}.csv"))
    print("Finished!")
} else {
    print("No results to save.")
}
