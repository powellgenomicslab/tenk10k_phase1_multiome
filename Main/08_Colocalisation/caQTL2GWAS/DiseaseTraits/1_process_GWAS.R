# Load libraries
library(data.table)
library(glue)
library(fst)

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
condition_index <- as.integer(args[1])

# Define condition names
pheno_names <- c(
    "alzheimer_GCST90027158",
    "breastca_GCST004988",
    "covid_GCST011071",
    "ibd_liu2023",
    "NHL_GCST90011819",
    "lungca_GCST004748",
    "lymphoma_GCST90018878",
    "parkinson_GCST009325",
    "prostateca_GCST90274713",
    "ra_GCST90132223",
    "sle_GCST003156",
    "myeloproliferative_GCST90000032",
    "lymphocytic_leukemia_GCST90011814",
    "nephrotic_GCST90258619",
    "kiryluk_IgAN"
)

# Define corresponding local filenames
pheno_files <- c(
    "GCST90027158.h_parsed.tsv",
    "GCST004988.h_parsed.tsv",
    "GCST011071_parsed.tsv",
    "ibd_EAS_EUR_SiKJEF_meta_IBD.tsv",
    "NHL_GCST90011819_parsed.tsv",
    "GCST004748.h_parsed.tsv",
    "GCST90018878.h_parsed.tsv",
    "GCST009325.h_parsed.tsv",
    "GCST90274713.h_parsed.tsv",
    "GCST90132223_parsed.tsv",
    "bentham_2015_26502338_sle_parsed.tsv",
    "myeloproliferative_GCST90000032_parsed.tsv",
    "lymphocytic_leukemia_GCST90011814_parsed.tsv",
    "nephrotic_GCST90258619_parsed.tsv",
    "Kiryluk_IgAN_parsed.tsv"
)
pheno_files <- file.path("/g/data/ei56/jf1058/TenK10K/Multiome/GWAS", pheno_files)

# Construct full paths
pheno_map <- setNames(pheno_names, pheno_files)

# Load SNP mapping file
print("Loading SNP mapping file ...")
snp_map <- read_fst("/g/data/fy54/od8037/TenK10K/ColocNew/snp_id_map.fst", as.data.table = TRUE)
setkey(snp_map, swapped_style)

# Get the condition name from the index
condition_path <- pheno_files[condition_index]
condition_name <- pheno_map[condition_path]
print(condition_name)

# Load GWAS data for condition
gwas_df <- fread(condition_path)[grepl("_[ACGT]_[ACGT]$", snp) & varbeta != 0]

# Re-map SNP IDs in GWAS data to PLINK format
print("Re-mapping SNP IDs in GWAS data to PLINK format ...")
setnames(gwas_df, "snp", "old_snp")
gwas_df[, snp := gsub("_", ":", sub("^chr", "", old_snp))]
gwas_df[snp_map, bim_style := i.bim_style, on = .(snp = swapped_style)]
gwas_df[!is.na(bim_style), `:=`(
    beta = -beta,
    snp = bim_style
)]
gwas_df <- na.omit(gwas_df[, .(chromosome, beta, varbeta, position, snp, p_value)])[!duplicated(snp)]

# Save the processed GWAS data for each chromosome
print("Saving GWAS data ...")
for (chr_num in 1:22) {
    print(glue("Chromosome {chr_num}"))
    dir_path <- glue("/g/data/fy54/od8037/TenK10K/ColocNew/DiseaseTraits/GWAS/{condition_name}")
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
    }
    gwas_df_subset <- gwas_df[chromosome == glue("chr{chr_num}"), .(beta, varbeta, position, snp)]
    write_fst(gwas_df_subset, glue("{dir_path}/chr{chr_num}.fst"))
}

print("Finished!")