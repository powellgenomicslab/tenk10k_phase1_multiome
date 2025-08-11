
### This part should be used when the bash file designates which cell type to run, and we're paralleling different chromosomes.
#!/g/data1a/ei56/as8574/micromamba/envs/r_matrixeqtl/bin Rscript
args = commandArgs(trailingOnly=TRUE)
cell_type <- args[1]
chr <- args[2]

library(data.table)
library(dplyr)
library(stringr)
setDTthreads(4)

# Directories
proj_dir <- paste0("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/Final_caQTL/NewGenotypes/Runs/", cell_type, "/Results/")

#for (cell_type in cell_type_list) {
# caqtl_dir <- paste0(proj_dir, cell_type)
caqtl_dir <- proj_dir

preprocessing_dir <- paste0("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/preprocessing_caQTL_rawID/", cell_type)
matrix_caqtl_dir <- paste0(preprocessing_dir, "/matrixQTL/")
dir.create(preprocessing_dir, showWarnings = FALSE)
dir.create(matrix_caqtl_dir, showWarnings = FALSE)

genotype_dir <- "/directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/"

caqtl_filepath <- paste0(caqtl_dir, "/TenK10K.cis_qtl_pairs.chr", chr, ".csv")
# if(!file.exists(caqtl_filepath)) {
#  next  # Skip to the next chromatin
# }
input_df <- fread(caqtl_filepath)
input_df <- input_df %>%
  filter(
    str_length(str_split_fixed(variant_id, ":", 4)[,3]) == 1 & 
      str_length(str_split_fixed(variant_id, ":", 4)[,4]) == 1
  )

input_df <- input_df[af > 0.01]
input_df[, Tstat := slope / slope_se]

# Need to reformat the SNPID 
output_df <- input_df[ , .(variant_id, phenotype_id, slope, Tstat, pval_nominal)]
output_df <- output_df[output_df$pval_nominal < 0.05]

### Do not need as the new caQTL using WGS already presents the correct SNP and genomic location
# Merge with genotype to get allele information
bim <- fread(paste0(genotype_dir, "TenK10K_TOB_ATAC_renamed_chr", chr, "_common_variants_qced.bim"), col.names = c("chr", "snp_id", "cm", "pos", "allele1", "allele2"))

output_df <- merge(output_df, bim, by.x = 'variant_id', by.y = 'snp_id', all.x = TRUE)
output_df <- output_df[!is.na(chr)]

### Use original snp
## output_df[, variant_id := paste(variant_id, allele2, allele1, sep=":")]
### Use hg38 position as snp id
#output_df[, variant_id := paste(chr, pos, allele2, allele1, sep=":")]

output_df <- output_df[, .(variant_id, phenotype_id, slope, Tstat, pval_nominal)]
output_df <- unique(output_df, by = c("variant_id", "phenotype_id"))
colnames(output_df) <- c("SNP", "peak", "beta", "t-stat", "p-value")

matrix_caqtl_filename <- paste0(matrix_caqtl_dir, "Chr", chr, "_MatrixcaQTL.tsv")
fwrite(output_df, matrix_caqtl_filename, sep = "\t", na = "NA", quote = FALSE)
#}