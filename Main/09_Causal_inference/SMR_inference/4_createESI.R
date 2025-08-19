
# Use this part when we want to parallel different chromosome's run

args = commandArgs(trailingOnly=TRUE)
cell_type <- args[1]
chr_given <- args[2]

# Upload libraries
library(data.table)
library(dplyr)
library(stringr)
setDTthreads(4)

# This is how Anne's output look like
# > output_df
# chromosome       SNP        genetic_distance       position      effect_allele     other_allele      freq
# <int>           <char>            <num>    <int>        <char>.    <char>     <num>
# 1:         10  10:45274190:C:T          0      45274190        T        C       0.5433400
# 2:         10  10:45274406:G:T          0      45274406        T        G       0.1054440
# 3:         10  10:45274474:G:A          0      45274474        A        G       0.5433400

work_dir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/"
# cell_type_list <- sort(list.dirs(paste0(work_dir, "output/main_analyses/preprocessing_caQTL_rawID"), full.names = FALSE, recursive = FALSE))

# for (cell_type in cell_type_list) {
besd_dir <- paste0(work_dir, "output/main_analyses_NewGenotypes/preprocessing_caQTL_rawID/", cell_type, "/BESD/")
for (chrNumber in seq(1,22)){
  if (chrNumber != chr_given) {
    next
  }
  caqtl_filepath <- paste0("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/Final_caQTL/NewGenotypes/Runs/", cell_type, "/Results/TenK10K.cis_qtl_pairs.chr", chrNumber, ".csv")
  
  # Plink file formatting: 21      21:5031125:CT:C 0       5031125 C       CT
  # Chr,Chr:POS:REF:ALT,Dist,POS,ALT,REF
  # Little ESI
  little_esi_filepath <- paste0(besd_dir, "Chr", chrNumber, ".esi")
  esi_df <- fread(little_esi_filepath)
  caqtl_df <- fread(caqtl_filepath)
  caqtl_df <- caqtl_df %>%
    filter(
      str_length(str_split_fixed(variant_id, ":", 4)[,3]) == 1 &
        str_length(str_split_fixed(variant_id, ":", 4)[,4]) == 1
    )
  caqtl_df <- caqtl_df[caqtl_df$pval_nominal < 0.05]
  caqtl_df <- caqtl_df[caqtl_df$af > 0.01]
  caqtl_df <- caqtl_df[, .(variant_id, af)]

  genotype_dir <- "/directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/"
  bim <- fread(paste0(genotype_dir, "TenK10K_TOB_ATAC_renamed_chr", chrNumber, "_common_variants_qced.bim"), col.names = c("chr", "snp_id", "cm", "pos", "allele1", "allele2"))

  # Match ESI
  # caqtl_df <- caqtl_df[variant_id %in% esi_df$V2]
  caqtl_df <- unique(caqtl_df, by = c("variant_id"))

  output_df <- merge(caqtl_df, bim, by.x = 'variant_id', by.y = 'snp_id', all.x = TRUE)
  ### Use original snp
  # output_df[, variant_id := paste(variant_id, allele2, allele1, sep=":")]
  ### Use hg38 position as snp id
  #output_df[, variant_id := paste(chrNumber, pos, allele2, allele1, sep=":")]

  colnames(output_df) <- c("SNP", "freq", "chromosome", "genetic_distance", "position", "effect_allele", "other_allele")
  output_df <- output_df[!duplicated(output_df$SNP), ]
  # Reorder columns for output
  output_df <- output_df[, .(chromosome, SNP, genetic_distance, position, effect_allele, other_allele, freq)]
  fwrite(output_df, little_esi_filepath, sep="\t", quote=F, col.names=F, na = "NA")
}
# }
