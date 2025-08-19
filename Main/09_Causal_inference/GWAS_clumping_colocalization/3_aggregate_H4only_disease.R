
library(data.table)

pheno_names <- c(
  "alzheimer_GCST90027158",
  "breastca_GCST004988",
  "covid_GCST011071",
  "lungca_GCST004748",
  "lymphoma_GCST90018878",
  "parkinson_GCST009325",
  "prostateca_GCST90274713",
  "ra_GCST90132223",
  "sle_GCST003156",
  "kiryluk_IgAN",
  "asthma",
  "ms",
  "t1dm",
  "CD_EAS_EUR",
  "IBD_EAS_EUR",
  "UC_EAS_EUR"
)

save_dir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/disease_traits/table_H4_newGWAS/"

coloc_table <- list()
smr_table <- list()
row_count <- 1
for (trait in pheno_names) {
  coloc_df <- fread(paste0(save_dir, "coloc_", trait, ".csv"))
  coloc_df <- coloc_df[V1 == trait]
  coloc_table[[row_count]] <- coloc_df
  smr_df <- fread(paste0(save_dir, "smr_", trait, ".csv"))
  smr_df <- smr_df[V1 == trait]
  smr_table[[row_count]] <- smr_df
  row_count <- row_count + 1
}

coloc_table <- rbindlist(coloc_table)
smr_table <- rbindlist(smr_table)

write.csv(coloc_table[, -"V1"], paste0(save_dir, "coloc_table.csv"), row.names = coloc_table$V1)
write.csv(smr_table[, -"V1"], paste0(save_dir, "smr_table.csv"), row.names = smr_table$V1)

dir.create(paste0(save_dir, "trait_specific/"), showWarnings = FALSE)

for (trait in pheno_names) {
  file.rename(paste0(save_dir, "coloc_", trait, ".csv"), paste0(save_dir, "trait_specific/coloc_", trait, ".csv"))
  file.rename(paste0(save_dir, "smr_", trait, ".csv"), paste0(save_dir, "trait_specific/smr_", trait, ".csv"))
}