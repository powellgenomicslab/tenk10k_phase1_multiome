
library(data.table)

pheno_names <- c(
  "alanine_aminotransferase", 
  "albumin", 
  "alkaline_phosphatase", 
  "apolipoprotein_a", 
  "apolipoprotein_b", 
  "aspartate_aminotransferase", 
  "c_reactive_protein", 
  "calcium", 
  "cholesterol", 
  "creatinine", 
  "cystatin_c", 
  "eosinophil_count", 
  "eosinophil_percent", 
  "gamma_glutamyltransferase", 
  "glucose", 
  "glycated_haemoglobin", 
  "haematocrit", 
  "haemoglobin_concentration", 
  "hdl_cholesterol", 
  "igf_1", 
  "ldl_cholesterol_direct", 
  "lymphocyte_count", 
  "lymphocyte_percent", 
  "mean_corpuscular_haemoglobin_concentration", 
  "mean_corpuscular_haemoglobin", 
  "mean_corpuscular_volume", 
  "mean_platelet_volume", 
  "mean_sphered_cell_volume", 
  "neutrophil_count", 
  "neutrophil_percent", 
  "phosphate", 
  "platelet_count", 
  "platelet_crit", 
  "platelet_distribution_width", 
  "red_blood_cell_count", 
  "red_blood_cell_distribution_width", 
  "shbg", 
  "total_bilirubin", 
  "total_protein", 
  "triglycerides", 
  "urate", 
  "urea", 
  "vitamin_d", 
  "white_blood_cell_count"
)

save_dir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/blood_traits/table_H4/"

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