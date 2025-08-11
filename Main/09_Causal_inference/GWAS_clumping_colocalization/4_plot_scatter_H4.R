
library(data.table)
library(ggplot2)
library(dplyr)

pheno_abbs <- c(
  "alzheimer",
  "breastca",
  "colorectalca",
  "covid",
  "ibd",
  "NHL",
  "lungca",
  "lymphoma",
  "parkinson",
  "prostateca",
  "ra",
  "sle",
  "myeloproliferative",
  "lymphocytic_leukemia",
  "nephrotic",
  "kiryluk_IgAN",
  "asthma",
  "ms",
  "t1dm",
  "crohns",
  "ibd",
  "uc",
  "alanine_amino", 
  "albumin", 
  "alkaline_pho", 
  "apolipoprotein_a", 
  "apolipoprotein_b", 
  "aspartate_amino", 
  "c_reactive_protein", 
  "calcium", 
  "cholesterol", 
  "creatinine", 
  "cystatin_c", 
  "eosinophil_count", 
  "eosinophil_percent", 
  "gamma_glutamyl", 
  "glucose", 
  "glycated_haem", 
  "haematocrit", 
  "haemoglobin_conc", 
  "hdl_cholesterol", 
  "igf_1", 
  "ldl_cholesterol_direct", 
  "lymphocyte_count", 
  "lymphocyte_percent", 
  "mean_corpus_haem_conc", 
  "mean_corpus_haem", 
  "mean_corpus_volume", 
  "mean_platelet_volume", 
  "mean_sphered_cell_volume", 
  "neutrophil_count", 
  "neutrophil_percent", 
  "phosphate", 
  "platelet_count", 
  "platelet_crit", 
  "platelet_dist_width", 
  "red_blood_cell_count", 
  "red_blood_cell_dist_width", 
  "shbg", 
  "total_bilirubin", 
  "total_protein", 
  "triglycerides", 
  "urate", 
  "urea", 
  "vitamin_d", 
  "white_blood_cell_count"
)


pheno_names <- c(
  "alzheimer_GCST90027158",
  "breastca_GCST004988",
  "colorectalca_GCST90129505",
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
  "kiryluk_IgAN",
  "asthma",
  "ms",
  "t1dm",
  "CD_EAS_EUR",
  "IBD_EAS_EUR",
  "UC_EAS_EUR",
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
mapping <- setNames(pheno_abbs, pheno_names)
###################################
coloc_table_disease <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/disease_traits/table_H4_newGWAS/coloc_table.csv")
coloc_table_disease <- coloc_table_disease[!V1 %in% "nephrotic_GCST90258619"]
coloc_table_disease$eQTL_prop <- (coloc_table_disease$Union_eQTL / coloc_table_disease$Num_GWAS_loci) * 100
coloc_table_disease$caQTL_prop <- (coloc_table_disease$Union_caQTL / coloc_table_disease$Num_GWAS_loci) * 100

smr_table_disease <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/disease_traits/table_H4_newGWAS/smr_table.csv")
smr_table_disease <- smr_table_disease[!V1 %in% "nephrotic_GCST90258619"]
smr_table_disease$eQTL_prop <- (smr_table_disease$Union_eQTL / smr_table_disease$Num_GWAS_loci) * 100
smr_table_disease$caQTL_prop <- (smr_table_disease$Union_caQTL / smr_table_disease$Num_GWAS_loci) * 100

###################################
coloc_table_blood <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/blood_traits/table_H4/coloc_table.csv")
coloc_table_blood$eQTL_prop <- (coloc_table_blood$Union_eQTL / coloc_table_blood$Num_GWAS_loci) * 100
coloc_table_blood$caQTL_prop <- (coloc_table_blood$Union_caQTL / coloc_table_blood$Num_GWAS_loci) * 100

smr_table_blood <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/blood_traits/table_H4/smr_table.csv")
smr_table_blood$eQTL_prop <- (smr_table_blood$Union_eQTL / smr_table_blood$Num_GWAS_loci) * 100
smr_table_blood$caQTL_prop <- (smr_table_blood$Union_caQTL / smr_table_blood$Num_GWAS_loci) * 100


disease_names <- c("alzheimer_GCST90027158", "breastca_GCST004988", "colorectalca_GCST90129505", "covid_GCST011071", "ibd_liu2023", 
                   "NHL_GCST90011819", "lungca_GCST004748", "lymphoma_GCST90018878", "parkinson_GCST009325", "prostateca_GCST90274713", 
                   "ra_GCST90132223", "sle_GCST003156", "myeloproliferative_GCST90000032", "lymphocytic_leukemia_GCST90011814", 
                   "nephrotic_GCST90258619", "kiryluk_IgAN", "asthma", "ms", "t1dm", "CD_EAS_EUR", "IBD_EAS_EUR", "UC_EAS_EUR")
  
blood_count_names <- c("eosinophil_count", "eosinophil_percent", "haemoglobin_concentration", "lymphocyte_count", "lymphocyte_percent", "mean_corpuscular_haemoglobin_concentration", "mean_corpuscular_volume", "mean_platelet_volume", "mean_sphered_cell_volume", "neutrophil_count", "neutrophil_percent", "platelet_count", "platelet_distribution_width", "red_blood_cell_count", "red_blood_cell_distribution_width", "white_blood_cell_count")

###################################
coloc_table <- rbind(coloc_table_disease, coloc_table_blood)
coloc_table <- coloc_table %>% mutate(group = case_when(V1 %in% disease_names ~ "disease", V1 %in% blood_count_names ~ "blood_count", TRUE ~ "blood_serum"))

coloc_table$label <- mapping[coloc_table$V1]
p <- ggplot(coloc_table, aes(x = eQTL_prop, y = caQTL_prop, color = group)) + geom_point() + scale_color_manual(values = c("disease" = "#ecc15a", "blood_count" = "#57a7b5", "blood_serum" = "#417e41")) +
  geom_text(aes(label = label), vjust = -0.5, color = "black") + geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + theme_bw() +
  theme(panel.grid = element_blank()) + labs(x = "GWAS loci proportion - eQTL %", y = "GWAS loci proportion - caQTL %") + theme(legend.position = "none")
ggsave("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/figure_H4/coloc_scatter.png", p, width = 7, height = 6.5, dpi = 300)

coloc_table$label <- ifelse(coloc_table$V1 %in% c("parkinson_GCST009325", "covid_GCST011071", "IBD_EAS_EUR", "lymphocyte_percent", "c_reactive_protein", "neutrophil_percent"), mapping[coloc_table$V1], "")

p <- ggplot(coloc_table, aes(x = eQTL_prop, y = caQTL_prop, color = group)) + geom_point() + scale_color_manual(values = c("disease" = "#ecc15a", "blood_count" = "#57a7b5", "blood_serum" = "#417e41"), 
  labels = c("disease" = "disease traits", "blood_count" = "blood_count traits", "blood_serum" = "blood_serum traits")) + geom_text(aes(label = label), vjust = -0.5, color = "black", size = 5) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
  plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + xlim(0, 60) + ylim(0, 72) +
  labs(title = "coloc", x = "GWAS loci proportion - eQTL %", y = "GWAS loci proportion - caQTL %", color = "Group")
ggsave("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/figure_H4_newGWAS/coloc_scatter_label.png", p, width = 5.5, height = 4, dpi = 300)

p <- ggplot(coloc_table, aes(x = eQTL_prop, y = caQTL_prop, color = group)) + geom_point() + scale_color_manual(values = c("disease" = "#ecc15a", "blood_count" = "#57a7b5", "blood_serum" = "#417e41"), 
  labels = c("disease" = "disease traits", "blood_count" = "blood_count traits", "blood_serum" = "blood_serum traits")) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
  plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + xlim(0, 60) + ylim(0, 72) +
  labs(title = "coloc", x = "GWAS loci proportion - eQTL %", y = "GWAS loci proportion - caQTL %", color = "Group")
ggsave("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/figure_H4_newGWAS/coloc_scatter_label_clean.png", p, width = 5.5, height = 4, dpi = 300)

###################################
smr_table <- rbind(smr_table_disease, smr_table_blood)
smr_table <- smr_table %>% mutate(group = case_when(V1 %in% disease_names ~ "disease", V1 %in% blood_count_names ~ "blood_count", TRUE ~ "blood_serum"))

smr_table$label <- mapping[smr_table$V1]
p <- ggplot(smr_table, aes(x = eQTL_prop, y = caQTL_prop, color = group)) + geom_point() + scale_color_manual(values = c("disease" = "#ecc15a", "blood_count" = "#57a7b5", "blood_serum" = "#417e41")) +
  geom_text(aes(label = label), vjust = -0.5, color = "black") + geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + theme_bw() +
  theme(panel.grid = element_blank()) + labs(x = "GWAS loci proportion - eQTL %", y = "GWAS loci proportion - caQTL %") + theme(legend.position = "none")
ggsave("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/figure_H4/smr_scatter.png", p, width = 7, height = 6.5, dpi = 300)

smr_table$label <- ifelse(smr_table$V1 %in% c("parkinson_GCST009325", "covid_GCST011071", "IBD_EAS_EUR", "lymphocyte_percent", "c_reactive_protein", "neutrophil_percent"), mapping[smr_table$V1], "")
p <- ggplot(smr_table, aes(x = eQTL_prop, y = caQTL_prop, color = group)) + geom_point() + scale_color_manual(values = c("disease" = "#ecc15a", "blood_count" = "#57a7b5", "blood_serum" = "#417e41"), 
  labels = c("disease" = "disease traits", "blood_count" = "blood_count traits", "blood_serum" = "blood_serum traits")) + geom_text(aes(label = label), vjust = -0.5, color = "black", size = 4) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
  plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + xlim(0, 60) + ylim(0, 72) +
  labs(title = "SMR", x = "GWAS loci proportion - eQTL %", y = "GWAS loci proportion - caQTL %", color = "Group")
ggsave("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/figure_H4_newGWAS/smr_scatter_label.png", p, width = 5.5, height = 4, dpi = 300)

p <- ggplot(smr_table, aes(x = eQTL_prop, y = caQTL_prop, color = group)) + geom_point() + scale_color_manual(values = c("disease" = "#ecc15a", "blood_count" = "#57a7b5", "blood_serum" = "#417e41"), 
  labels = c("disease" = "disease traits", "blood_count" = "blood_count traits", "blood_serum" = "blood_serum traits")) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
  plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + xlim(0, 60) + ylim(0, 72) +
  labs(title = "SMR", x = "GWAS loci proportion - eQTL %", y = "GWAS loci proportion - caQTL %", color = "Group")
ggsave("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/figure_H4_newGWAS/smr_scatter_label_clean.png", p, width = 5.5, height = 4, dpi = 300)

