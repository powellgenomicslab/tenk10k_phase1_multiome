
library(data.table)
library(ggplot2)
library(tidyverse)

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


disease_names <- c("alzheimer_GCST90027158", "breastca_GCST004988", "colorectalca_GCST90129505", "covid_GCST011071", 
                   "NHL_GCST90011819", "lungca_GCST004748", "lymphoma_GCST90018878", "parkinson_GCST009325", "prostateca_GCST90274713", 
                   "ra_GCST90132223", "sle_GCST003156", "myeloproliferative_GCST90000032", "lymphocytic_leukemia_GCST90011814", 
                   "kiryluk_IgAN", "asthma", "ms", "t1dm", "CD_EAS_EUR", "IBD_EAS_EUR", "UC_EAS_EUR")

blood_count_names <- c("eosinophil_count", "eosinophil_percent", "haemoglobin_concentration", "lymphocyte_count", "lymphocyte_percent", "mean_corpuscular_haemoglobin_concentration", "mean_corpuscular_volume", "mean_platelet_volume", "mean_sphered_cell_volume", "neutrophil_count", "neutrophil_percent", "platelet_count", "platelet_distribution_width", "red_blood_cell_count", "red_blood_cell_distribution_width", "white_blood_cell_count")


coloc_table_disease <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/disease_traits/table_H4_newGWAS/coloc_table.csv")
coloc_table_disease <- coloc_table_disease[!V1 %in% "nephrotic_GCST90258619"]
coloc_table_disease[, None := Num_GWAS_loci - Shared - Specific_eQTL - Specific_caQTL]
coloc_table_disease <- coloc_table_disease[, c("V1", "None", "Specific_caQTL", "Shared", "Specific_eQTL")]

coloc_table_blood <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/blood_traits/table_H4/coloc_table.csv")
coloc_table_blood[, None := Num_GWAS_loci - Shared - Specific_eQTL - Specific_caQTL]
coloc_table_blood <- coloc_table_blood[, c("V1", "None", "Specific_caQTL", "Shared", "Specific_eQTL")]

coloc_table <- rbind(coloc_table_disease, coloc_table_blood)
coloc_table$V1 <- mapping[coloc_table$V1]

df_long <- coloc_table %>% pivot_longer(cols = -V1, names_to = "GWAS_loci_type", values_to = "Count")

df_long <- df_long %>% group_by(V1) %>% mutate(Proportion = Count / sum(Count)) %>% ungroup()

df_long$V1 <- factor(df_long$V1, levels = unique(df_long$V1))
df_long$GWAS_loci_type <- factor(df_long$GWAS_loci_type, levels = c("None", "Specific_caQTL", "Shared", "Specific_eQTL"))

p <- ggplot(df_long, aes(x = V1, y = Proportion, fill = GWAS_loci_type)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 16.5, linetype = "dashed", color = "black") +  # Line between 3rd and 4th bar
  ylab("Proportion") +
  xlab("") +
  scale_fill_manual(
    values = c(
      "None" = "#D3D3D3",
      "Specific_caQTL" = "#5c96a5",
      "Shared" = "#5c5ca5",
      #"Specific_eQTL" = "#f5c767",
      "Specific_eQTL" = "#a55c8e"
    )
  ) +
  # scale_fill_viridis_d() +  # Automatic viridis colors for discrete data
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(t = 5, r = 5, b = 30, l = 30, unit = "pt")
  )
ggsave("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/figure_H4_newGWAS/coloc_stackedbar.png", p, bg = "white")




smr_table_disease <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/disease_traits/table_H4_newGWAS/smr_table.csv")
smr_table_disease <- smr_table_disease[!V1 %in% "nephrotic_GCST90258619"]
smr_table_disease[, None := Num_GWAS_loci - Shared - Specific_eQTL - Specific_caQTL]
smr_table_disease <- smr_table_disease[, c("V1", "None", "Specific_caQTL", "Shared", "Specific_eQTL")]

smr_table_blood <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/blood_traits/table_H4/smr_table.csv")
smr_table_blood[, None := Num_GWAS_loci - Shared - Specific_eQTL - Specific_caQTL]
smr_table_blood <- smr_table_blood[, c("V1", "None", "Specific_caQTL", "Shared", "Specific_eQTL")]

smr_table <- rbind(smr_table_disease, smr_table_blood)
smr_table$V1 <- mapping[smr_table$V1]

df_long <- smr_table %>% pivot_longer(cols = -V1, names_to = "GWAS_loci_type", values_to = "Count")

df_long <- df_long %>% group_by(V1) %>% mutate(Proportion = Count / sum(Count)) %>% ungroup()

df_long$V1 <- factor(df_long$V1, levels = unique(df_long$V1))
df_long$GWAS_loci_type <- factor(df_long$GWAS_loci_type, levels = c("None", "Specific_caQTL", "Shared", "Specific_eQTL"))

p <- ggplot(df_long, aes(x = V1, y = Proportion, fill = GWAS_loci_type)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 16.5, linetype = "dashed", color = "black") +  # Line between 3rd and 4th bar
  ylab("Proportion") +
  xlab("") +
  scale_fill_manual(
    values = c(
      "None" = "#D3D3D3",
      "Specific_caQTL" = "#5c96a5",
      "Shared" = "#5c5ca5",
      #"Specific_eQTL" = "#f5c767",
      "Specific_eQTL" = "#a55c8e"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(t = 5, r = 5, b = 30, l = 30, unit = "pt")
  )
ggsave("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/figure_H4_newGWAS/smr_stackedbar.png", p, bg = "white")