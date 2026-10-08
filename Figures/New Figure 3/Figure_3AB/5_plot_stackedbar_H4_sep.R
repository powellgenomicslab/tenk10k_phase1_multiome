
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
  "rheumatoid arthritis",
  "sle",
  "myeloproliferative",
  "lymphocytic_leukemia",
  "nephrotic",
  "IgA nephropathy",
  "asthma",
  "multiple sclerosis",
  "t1dm",
  "crohns",
  "ibd",
  "ulcerative colitis",
  "alanine amino", 
  "albumin", 
  "alkaline pho", 
  "apolipoprotein a", 
  "apolipoprotein b", 
  "aspartate amino", 
  "c reactive protein", 
  "calcium", 
  "cholesterol", 
  "creatinine", 
  "cystatin c", 
  "eosinophil count", 
  "eosinophil percent", 
  "gamma glutamyl", 
  "glucose", 
  "glycated haem", 
  "haematocrit", 
  "haemoglobin conc", 
  "hdl cholesterol", 
  "igf 1", 
  "ldl cholesterol direct", 
  "lymphocyte count", 
  "lymphocyte percent", 
  "mean corpus haem conc", 
  "mean corpus haem", 
  "mean corpus volume", 
  "mean platelet volume", 
  "mean sphered cell volume", 
  "neutrophil count", 
  "neutrophil percent", 
  "phosphate", 
  "platelet count", 
  "platelet crit", 
  "platelet dist width", 
  "red blood cell count", 
  "red blood cell dist width", 
  "shbg", 
  "total bilirubin", 
  "total protein", 
  "triglycerides", 
  "urate", 
  "urea", 
  "vitamin d", 
  "white blood cell count"
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
  "cd",
  "ibd",
  "uc",
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
                   "kiryluk_IgAN", "asthma", "ms", "t1dm", "cd", "ibd", "uc")

blood_count_names <- c("eosinophil_count", "eosinophil_percent", "haemoglobin_concentration", "lymphocyte_count", "lymphocyte_percent", "mean_corpuscular_haemoglobin_concentration", "mean_corpuscular_volume", "mean_platelet_volume", "mean_sphered_cell_volume", "neutrophil_count", "neutrophil_percent", "platelet_count", "platelet_distribution_width", "red_blood_cell_count", "red_blood_cell_distribution_width", "white_blood_cell_count")


coloc_table_disease <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final_1Mb/disease_traits/table_H4_newGWAS/coloc_table.csv")

coloc_table_disease[, None := Num_GWAS_loci - Shared - Specific_eQTL - Specific_caQTL]
coloc_table_disease <- coloc_table_disease[, c("V1", "None", "Specific_caQTL", "Shared", "Specific_eQTL")]

coloc_table_blood <- fread("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final_1Mb/blood_traits/table_H4/coloc_table.csv")
coloc_table_blood[, None := Num_GWAS_loci - Shared - Specific_eQTL - Specific_caQTL]
coloc_table_blood <- coloc_table_blood[, c("V1", "None", "Specific_caQTL", "Shared", "Specific_eQTL")]

coloc_table <- rbind(coloc_table_disease, coloc_table_blood)

# Add a column to annotate each row is from which group
coloc_table_all <- coloc_table %>% mutate(group = case_when(V1 %in% disease_names ~ "disease", V1 %in% blood_count_names ~ "blood_count", TRUE ~ "blood_serum"))
coloc_table <- coloc_table_all[group == "disease"]
coloc_table[, group := NULL]

coloc_table$V1 <- mapping[coloc_table$V1]

df_long <- coloc_table %>% pivot_longer(cols = -V1, names_to = "GWAS_loci_type", values_to = "Count")

df_long <- df_long %>% group_by(V1) %>% mutate(Proportion = Count / sum(Count)) %>% ungroup()

df_long$V1 <- factor(df_long$V1, levels = unique(df_long$V1))
df_long$GWAS_loci_type <- factor(df_long$GWAS_loci_type, levels = c("None", "Specific_caQTL", "Shared", "Specific_eQTL"))

p <- ggplot(df_long, aes(x = V1, y = Proportion, fill = GWAS_loci_type)) +
  geom_bar(stat = "identity") +
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
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(t = 5, r = 5, b = 30, l = 30, unit = "pt")
  )
ggsave("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final_1Mb/figure_H4_newGWAS/coloc_stackedbar_disease.png", p, width = 6, height = 5, dpi = 300, bg = "white")
