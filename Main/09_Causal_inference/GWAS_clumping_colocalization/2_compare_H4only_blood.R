
library(data.table)

args = commandArgs(trailingOnly=TRUE)
trait <- args[1]
# cell_type_input <- args[2]

# Define condition names
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

plink_result_dir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/plink_output/blood_traits_newGWAS/"
eQTL_result_basedir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_OldGenotypes_not_used_for_final_manuscript/preprocessing_eQTL/"
caQTL_result_basedir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/preprocessing_caQTL_rawID/"

eQTL_coloc_result_basedir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Coloc_finalresults/blood_eQTL_GWAS/"
caQTL_coloc_result_basedir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/Coloc_finalresults/blood_caQTL_GWAS/"
eQTL_SMR_result_basedir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_eQTL_newUKBBGWAS/"
caQTL_SMR_result_basedir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_run_rawID_newUKBBGWAS/"

trait_list <- pheno_names
cell_type_list <- sort(list.files(caQTL_SMR_result_basedir))
cell_type_list_dual <- paste0(rep(cell_type_list, each = 2), c("_eQTL", "_caQTL"))
cell_type_list_dual <- c("Union_eQTL", "Union_caQTL", "Num_GWAS_loci", "Shared", "Specific_eQTL", "Specific_caQTL", cell_type_list_dual)

# Create empty matrix with trait_list as rows and cell_type_list as columns
coloc_table <- matrix(NA, nrow = length(trait_list), ncol = length(cell_type_list_dual))
rownames(coloc_table) <- pheno_names
colnames(coloc_table) <- cell_type_list_dual
coloc_table <- as.data.frame(coloc_table)
smr_table <- matrix(NA, nrow = length(trait_list), ncol = length(cell_type_list_dual))
rownames(smr_table) <- pheno_names
colnames(smr_table) <- cell_type_list_dual
smr_table <- as.data.frame(smr_table)

#for (trait in trait_list) {

sig_SNP_union <- list()
eQTL_coloc_union <- list()
caQTL_coloc_union <- list()
eQTL_smr_union <- list()
caQTL_smr_union <- list()

# cell_type_list <- c("CD4_Naive", "CD14_Mono")

trait_fullname <- paste0("white_british_", trait, "_snp_gwas_results_hg38")

for (cell_type in cell_type_list) {
  
  eQTL_result_dir <- paste0(eQTL_result_basedir, cell_type, "/matrixQTL/")
  caQTL_result_dir <- paste0(caQTL_result_basedir, cell_type, "/matrixQTL/")
  
  eQTL_SMR_result_dir <- paste0(eQTL_SMR_result_basedir, cell_type, "/blood_traits/", trait_fullname, "/")
  caQTL_SMR_result_dir <- paste0(caQTL_SMR_result_basedir, cell_type, "/blood_traits/", trait_fullname, "/")
  
  sig_SNP_union_ct <- list()
  eQTL_coloc_union_ct <- list()
  caQTL_coloc_union_ct <- list()
  eQTL_smr_union_ct <- list()
  caQTL_smr_union_ct <- list()
  
  for (chr in seq(1, 22)) {
  
    print(paste0(trait, cell_type, chr))
    # Read the significant SNP list from PLINK clumping outputs
    if (file.exists(file.path(plink_result_dir, trait_fullname, paste0("clumped_chr", chr, ".clumped"))) == FALSE) {
      next
    }
    plink_df <- fread(file.path(plink_result_dir, trait_fullname, paste0("clumped_chr", chr, ".clumped")))
    sig_SNP_list <- unique(plink_df$SNP)
    sig_SNP_union_ct <- c(sig_SNP_union_ct, sig_SNP_list)
    
    # Read eQTL-GWAS coloc results
    if (file.exists(file.path(eQTL_coloc_result_basedir, trait, cell_type, paste0("chr", chr, ".csv"))) == FALSE) {
      next
    }
    eQTL_coloc_df <- fread(file.path(eQTL_coloc_result_basedir, trait, cell_type, paste0("chr", chr, ".csv")))
    eQTL_coloc_df[, PP.H3H4 := PP.H3.abf + PP.H4.abf]
    # eQTL_coloc_df <- eQTL_coloc_df[PP.H3H4 > 0.8]
    eQTL_coloc_df <- eQTL_coloc_df[PP.H4.abf > 0.8]
    
    # Read caQTL-GWAS coloc results
    if (file.exists(file.path(caQTL_coloc_result_basedir, trait, cell_type, paste0("chr", chr, ".csv"))) == FALSE) {
      next
    }
    caQTL_coloc_df <- fread(file.path(caQTL_coloc_result_basedir, trait, cell_type, paste0("chr", chr, ".csv")))
    caQTL_coloc_df[, PP.H3H4 := PP.H3.abf + PP.H4.abf]
    # caQTL_coloc_df <- caQTL_coloc_df[PP.H3H4 > 0.8]
    caQTL_coloc_df <- caQTL_coloc_df[PP.H4.abf > 0.8]
    
    # Read eQTL-GWAS SMR results
    if (file.exists(file.path(eQTL_SMR_result_dir, paste0("Chr", chr, "_results.smr"))) == FALSE) {
      next
    }
    eQTL_smr_df <- fread(file.path(eQTL_SMR_result_dir, paste0("Chr", chr, "_results.smr")))
    eQTL_smr_df <- eQTL_smr_df[eQTL_smr_df$p_SMR < (0.005 / nrow(eQTL_smr_df))]
    eQTL_smr_df <- eQTL_smr_df[eQTL_smr_df$p_HEIDI > 5e-8]
    
    # Read caQTL-GWAS SMR results
    if (file.exists(file.path(caQTL_SMR_result_dir, paste0("Chr", chr, "_results.smr"))) == FALSE) {
      next
    }
    caQTL_smr_df <- fread(file.path(caQTL_SMR_result_dir, paste0("Chr", chr, "_results.smr")))
    caQTL_smr_df <- caQTL_smr_df[caQTL_smr_df$p_SMR < (0.005 / nrow(caQTL_smr_df))]
    caQTL_smr_df <- caQTL_smr_df[caQTL_smr_df$p_HEIDI > 5e-8]
    
    # Find genes holding significant eQTLs with PLINK sig SNPs
    eQTL_df <- fread(file.path(eQTL_result_dir, paste0("Chr", chr, "_MatrixcaQTL.tsv")))
    eQTL_df <- eQTL_df[eQTL_df$'p-value' < 5e-8]
    eQTL_df <- eQTL_df[SNP %in% sig_SNP_list]
    sig_genes_eQTL <- unique(eQTL_df$gene)
    
    # Find peaks holding significant caQTLs with PLINK sig SNPs
    caQTL_df <- fread(file.path(caQTL_result_dir, paste0("Chr", chr, "_MatrixcaQTL.tsv")))
    caQTL_df <- caQTL_df[caQTL_df$'p-value' < 5e-8]
    caQTL_df <- caQTL_df[SNP %in% sig_SNP_list]
    sig_peaks_caQTL <- unique(caQTL_df$peak)
    
    # Get the list of sig SNPs colocalized with eQTL_coloc
    eQTL_coloc_sig_genes <- unique(eQTL_coloc_df[gene %in% sig_genes_eQTL]$gene)
    eQTL_coloc_sig_SNPs <- unique(eQTL_df[gene %in% eQTL_coloc_sig_genes]$SNP)
    eQTL_coloc_union_ct <- c(eQTL_coloc_union_ct, eQTL_coloc_sig_SNPs)
    
    # Get the list of sig SNPs colocalized with caQTL_coloc
    caQTL_coloc_sig_peaks <- unique(caQTL_coloc_df[peak %in% sig_peaks_caQTL]$peak)
    caQTL_coloc_sig_SNPs <- unique(caQTL_df[peak %in% caQTL_coloc_sig_peaks]$SNP)
    caQTL_coloc_union_ct <- c(caQTL_coloc_union_ct, caQTL_coloc_sig_SNPs)
    
    # Get the list of sig SNPs colocalized with eQTL_SMR
    eQTL_smr_sig_genes <- unique(eQTL_smr_df[probeID %in% sig_genes_eQTL]$probeID)
    eQTL_smr_sig_SNPs <- unique(eQTL_df[gene %in% eQTL_smr_sig_genes]$SNP)
    eQTL_smr_union_ct <- c(eQTL_smr_union_ct, eQTL_smr_sig_SNPs)
    
    # Get the list of sig SNPs colocalized with caQTL_SMR
    caQTL_smr_sig_peaks <- unique(caQTL_smr_df[probeID %in% sig_peaks_caQTL]$probeID)
    caQTL_smr_sig_SNPs <- unique(caQTL_df[peak %in% caQTL_smr_sig_peaks]$SNP)
    caQTL_smr_union_ct <- c(caQTL_smr_union_ct, caQTL_smr_sig_SNPs)
  }
  
  sig_SNP_union <- c(sig_SNP_union, sig_SNP_union_ct)
  eQTL_coloc_union <- c(eQTL_coloc_union, eQTL_coloc_union_ct)
  caQTL_coloc_union <- c(caQTL_coloc_union, caQTL_coloc_union_ct)
  eQTL_smr_union <- c(eQTL_smr_union, eQTL_smr_union_ct)
  caQTL_smr_union <- c(caQTL_smr_union, caQTL_smr_union_ct)
  
  coloc_table[trait, paste0(cell_type, "_eQTL")] <- length(eQTL_coloc_union_ct)
  coloc_table[trait, paste0(cell_type, "_caQTL")] <- length(caQTL_coloc_union_ct)
  smr_table[trait, paste0(cell_type, "_eQTL")] <- length(eQTL_smr_union_ct)
  smr_table[trait, paste0(cell_type, "_caQTL")] <- length(caQTL_smr_union_ct)
}

coloc_table[trait, "Union_eQTL"] <- length(unique(eQTL_coloc_union))
coloc_table[trait, "Union_caQTL"] <- length(unique(caQTL_coloc_union))
coloc_table[trait, "Num_GWAS_loci"] <- length(unique(sig_SNP_union))
smr_table[trait, "Union_eQTL"] <- length(unique(eQTL_smr_union))
smr_table[trait, "Union_caQTL"] <- length(unique(caQTL_smr_union))
smr_table[trait, "Num_GWAS_loci"] <- length(unique(sig_SNP_union))

# add qtl-specific and shared
shared_coloc_union <- intersect(unique(eQTL_coloc_union), unique(caQTL_coloc_union))
eQTL_coloc_union_spec <- setdiff(unique(eQTL_coloc_union), unique(caQTL_coloc_union))
caQTL_coloc_union_spec <- setdiff(unique(caQTL_coloc_union), unique(eQTL_coloc_union))
shared_smr_union <- intersect(unique(eQTL_smr_union), unique(caQTL_smr_union))
eQTL_smr_union_spec <- setdiff(unique(eQTL_smr_union), unique(caQTL_smr_union))
caQTL_smr_union_spec <- setdiff(unique(caQTL_smr_union), unique(eQTL_smr_union))
# union(eQTL_coloc_union, caQTL_coloc_union)

coloc_table[trait, "Shared"] <- length(unique(shared_coloc_union))
coloc_table[trait, "Specific_eQTL"] <- length(unique(eQTL_coloc_union_spec))
coloc_table[trait, "Specific_caQTL"] <- length(unique(caQTL_coloc_union_spec))
smr_table[trait, "Shared"] <- length(unique(shared_smr_union))
smr_table[trait, "Specific_eQTL"] <- length(unique(eQTL_smr_union_spec))
smr_table[trait, "Specific_caQTL"] <- length(unique(caQTL_smr_union_spec))
  
# }

write.csv(coloc_table, paste0("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/blood_traits/table_H4/coloc_", trait, ".csv"), row.names = TRUE)
write.csv(smr_table, paste0("/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/sta_output_final/blood_traits/table_H4/smr_", trait, ".csv"), row.names = TRUE)

