
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/lib
smr_tool=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/software/smr-1.3.2-linux-x86_64/smr

disease=$1
cell_type=$2

# mkdir -p /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_run_rawID

CELL_TYPES_DIR=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/preprocessing_caQTL_rawID/
CELL_TYPES=($(ls -1 "$CELL_TYPES_DIR"))

# for cell_type in "${CELL_TYPES[@]}"
# do

echo "==============================================="
echo "Processing cell type: $cell_type"
echo "==============================================="

mv /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_run_rawID/${cell_type}/disease_traits/${disease} /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_run_rawID/${cell_type}/disease_traits/${disease}_wrongGWAS
# mkdir -p /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_run_rawID/${cell_type}/disease_traits/
mkdir -p /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_run_rawID/${cell_type}/disease_traits/${disease}
  
# Loop through chromosomes 1-22
for chr in {1..22}
do
  echo "Processing chromosome ${chr}..."

  ${smr_tool} --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr${chr}_common_variants_qced \
    --gwas-summary /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/GWAS/Additional_Set/ma_no_INDEL_RARE_REPEAT/${disease}/chr${chr}.ma \
    --beqtl-summary /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/preprocessing_caQTL_rawID/$cell_type/BESD/Chr${chr} \
    --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_run_rawID/${cell_type}/disease_traits/${disease}/Chr${chr}_results \
    >> /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/GWAS_run_rawID/${cell_type}/disease_traits/${disease}/SMR.log 2>&1
  
done

# done
  
echo "Completed processing for trait: $disease cell type: $cell_type"
echo ""
# done

echo "All cell types processed."

# --gwas-summary /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/GWAS/Flagship_disease_trait/ma_no_INDEL_RARE_REPEAT/${disease}/chr${chr}.ma \