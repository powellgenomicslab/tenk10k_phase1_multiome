
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/lib
smr_tool=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/software/smr-1.3.2-linux-x86_64/smr

mkdir -p /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/SMR_run_rawID

# for cell_type in "${CELL_TYPES[@]}"
# do

cell_type=$1
echo "==============================================="
echo "Processing cell type: $cell_type"
echo "==============================================="

mkdir -p /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/SMR_run_rawID/${cell_type}/

# Loop through chromosomes 1-22
for chr in {1..22}
do
  echo "Processing chromosome ${chr}..."

  ${smr_tool} --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr${chr}_common_variants_qced \
    --beqtl-summary /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/preprocessing_caQTL_rawID/$cell_type/BESD/Chr${chr} \
    --beqtl-summary /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses/preprocessing_eQTL/$cell_type/BESD/Chr${chr} \
    --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/SMR_run_rawID/${cell_type}/Chr${chr}_results \
    >> /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/SMR_run_rawID/${cell_type}/SMR.log 2>&1

done

echo "Completed processing for cell type: $cell_type"
echo ""

# done

echo "All cell types processed."
