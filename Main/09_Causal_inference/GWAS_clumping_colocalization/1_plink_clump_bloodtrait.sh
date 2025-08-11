
plink_tool=/directflow/SCCGGroupShare/projects/jayfan/Software/PLINK/plink

disease=$1

mkdir -p /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/plink_output/blood_traits_newGWAS/${disease}/

for chr in {1..22}
do
  
  sed -i '1s/\bp\b/P/' /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/GWAS/Flagship_blood_trait/ma_no_INDEL_RARE_REPEAT_newGWAS_cp/${disease}_chr${chr}.ma

  ${plink_tool} --bfile /directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered/TenK10K_TOB_ATAC_renamed_chr${chr}_common_variants_qced \
    --clump /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/GWAS/Flagship_blood_trait/ma_no_INDEL_RARE_REPEAT_newGWAS_cp/${disease}_chr${chr}.ma \
    --clump-p1 5e-8 --clump-p2 1e-5 --clump-r2 0.01 --clump-kb 1000 \
    --out /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/plink_output/blood_traits_newGWAS/${disease}/clumped_chr${chr} \
    >> /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/coloc_SMR_compare/PLINK_clumping/plink_output/blood_traits_newGWAS/${disease}/clump.log 2>&1

done