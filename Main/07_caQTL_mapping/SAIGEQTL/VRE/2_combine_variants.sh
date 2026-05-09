VRE_DIR=/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/VRE/FinalVRE

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

# Merge genotype files from each chromosome
/directflow/SCCGGroupShare/projects/oscdon/software/PLINK/plink \
    --bfile /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/VRE/PrunedVariants/chr1_pruned \
    --merge-list $VRE_DIR/merge_list.txt \
    --make-bed \
    --keep-allele-order \
    --out /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/VRE/PrunedVariants/all_pruned

# Randomly select 2000 variants from the merged file
shuf --random-source=<(get_seeded_random 123) -n 2000 /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/VRE/PrunedVariants/all_pruned.bim \
    | awk '{print $2}' \
    > $VRE_DIR/VRE_list.txt

# Export 2000 variants into final file
/directflow/SCCGGroupShare/projects/oscdon/software/PLINK2/plink2 \
    --bfile /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/VRE/PrunedVariants/all_pruned \
    --extract $VRE_DIR/VRE_list.txt \
    --make-bed \
    --out $VRE_DIR/VRE_geno
