## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -o /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/VRE/Logs/1_get_VRE_plink.OU
#$ -e /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/VRE/Logs/1_get_VRE_plink.ER
#$ -N get_VRE_plink
#$ -q short.q
#$ -pe smp 1
#$ -l mem_requested=40G
#$ -l tmp_requested=30G
#$ -r yes
#$ -t 1-22

CHR=${SGE_TASK_ID}
GENODIR=/directflow/SCCGGroupShare/projects/angxue/data/tenk10k_phase1/genotype/from_wgs/filtered
OUTDIR=/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/VRE/PrunedVariants

# 1. Initial filtering
/directflow/SCCGGroupShare/projects/oscdon/software/PLINK2/plink2 \
    --bfile ${GENODIR}/TenK10K_TOB_ATAC_renamed_chr${CHR}_common_variants_qced \
    --snps-only \
    --max-alleles 2 \
    --mac 20 \
    --make-bed \
    --out ${TMPDIR}/temp_filtered

# 2. LD pruning
/directflow/SCCGGroupShare/projects/oscdon/software/PLINK2/plink2 \
    --bfile ${TMPDIR}/temp_filtered \
    --indep-pairwise 500000 0.2 \
    --out ${TMPDIR}/pruning_parameters

# 3. Extract pruned variants
/directflow/SCCGGroupShare/projects/oscdon/software/PLINK2/plink2 \
    --bfile ${TMPDIR}/temp_filtered \
    --extract ${TMPDIR}/pruning_parameters.prune.in \
    --make-bed \
    --out ${OUTDIR}/chr${CHR}_pruned
