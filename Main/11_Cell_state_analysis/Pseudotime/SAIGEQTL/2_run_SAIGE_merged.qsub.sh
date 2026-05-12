#!/bin/bash
#PBS -P fy54
#PBS -q normal
#PBS -N merge_chrNumber
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=10GB
#PBS -l jobfs=2GB
#PBS -l ncpus=2
#PBS -r y
#PBS -l wd
#PBS -o /g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Logs/3_run_SAIGE_merged/chr_chrNumber.OU
#PBS -e /g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Logs/3_run_SAIGE_merged/chr_chrNumber.ER

module load singularity
module load parallel

export CHR_NUM=chrNumber
export ITER=JOB_INDEX
echo "Processing merged objects for chromosome ${CHR_NUM}"

SAIGE_DIR=/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL
export VRE_FILE=${SAIGE_DIR}/Common/VRE/FinalVRE/VRE_geno
export PHENO_DIR=/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Phenotypes/FinalTSVs_merged
export PLINK_GENO_DIR=/g/data/ei56/ax3061/proj/tenk10k/caQTL/data/genotype
export REGION_FILE=${SAIGE_DIR}/Dynamic/Regions/chr${CHR_NUM}.csv
export LOGFILE=/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Logs/3_run_SAIGE_merged/chr${CHR_NUM}.log
export OUT_DIR=/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/chr${CHR_NUM}
export TMP_DIR=/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/chr${CHR_NUM}_tmp

mkdir -p "$(dirname "$LOGFILE")"
mkdir -p "$OUT_DIR"
mkdir -p "$TMP_DIR"
touch "$LOGFILE"
> "$LOGFILE"

# Define function to process a single peak
run_saige_peak() {
    local PEAK_ID=$1
    local CHR=$2
    local START=$3
    local END=$4

    local TASK_TMP="${TMP_DIR}/${PEAK_ID}"
    mkdir -p "$TASK_TMP"
    mkdir -p "${OUT_DIR}/${PEAK_ID}"

    echo "[$(date '+%F %T')] $PEAK_ID" >> "$LOGFILE"

    # Step 1: VRE
    singularity exec --bind /g/data:/g/data /g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/software/saigeqtldynamic0.2.5.1_latest.sif step1_fitNULLGLMM_qtl.R  \
        --phenoFile=${PHENO_DIR}/${CELLTYPE}/chr${CHR_NUM}.tsv \
        --plinkFile=${VRE_FILE} \
        --outputPrefix=${TASK_TMP}/${PEAK_ID} \
        --phenoCol=${PEAK_ID} \
        --covarColList=repeat_num,nCount_peaks,sex,geno_pc1,geno_pc2,geno_pc3,geno_pc4,geno_pc5,geno_pc6,age,PC1,PC2,PC3,PC4,PC5,pseudotime \
        --sampleCovarColList=sex,geno_pc1,geno_pc2,geno_pc3,geno_pc4,geno_pc5,geno_pc6,age \
        --dynamicCovarColList=pseudotime \
        --sampleIDColinphenoFile=donor_id \
        --traitType=count \
        --skipVarianceRatioEstimation=FALSE \
        --isRemoveZerosinPheno=FALSE \
        --useSparseGRMtoFitNULL=FALSE \
        --useGRMtoFitNULL=FALSE \
        --isCovariateOffset=FALSE \
        --isCovariateTransform=TRUE \
        --skipModelFitting=FALSE \
        --tol=1e-6 \
        --IsOverwriteVarianceRatioFile=TRUE \
        --isStoreSigma=TRUE \
        --tauInit=1,0.1,0 \
        --nThreads=2 > /dev/null 2>&1

    # Step 1.5: Make Region File
    echo -e "${CHR}\t${START}\t${END}" > "${TASK_TMP}/${PEAK_ID}_region.txt"

    # Step 2: caQTL Mapping
    singularity exec --bind /g/data:/g/data /g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/software/saigeqtldynamic0.2.5.1_latest.sif step2_tests_qtl.R \
        --bedFile=${PLINK_GENO_DIR}/TenK10K_TOB_ATAC_renamed_chr${CHR_NUM}_common_variants_qced.bed \
        --bimFile=${PLINK_GENO_DIR}/TenK10K_TOB_ATAC_renamed_chr${CHR_NUM}_common_variants_qced.bim \
        --famFile=${PLINK_GENO_DIR}/TenK10K_TOB_ATAC_renamed_chr${CHR_NUM}_common_variants_qced.fam \
        --SAIGEOutputFile=${OUT_DIR}/${PEAK_ID} \
        --chrom=${CHR_NUM} \
        --GMMATmodelFile=${TASK_TMP}/${PEAK_ID}.rda \
        --varianceRatioFile=${TASK_TMP}/${PEAK_ID}.varianceRatio.txt \
        --rangestoIncludeFile=${TASK_TMP}/${PEAK_ID}_region.txt \
        --vcfField=GT \
        --minMAF=0.01 \
        --minMAC=5 \
        --LOCO=FALSE \
        --markers_per_chunk=10000 \
        --SPAcutoff=2 \
        --pval_cutoff_for_gxe=1 \
        --is_permute_e=FALSE \
        --is_permute_ginge=FALSE > /dev/null 2>&1

    # Step 3: Compute p-values
    singularity exec --bind /g/data:/g/data /g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/software/saigeqtldynamic0.2.5.1_latest.sif step3_gene_pvalue_qtl.R \
        --assocFile=${OUT_DIR}/${PEAK_ID} \
        --geneName=${PEAK_ID} \
        --genePval_outputFile=${OUT_DIR}/${PEAK_ID}_ACAT > /dev/null 2>&1

    # Remove unnecessary files
    rm -rf "$TASK_TMP"
    rm -f "${OUT_DIR}/${PEAK_ID}.index"
}

export -f run_saige_peak

TOTAL_LINES=$(wc -l < "${REGION_FILE}")
LINES_PER_JOB=$((TOTAL_LINES / 40))

START=$(( (ITER - 1) * LINES_PER_JOB + 1 ))

if [ $ITER -eq 40 ]; then
    END=$TOTAL_LINES
else
    END=$(( ITER * LINES_PER_JOB ))
fi

sed -n "${START},${END}p" "${REGION_FILE}" | while IFS=$'\t' read -r col1 col2 col3 col4; do
    run_saige_peak "$col1" "$col2" "$col3" "$col4"
done

echo "Finished chromosome ${CHR_NUM} job id ${ITER}"
