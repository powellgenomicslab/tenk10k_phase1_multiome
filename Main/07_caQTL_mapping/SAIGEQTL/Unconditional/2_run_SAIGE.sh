#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N celltype_chrnum
#PBS -l walltime=9:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l ncpus=28
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/OutputLogs/2_run_SAIGE/celltype/chrchrnum.OU
#PBS -e /g/data/fy54/od8037/Brenner/SAIGEQTL/Rare/OutputLogs/2_run_SAIGE/celltype/chrchrnum.ER

eval "$(pixi shell-hook --manifest-path=/g/data/fy54/od8037/software/SAIGEQTL/pixi.toml -s bash)"
export CONDA_OVERRIDE_GLIBC=2.28

export CELLTYPE=celltype
export CHR_NUM=chrnum
echo "Processing ${CELLTYPE} and chromosome ${CHR_NUM}"

SAIGE_DIR=/g/data/fy54/od8037/Brenner/SAIGEQTL
export VRE_FILE=${SAIGE_DIR}/VRE/FinalVRE/VRE_geno
export PHENO_DIR=${SAIGE_DIR}/Phenotypes/FinalTSVs
export GENO_DIR=${SAIGE_DIR}/Rare/Genotypes
export REGION_FILE=${SAIGE_DIR}/Rare/Regions/${CELLTYPE}/chr${CHR_NUM}.csv
export LOGFILE=${SAIGE_DIR}/Rare/Logs/2_run_SAIGE/${CELLTYPE}/chr${CHR_NUM}.log
export OUT_DIR=${SAIGE_DIR}/Rare/Results/${CELLTYPE}/chr${CHR_NUM}

mkdir -p "$(dirname "$LOGFILE")"
mkdir -p "$OUT_DIR"
touch "$LOGFILE"
> "$LOGFILE"

# Define function to process a single peak
run_saige_peak() {
    local PEAK_ID=$1
    local CHR=$2
    local START=$3
    local END=$4
    
    local TASK_TMP="${TMPDIR}/${PEAK_ID}"
    mkdir -p "$TASK_TMP"

    # Step 1: VRE
    Rscript /g/data/fy54/od8037/software/SAIGEQTL/extdata/step1_fitNULLGLMM_qtl.R \
        --phenoFile=${PHENO_DIR}/${CELLTYPE}/chr${CHR_NUM}.tsv \
        --plinkFile=${VRE_FILE} \
        --outputPrefix=${TASK_TMP}/${PEAK_ID} \
        --phenoCol=${PEAK_ID} \
        --covarColList=nCount_peaks,sex,geno_pc1,geno_pc2,geno_pc3,geno_pc4,geno_pc5,geno_pc6,age,PC1,PC2,PC3,PC4,PC5 \
        --sampleCovarColList=sex,geno_pc1,geno_pc2,geno_pc3,geno_pc4,geno_pc5,geno_pc6,age \
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
        --nThreads=1 > /dev/null 2>&1 || { echo "[$(date '+%F %T')] FAILED Step 1: $PEAK_ID" >> "$LOGFILE"; rm -rf "$TASK_TMP"; return 1; }

    # Step 1.5: Make Group File
    echo -e "${CHR}\t${START}\t${END}" > "${TASK_TMP}/${PEAK_ID}_region.txt"
    
    Rscript /g/data/fy54/od8037/software/SAIGEQTL/extdata/makeGroupFile.R \
        --vcfFile=${GENO_DIR}/TenK10K_TOB_ATAC_renamed_chr${CHR}_combined_lowfreq_rare_maf005.vcf.gz \
        --vcfFileIndex=${GENO_DIR}/TenK10K_TOB_ATAC_renamed_chr${CHR}_combined_lowfreq_rare_maf005.vcf.gz.csi \
        --vcfField=GT \
        --regionFile=${TASK_TMP}/${PEAK_ID}_region.txt \
        --outputPrefix=${TASK_TMP}/${PEAK_ID}.grp > /dev/null 2>&1 || { echo "[$(date '+%F %T')] FAILED Step 1.5: $PEAK_ID" >> "$LOGFILE"; rm -rf "$TASK_TMP"; return 1; }

    # Step 2: Set-Based Tests
    Rscript /g/data/fy54/od8037/software/SAIGEQTL/extdata/step2_tests_qtl.R \
        --vcfFile=${GENO_DIR}/TenK10K_TOB_ATAC_renamed_chr${CHR_NUM}_combined_lowfreq_rare_maf005.vcf.gz \
        --vcfFileIndex=${GENO_DIR}/TenK10K_TOB_ATAC_renamed_chr${CHR_NUM}_combined_lowfreq_rare_maf005.vcf.gz.csi \
        --SAIGEOutputFile=${OUT_DIR}/${PEAK_ID} \
        --chrom=${CHR_NUM} \
        --GMMATmodelFile=${TASK_TMP}/${PEAK_ID}.rda \
        --varianceRatioFile=${TASK_TMP}/${PEAK_ID}.varianceRatio.txt \
        --groupFile=${TASK_TMP}/${PEAK_ID}.grp_${CHR_NUM} \
        --annotation_in_groupTest=null \
        --vcfField=GT \
        --maxMAF_in_groupTest=0.05 \
        --minMAF_in_groupTest_Exclude=0 \
        --MACCutoff_to_CollapseUltraRare=10 \
        --is_single_in_groupTest=TRUE \
        --is_equal_weight_in_groupTest=TRUE \
        --LOCO=FALSE \
        --markers_per_chunk=10000 \
        --SPAcutoff=2 > /dev/null 2>&1 || { echo "[$(date '+%F %T')] FAILED Step 2: $PEAK_ID" >> "$LOGFILE"; rm -rf "$TASK_TMP"; return 1; }

    echo "[$(date '+%F %T')] SUCCESS: $PEAK_ID" >> "$LOGFILE"

    # Remove unnecessary files
    rm -rf "$TASK_TMP"
    rm -f "${OUT_DIR}/${PEAK_ID}.index"
}
export -f run_saige_peak

module load parallel
cat "${REGION_FILE}" | parallel --colsep '\t' -j 28 run_saige_peak {1} {2} {3} {4}

echo "Finished chromosome ${CHR_NUM}"
