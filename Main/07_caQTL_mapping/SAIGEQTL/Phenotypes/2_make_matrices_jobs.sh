# Define file paths
TEMPLATE="/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/2_make_matrices.qsub.sh"
OUT_DIR="/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Scripts/2_make_matrices"
CELLTYPES="/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/UNFINISHED.txt"

# Create jobs
for i in {1..7}; do

    # Extract the i-th cell type from the file
    CELLTYPE=$(sed -n "${i}p" "$CELLTYPES")
    echo "${CELLTYPE}"

    mkdir -p "/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Logs/2_make_matrices/${CELLTYPE}"

    OUT_FILE="$OUT_DIR/${CELLTYPE}.qsub.sh"
    sed -e "s/celltype/${CELLTYPE}/g" "$TEMPLATE" > "$OUT_FILE"
done

# Submit jobs
for i in {1..7}; do

    # Extract the i-th cell type from the file
    CELLTYPE=$(sed -n "${i}p" "$CELLTYPES")

    QSUB_SCRIPT="$OUT_DIR/${CELLTYPE}.qsub.sh"
    if [[ -f "$QSUB_SCRIPT" ]]; then
        echo "Submitting job for ${CELLTYPE}"
        qsub "$QSUB_SCRIPT"
    else
        echo "Warning: qsub script not found for ${CELLTYPE}: $QSUB_SCRIPT"
    fi
done
