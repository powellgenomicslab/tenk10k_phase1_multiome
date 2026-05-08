#!/bin/bash
# Define file paths
TEMPLATE="/g/data/fy54/od8037/TenK10K/ColocNew/DiseaseTraits/2_coloc.qsub.sh"
BASE_OUT_DIR="/g/data/fy54/od8037/TenK10K/ColocNew/DiseaseTraits/Scripts/2_coloc"
CELLTYPES="/g/data/ei56/od8037/Final_caQTL/celltype_names.txt"
CONDITIONS="/g/data/fy54/od8037/TenK10K/ColocNew/DiseaseTraits/conditions.txt"

# Loop through cell types
for i in {1..26}; do
    # Extract the i-th cell type from the file
    CELLTYPE=$(sed -n "${i}p" "$CELLTYPES")
    echo "${CELLTYPE}"

    # Loop through condition numbers
    for j in {1..15}; do
        # Extract the j-th condition from the file
        CONDITION=$(sed -n "${j}p" "$CONDITIONS")
        mkdir -p "/g/data/fy54/od8037/TenK10K/ColocNew/DiseaseTraits/OutputLogs/2_coloc/$CONDITION"
        mkdir -p "/g/data/fy54/od8037/TenK10K/ColocNew/DiseaseTraits/Logs/2_coloc/$CONDITION"
        mkdir -p "$BASE_OUT_DIR/$CONDITION"

        # Define output file name
        OUT_FILE="$BASE_OUT_DIR/${CONDITION}/${CELLTYPE}.qsub.sh"

        # Replace placeholders and generate qsub script
        sed -e "s/ct_name/${CELLTYPE}/g" \
            -e "s/condition/${CONDITION}/g" \
            "$TEMPLATE" > "$OUT_FILE"
    done
done

# Loop through all conditions 
for j in {1..15}; do
    CONDITION=$(sed -n "${j}p" "$CONDITIONS")

    # Loop through all cell types (1 to 26)
    for i in {1..26}; do
        CELLTYPE=$(sed -n "${i}p" "$CELLTYPES")

        # Construct expected qsub script path
        QSUB_SCRIPT="${BASE_OUT_DIR}/${CONDITION}/${CELLTYPE}.qsub.sh"

        # Submit if the script exists
        if [[ -f "$QSUB_SCRIPT" ]]; then
            echo "Submitting: $QSUB_SCRIPT"
            qsub "$QSUB_SCRIPT"
        else
            echo "Warning: Not found - $QSUB_SCRIPT"
        fi
    done
done

