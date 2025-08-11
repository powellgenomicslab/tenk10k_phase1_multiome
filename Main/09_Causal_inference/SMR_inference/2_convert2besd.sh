#!/bin/bash

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/lib
# smr_tool=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/software/smr-1.3.2-linux-x86_64/smr
smr_tool=/directflow/SCCGGroupShare/projects/angxue/software/smr-1.3.1-linux-x86_64/smr-1.3.1

CELL_TYPES_DIR="/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/output/main_analyses_NewGenotypes/preprocessing_caQTL_rawID"

# for cell_type in "${CELL_TYPES[@]}"
# do
cell_type=$1
echo "==============================================="
echo "Processing cell type: $cell_type"
echo "==============================================="

# Loop through chromosomes 1-22
# for chr in {1..22}
# do

chr=$2
CIS_CAQTL_FILE="$CELL_TYPES_DIR/$cell_type/matrixQTL/Chr${chr}_MatrixcaQTL.tsv"
BESD_FILE="$CELL_TYPES_DIR/$cell_type/BESD/Chr${chr}"
  
# Check if file exists
if [ ! -f "$CIS_CAQTL_FILE" ]; then
  echo "File for chromosome ${chr} in $cell_type does not exist. Skipping..."
  continue
fi
  
echo "Processing $cell_type - chromosome ${chr}..."
mkdir -p "$CELL_TYPES_DIR/$cell_type/BESD/"
${smr_tool} --eqtl-summary $CIS_CAQTL_FILE --matrix-eqtl-format --make-besd --out $BESD_FILE
  
if [ $? -eq 0 ]; then
  echo "Chromosome ${chr} for $cell_type completed successfully"
else
  echo "Error processing chromosome ${chr} for $cell_type"
fi

# done

echo "Completed processing for cell type: $cell_type"
echo ""
# done

echo "All cell types processed."