import duckdb
import os
import sys
import re

inputArgs = sys.argv[1:]
celltype = inputArgs[0]
chr_num = int(inputArgs[1])
print(f"Processing {celltype} for chromosome {chr_num}")

# Define input and output directories
source_path = f"/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Parquets/{celltype}/*.parquet"
out_path = f"/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/FinalTSVs/{celltype}"
os.makedirs(out_path, exist_ok = True)

# Define paths to metadata files
donor_cov = "/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Covariates/donor_covariates.parquet" 
cell_cov = "/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Covariates/cell_covariates.parquet"

# Define path to PCs
pc_path = f"/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/PCA/{celltype}.csv" 

# Define columns of interest
donor_cols = ["sex", "geno_pc1", "geno_pc2", "geno_pc3", "geno_pc4", "geno_pc5", "geno_pc6", "age"]
cell_cols = ["repeat_num", "nCount_peaks"]
new_pc_cols = ["PC1", "PC2", "PC3", "PC4", "PC5"]

# Read schema
with open(f"/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/ValidPeaks/{celltype}.txt") as f:
    valid_peaks = [line.strip() for line in f if line.strip()]
all_columns = [re.sub(r'[:\-]', '_', x) for x in valid_peaks]

# Filter for chromosome being processed
target_cols = [c for c in all_columns if c.startswith(f"chr{chr_num}_")]
if not target_cols:
    print(f"No columns found for chromosome {chr_num}, exiting.")
    sys.exit(0)
print(f"Processing {len(target_cols)} columns")

# Initialise DuckDB
duckdb.sql("SET threads=8")
duckdb.sql("SET memory_limit='40GB'")

# Construct query
select_parts = ["t1.donor_id", "t1.barcode"]
select_parts.extend([f"t3.{col}" for col in cell_cols])
select_parts.extend([f"t2.{col}" for col in donor_cols])
select_parts.extend([f"t4.{col}" for col in new_pc_cols])
select_parts.extend([f"t1.{col}" for col in target_cols])

final_select = ", ".join(select_parts)
    
query = f"""
    COPY (
        SELECT {final_select} 
        FROM '{source_path}' AS t1
        LEFT JOIN '{donor_cov}' AS t2 ON t1.donor_id = t2.donor_id
        LEFT JOIN '{cell_cov}' AS t3 ON t1.barcode = t3.barcode
        LEFT JOIN '{pc_path}' AS t4 ON t1.barcode = t4.barcode
    ) 
    TO '{out_path}/chr{chr_num}.tsv' 
    (HEADER, DELIMITER '\t')
"""
    
duckdb.sql(query)

print("Finished!")