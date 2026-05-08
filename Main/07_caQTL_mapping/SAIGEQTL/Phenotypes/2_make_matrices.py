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

# Read schema to filter for columns relevant to this chromosome
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
select_parts.extend([f"t1.{col}" for col in target_cols])

final_select = ", ".join(select_parts)
    
# Removed LEFT JOINs; strictly reads from source parquets
query = f"""
    COPY (
        SELECT {final_select} 
        FROM '{source_path}' AS t1
    ) 
    TO '{out_path}/chr{chr_num}_no_PCs.tsv' 
    (HEADER, DELIMITER '\t')
"""
    
duckdb.sql(query)

print("Finished!")