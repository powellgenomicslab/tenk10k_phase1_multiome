# Import packages
import snapatac2 as snap
import pandas as pd
import os
import sys

# Convert the argument to integer
inputArgs = sys.argv[1:]
celltype_num = int(inputArgs[0])-1
print(celltype_num)

# Get celltype for this run
with open("celltype_names.txt", "r") as file:
    celltype_list = [line.strip() for line in file.readlines()]
celltype = celltype_list[celltype_num]
print(f"Starting {celltype}")

df_list = []
mat_path = "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/ATAC_Final/output/New_Peak_scanpy"
for file in sorted(os.listdir(mat_path)): 
    if not file.endswith("_1.h5ad"):
        continue
    
    print(f"Reading: {file}")
    data = snap.read(f"{mat_path}/{file}", backed = None)

    print(f"Filtering for '{celltype}'")
    data = data[data.obs["predicted.id"] == celltype]

    if data.n_obs == 0:
        print(f"No cells in {celltype}")
        continue

    print("Aggregating counts")
    data_agg = snap.tl.aggregate_X(data, groupby = "donor_id")
    df = pd.DataFrame.transpose(data_agg.to_df())
    df_list.append(df)

print("Concatenating")
if len(df_list) == 0:
    print("No dataframes to concatenate. Exiting.")
    sys.exit()
comb_df = pd.concat(df_list, axis = 1)

print("Saving")
comb_df.to_csv(f"PseudobulkMatrices/{celltype}.csv")
print("Finished!")
