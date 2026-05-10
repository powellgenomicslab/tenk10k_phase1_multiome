
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--chr', type=int, required=True, help='Chromosome number (1-22)')
args = parser.parse_args()

chr = args.chr
input_dir = "/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Phenotypes/FinalTSVs/"
output_dir = "/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Phenotypes/FinalTSVs_merged/"
region_dir = "/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Common/Regions/"
region_output_dir = "/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Regions/"

for chr in range(1, 23):
    dfs = []
    for celltype in sorted(os.listdir(input_dir)):
        celltype_path = os.path.join(input_dir, celltype)
        if os.path.isdir(celltype_path):
            file_path = os.path.join(celltype_path, f"chr{chr}.tsv")
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, sep='\t')
                df['celltype'] = celltype
                dfs.append(df)
                print(f"{celltype} added")

    merged_df = pd.concat(dfs, ignore_index=True)
    merged_df = merged_df.dropna(axis=1)  # Drop columns with any NaN
    merged_df.to_csv(os.path.join(output_dir, f"chr{chr}.tsv"), sep='\t', index=False)
    print(f"Merged {len(dfs)} files for chr{chr}")


    merged_df = pd.read_csv(os.path.join(output_dir, f"chr{chr}.tsv"), sep='\t')
    region_dfs = []
    for celltype in sorted(os.listdir(input_dir)):
        celltype_path = os.path.join(region_dir, celltype)
        if os.path.isdir(celltype_path):
            file_path = os.path.join(celltype_path, f"chr{chr}.csv")
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, sep='\t', header=None)
                region_dfs.append(df)
                print(f"{celltype} added")
    merged_region_df = pd.concat(region_dfs, ignore_index=True)
    merged_region_df.columns = ['peak', 'chr', 'start', 'end']
    merged_region_df.drop_duplicates(subset='peak', inplace=True)
    merged_region_df = merged_region_df[merged_region_df['peak'].isin(merged_df.columns)]
    merged_region_df.to_csv(os.path.join(region_output_dir, f"chr{chr}.csv"), sep='\t', index=False, header=False)