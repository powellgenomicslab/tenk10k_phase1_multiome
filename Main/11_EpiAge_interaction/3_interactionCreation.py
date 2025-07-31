import anndata as ad
import os
import pandas as pd
import glob

# Load the AnnData object
adata = ad.read_h5ad('output/20250630_celltype_subsets/repeat1_merged_dataset.h5ad')

# Extract the obs DataFrame
obs = adata.obs

# Display or use the obs DataFrame
print(obs)

obs.to_csv('output/20250630_celltype_subsets/repeat1_merged_dataset_meta.csv')

# For every *epitrace_age.csv file in the directory, read it and append to the obs DataFrame
epitrace_files = glob.glob('output/20250630_celltype_subsets/*epitrace_age.csv')
for epitrace_file in epitrace_files:
    # Read the epigenetic age csv
    epitrace_df = pd.read_csv(epitrace_file, index_col=0)

    # Extract the cell type from the file name
    cell_type = '_'.join(os.path.basename(epitrace_file).split('_')[:-2])

    # Merge with the obs DataFrame on the index
    obs_celltype = epitrace_df.merge(obs, left_index=True, right_index=True, how='left')

    # Calculate the median EpiTraceAge_iterative per donor_id into a new df with the donor_id as index
    median_epitrace = obs_celltype
    median_epitrace['EpiTraceAge_median'] = median_epitrace.groupby('donor_id')['EpiTraceAge_iterative'].transform('median')
    median_epitrace = median_epitrace.groupby('donor_id').agg({'EpiTraceAge_iterative': 'median'})
    median_epitrace.reset_index(inplace=True)
    median_epitrace.columns = ['individual', 'median_epiage']

    # Save the median_epitrace DataFrame to a text file
    output_path = f"data/20250630_epiages/{cell_type}_epiage.txt"
    median_epitrace.to_csv(output_path, sep='\t', index=False)

    # Calculate the standard deviation of EpiTraceAge_iterative per donor_id into a new df with the donor_id as index
    sd_epitrace = obs_celltype
    sd_epitrace['EpiTraceAge_std'] = sd_epitrace.groupby('donor_id')['EpiTraceAge_iterative'].transform('std')
    sd_epitrace = sd_epitrace.groupby('donor_id').agg({'EpiTraceAge_iterative': 'std'})
    sd_epitrace.reset_index(inplace=True)
    sd_epitrace.columns = ['individual', 'sd_epiage']

    # Save the median_epitrace DataFrame to a text file
    output_path = f"data/20250630_epiages/{cell_type}_epiageSD.txt"
    sd_epitrace.to_csv(output_path, sep='\t', index=False)
