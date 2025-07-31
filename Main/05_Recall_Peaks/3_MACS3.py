# Import packages
import snapatac2 as snap
import sys

# Extract celltype from argument
inputArgs = sys.argv[1:]
cur_celltype = inputArgs[0]
print(f"Cell type: {cur_celltype}")

# Obtain fragment file
fragment_file = f"/g/data/ei56/od8037/TenK10K/PeakCalling/CelltypeFragments/AllCombinedFragments/{cur_celltype}.tsv.gz"
print(f"Path to fragment file: {fragment_file}")

# Read in data
print("Importing data")
data = snap.pp.import_fragments(
    fragment_file, 
    chrom_sizes = snap.genome.hg38,
    sorted_by_barcode = False
)

# Call peaks using MACS3
print("Calling peaks")
snap.tl.macs3(data, inplace = True)

# Save the peaks
print("Saving peaks")
peaks_df = data.uns["macs3_pseudobulk"]
print(type(peaks_df))
peaks_df.to_csv(f"/g/data/ei56/od8037/TenK10K/PeakCalling/CelltypePeaks/{cur_celltype}.csv")
print("Finished!")