import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binomtest
import gzip
import re

input_dir = "output/4-regression-results"
file_paths = glob.glob(os.path.join(input_dir, "interaction_results_*_chr[0-9]*.csv"))
cell_types = list({re.sub(r".*_(results_|SNP)(.*)_chr\d+\.csv", r"\2", os.path.basename(fp)) for fp in file_paths})

for cell_type in cell_types:
    # Filter file paths for the current cell type
    cell_type_files = [fp for fp in file_paths if re.search(rf"interaction_results_{cell_type}_chr\d+\.csv", os.path.basename(fp))]

    # Import and combine all files into a single DataFrame
    combined_results = pd.concat([pd.read_csv(f) for f in cell_type_files], ignore_index=True)

    # Sort by adj_pval
    combined_results = combined_results.sort_values('adj_pval')

    # Save the combined results to a CSV file
    output_path = f"output/4-regression-results/combined_results_{cell_type}.csv"
    combined_results.to_csv(output_path, index=False)

# Path pattern for the CSV files
file_pattern = 'output/4-regression-results/combined*.csv'

# List to hold individual DataFrames
dfs = []

for filepath in glob.glob(file_pattern):
    # Extract filename
    filename = os.path.basename(filepath)
    # Extract cell type name
    cell_type = filename.replace('combined_results_', '').replace('.csv', '')
    # Read CSV
    df = pd.read_csv(filepath)
    # Add cell type column
    df['cell_type'] = cell_type
    dfs.append(df)

# Combine all DataFrames
combined_df = pd.concat(dfs, ignore_index=True)

# Filter for significant results
filtered_df = combined_df[combined_df['adj_pval'] < 0.05]

# Load color palette
palette_df = pd.read_csv("data/colour_palette_table.tsv", sep="\t")
palette_df.columns = palette_df.columns.str.strip()
palette_df = palette_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
palette_dict = dict(zip(palette_df["cell_type"], palette_df["color"]))

# Map cell_type_plot from palette_df to filtered_df
filtered_df = filtered_df.merge(
    palette_df[["wg2_scpred_prediction", "cell_type"]],
    left_on="cell_type",
    right_on="wg2_scpred_prediction",
    how="left",
)
filtered_df.rename(columns={"cell_type_y": "cell_type_plot"}, inplace=True)

# Count number of SNPs per cell type
snp_counts = filtered_df.groupby("cell_type_plot")["snp"].nunique().reset_index()

# Map colors to cell types
bar_colors = [palette_dict.get(ct, "#333333") for ct in snp_counts["cell_type_plot"]]

plt.figure(figsize=(10, 6))
bars = plt.bar(snp_counts["cell_type_plot"], snp_counts["snp"], color=bar_colors)

# Add counts on top of bars
# Desired cell type order
desired_order = [
    "CD4 TCM", "CD4 Naive", "CD4 TEM", "CD4 CTL", "Treg", "CD4 Proliferating",
    "gdT", "MAIT", "dnT", "ILC", "CD8 TEM", "CD8 Naive", "CD8 TCM", "CD8 Proliferating",
    "NK", "NK CD56bright", "NK Proliferating", "B naive", "B intermediate", "B memory",
    "Plasmablast", "CD14 Mono", "CD16 Mono", "cDC2", "pDC", "cDC1", "ASDC", "HSPC"
]

# Ensure all desired cell types are present, fill missing with 0
snp_counts_ordered = pd.DataFrame({
    "cell_type_plot": desired_order,
    "snp": [snp_counts.set_index("cell_type_plot").get("snp", {}).get(ct, 0) for ct in desired_order]
})

# Replace "B naive" with "B Naive" in cell_type_plot column
snp_counts_ordered["cell_type_plot"] = snp_counts_ordered["cell_type_plot"].replace("B naive", "B Naive")

# Map colors to cell types in desired order
bar_colors_ordered = [palette_dict.get(ct, "#333333") for ct in desired_order]

plt.figure(figsize=(14, 7))
bars = plt.bar(snp_counts_ordered["cell_type_plot"], snp_counts_ordered["snp"], color=bar_colors_ordered)

# Remove extra spacing by adjusting x-axis limits
plt.xlim(-0.5, len(snp_counts_ordered["cell_type_plot"]) - 0.5)

# Remove top and right spines for a classic axis-only look
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# Add counts on top of bars
for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        height,
        f"{int(height)}",
        ha="center",
        va="bottom",
        fontsize=12,
    )

plt.xlabel("Cell Type", fontsize=14)
plt.ylabel("Number of Significant Interactions", fontsize=14)
plt.xticks(rotation=90, fontsize=12)
plt.tight_layout()
plt.savefig("figures/3-regression/numberSig-caSNPs-by-celltype.png", dpi=300, bbox_inches="tight")
plt.close()

# Calculate unique and shared SNPs per cell type
snp_celltype = filtered_df.groupby("snp")["cell_type_plot"].nunique().reset_index()
snp_celltype.columns = ["snp", "n_celltypes"]

# Merge with filtered_df to annotate SNPs as unique/shared
filtered_df = filtered_df.merge(snp_celltype, on="snp", how="left")
filtered_df["snp_type"] = np.where(filtered_df["n_celltypes"] == 1, "unique", "shared")

# Count unique and shared SNPs per cell type
snp_type_counts = (
    filtered_df.groupby(["cell_type_plot", "snp_type"])["snp"]
    .nunique()
    .unstack(fill_value=0)
    .reindex(desired_order)
    .fillna(0)
)

# Ensure both columns exist
for col in ["unique", "shared"]:
    if col not in snp_type_counts.columns:
        snp_type_counts[col] = 0

# Calculate total SNPs per cell type (same as previous plot)
snp_counts_ordered = pd.DataFrame({
    "cell_type_plot": desired_order,
    "snp": [snp_type_counts.loc[ct, ["unique", "shared"]].sum() if ct in snp_type_counts.index else 0 for ct in desired_order]
})

# Plot stacked barplot
plt.figure(figsize=(14, 7))
bottom_vals = snp_type_counts["shared"].values
bars_shared = plt.bar(
    snp_type_counts.index,
    snp_type_counts["shared"],
    color="#D88C9A",
    label="Shared SNPs"
)
bars_unique = plt.bar(
    snp_type_counts.index,
    snp_type_counts["unique"],
    bottom=bottom_vals,
    color="#7BAFAF",
    label="Unique SNPs"
)

# Add total SNP counts on top of bars
for i, ct in enumerate(snp_type_counts.index):
    total = snp_counts_ordered[snp_counts_ordered["cell_type_plot"] == ct]["snp"].values[0]
    plt.text(
        i,
        snp_type_counts.loc[ct, ["unique", "shared"]].sum(),
        f"{int(total)}",
        ha="center",
        va="bottom",
        fontsize=12,
    )

plt.xlabel("Cell Type", fontsize=14)
plt.ylabel("Number of Significant Interactions", fontsize=14)
plt.xticks(rotation=90, fontsize=12)
plt.xlim(-0.5, len(snp_type_counts.index) - 0.5)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

plt.legend(title="SNP Type", fontsize=12, title_fontsize=13, loc="upper right")
plt.tight_layout()
plt.savefig("figures/3-regression/numberSig-caSNPs-unique-shared-by-celltype.png", dpi=300, bbox_inches="tight")
plt.close()

# Percentage of positive estimates vs negative estimates
positive_count = (filtered_df['Estimate'] > 0).sum()
negative_count = (filtered_df['Estimate'] < 0).sum()
total_count = positive_count + negative_count

if total_count > 0:
    positive_percentage = (positive_count / total_count) * 100
    negative_percentage = (negative_count / total_count) * 100
else:
    positive_percentage = 0
    negative_percentage = 0


# Perform a binomial test to see if positives are more likely than negatives
# Null hypothesis: probability of positive = 0.5

p_value = binomtest(positive_count, n=total_count, p=0.5, alternative='greater')
print(f"Positive estimates: {positive_count} ({positive_percentage:.2f}%)")
print(f"Negative estimates: {negative_count} ({negative_percentage:.2f}%)")
print(f"Binomial test p-value (positive > negative): {p_value.pvalue:.4g}")

# Identify the strongest positive and negative estimates with MAF greater than 0.12
filtered_df['cell_type'] = filtered_df['cell_type_x']

strongest_positive = filtered_df[(filtered_df['Estimate'] > 0) & (filtered_df['MAF'] > 0.12)].nlargest(50, 'Estimate')
# strongest_negative = filtered_df[(filtered_df['Estimate'] < 0) & (filtered_df['MAF'] > 0.12)].nsmallest(5, 'Estimate')
# strongest_positive = strongest_positive[strongest_positive['peak'].str.contains('chr6:31243369')]
# GWAS hits table (replace this with your actual DataFrame)
gwas_hits_table = pd.DataFrame({
    "cell_type": [
        "B_intermediate", "B_intermediate", "B_naive", "B_naive", "B_naive", "B_naive",
        "CD14_Mono", "CD4_TCM", "CD4_TCM", "CD4_TCM", "CD4_TCM",
        "B_naive", "CD4_TCM", "NK", "NK"
    ],
    "trait_name": [
        "breastca_GCST004988", "asthma", "CD_EAS_EUR", "CD_EAS_EUR", "IBD_EAS_EUR", "IBD_EAS_EUR",
        "asthma", "CD_EAS_EUR", "CD_EAS_EUR", "IBD_EAS_EUR", "IBD_EAS_EUR",
        "prostateca_GCST90274713", "IBD_EAS_EUR", "CD_EAS_EUR", "IBD_EAS_EUR"
    ],
    "Chr": [1, 5, 16, 16, 16, 16, 16, 6, 6, 6, 6, 7, 21, 21, 21],
    "Peak": [
        "chr1:121365304-121366060", "chr5:132378406-132379096", "chr16:50337619-50338016", "chr16:50337619-50338016",
        "chr16:50337619-50338016", "chr16:50337619-50338016", "chr16:27332015-27332584",
        "chr6:167020258-167020914", "chr6:167020258-167020914", "chr6:167020258-167020914", "chr6:167020258-167020914",
        "chr7:107560909-107561351", "chr21:44143489-44143752", "chr21:44143489-44143752", "chr21:44143489-44143752"
    ],
    "Gene": [
        None, None, None, None, None, None, None, None, None, None, None,
        "COG5", "TRAPPC10", "TRAPPC10", None
    ]
})

# Standardize cell type column names to match filtered_df
gwas_hits_table["cell_type"] = gwas_hits_table["cell_type"].str.replace("_", " ")

# Filter filtered_df to match cell_type and peak in gwas_hits_table
gwas_hits = filtered_df[
    filtered_df.apply(
        lambda row: ((row["cell_type"].replace("_", " ") in gwas_hits_table["cell_type"].values) and
                     (row["peak"] in gwas_hits_table[gwas_hits_table["cell_type"] == row["cell_type"].replace("_", " ")]["Peak"].values)),
        axis=1
    )
]

# def load_and_plot_interactions(
#     variants_df,
#     output_dir="figures/3-regression/caSNPxPeak-boxplots",
# ):
#     """Create boxplots for SNP x EpiAge interactions on chromatin accessibility, with t-test significance."""
#     os.makedirs(output_dir, exist_ok=True)
#     print(f"Output directory: {output_dir}")

#     from scipy.stats import ttest_ind

#     for i, (idx, row) in enumerate(variants_df.iterrows()):
#         print(f"\nProcessing variant {i+1}/{len(variants_df)}")
#         cell_type = row["cell_type_x"]
#         peak = row["peak"]
#         snp = row["snp"]
#         estimate = row["Estimate"]
#         print(f"Cell type: {cell_type}, Peak: {peak}, SNP: {snp}")

#         match = re.match(r"^chr(\d+)", peak)
#         if match:
#             chr_num = match.group(1)
#         else:
#             print(f"Could not extract chromosome number from peak: {peak}")
#             continue

#         try:
#             # Load expression data
#             expr_file = (
#                 f"data/expressionBeds/{cell_type}/ExpressionBeds/chr{chr_num}.bed.gz"
#             )
#             print(f"Loading expression data from {expr_file}")
#             expr_df = pd.read_csv(expr_file, sep="\t", compression="gzip", index_col=3)

#             # Load covariates and epiAge
#             cov_file = f"data/expressionBeds/{cell_type}/covariates.txt"
#             print(f"Loading covariates from {cov_file}")
#             cov_df = pd.read_csv(cov_file, index_col=0).T
#             cov_df.index = cov_df.index.str.replace("X", "")

#             epiage_file = f"data/20250630_epiages/{cell_type}_epiage.txt"
#             print(f"Loading epiAge data from {epiage_file}")
#             epiage_df = pd.read_csv(epiage_file, sep="\t").dropna(
#                 subset=["median_epiage"]
#             )
#             epiage_df["epiAge_binary"] = np.where(
#                 epiage_df["median_epiage"] > 0.5, "Old", "Young"
#             )

#             # Load genotype data efficiently - only read header first to find SNP column
#             geno_file = f"data/genotype/wgs/TenK10K_TOB_ATAC_renamed_chr{chr_num}_common_variants.raw"
#             print(f"Loading genotype data from {geno_file}")
            
#             # Read header to find which column contains our SNP
#             with open(geno_file, 'r') as f:
#                 header = f.readline().strip().split()
            
#             # Find SNP column (may have _A, _G, _T, _C suffix)
#             snp_cols = [i for i, col in enumerate(header) if col.startswith(snp)]
#             if not snp_cols:
#                 print(f"SNP {snp} not found in genotype file headers")
#                 continue
            
#             # Load only IID column and the SNP column
#             iid_col = header.index('IID')
#             snp_col = snp_cols[0]  # Take first match
#             usecols = [iid_col, snp_col]
            
#             geno_df = pd.read_csv(geno_file, sep=r"\s+", usecols=usecols, index_col='IID')
#             geno_df.columns = [header[snp_col].split('_')[0]]  # Clean column name

#             if peak not in expr_df.index:
#                 print(f"Peak {peak} not found in expression data for {cell_type}.")
#                 continue
            
#             # Use the cleaned column name (should be the SNP name)
#             snp_col_name = geno_df.columns[0]
#             if snp_col_name != snp:
#                 print(f"Expected SNP {snp} but got column {snp_col_name}")

#             print("Preparing data for plotting...")
#             peak_expr = expr_df.loc[peak].iloc[
#                 3:
#             ]  # Skip first 3 columns (chr, start, end)
#             snp_geno = geno_df.iloc[:, 0]  # Get the first (and only) column

#             # Ensure indices match between peak_expr and snp_geno
#             common_ids = peak_expr.index.intersection(snp_geno.index)
#             if len(common_ids) == 0:
#                 print(
#                     f"No overlapping individuals between expression and genotype for {cell_type}, {peak}, {snp}."
#                 )
#                 continue

#             peak_expr = peak_expr.loc[common_ids]
#             snp_geno = snp_geno.loc[common_ids]

#             # Merge data
#             plot_data = pd.DataFrame(
#                 {
#                     "ID": peak_expr.index,
#                     "expr": peak_expr.values,
#                     "genotype": snp_geno.astype(str),
#                 }
#             )

#             if "individual" not in epiage_df.columns:
#                 print(f"'individual' column missing in epiage file for {cell_type}.")
#                 continue

#             print("Merging with epiAge data...")
#             plot_data = plot_data.merge(
#                 epiage_df[["individual", "epiAge_binary"]],
#                 left_on="ID",
#                 right_on="individual",
#                 how="left",
#             )
#             plot_data = plot_data.dropna()

#             # Remove rows where genotype is NA before plotting
#             plot_data = plot_data[plot_data["genotype"].notna()]
#             plot_data = plot_data[plot_data["genotype"] != "nan"]

#             # Filter out chromatin accessibility values above threshold
#             plot_data = plot_data[plot_data["expr"] <= 10]

#             # Order epiAge_binary: Young first, then Old
#             plot_data["epiAge_binary"] = pd.Categorical(
#                 plot_data["epiAge_binary"], categories=["Young", "Old"], ordered=True
#             )

#             print("Plotting boxplot...")
#             plt.figure(figsize=(16, 6))
#             # Set classic theme: white background, no grid, black axes
#             plt.style.use("default")
#             ax = plt.gca()
#             ax.set_facecolor("white")
#             for spine in ax.spines.values():
#                 spine.set_visible(True)
#                 spine.set_color("black")
#             ax.grid(False)

#             # Draw boxplot with dodge for space between Young and Old
#             sns.boxplot(
#                 data=plot_data,
#                 x="genotype",
#                 y="expr",
#                 hue="epiAge_binary",
#                 palette={"Young": "#7FB2FF", "Old": "#D6C7A1"},
#                 order=sorted(plot_data["genotype"].unique()),
#                 hue_order=["Young", "Old"],
#                 dodge=True,
#                 width=0.6,
#                 fliersize=2,
#                 linewidth=1.5,
#                 boxprops=dict(edgecolor="black"),
#                 medianprops=dict(color="black"),
#                 whiskerprops=dict(color="black"),
#                 capprops=dict(color="black"),
#             )

#             # Add sample counts and t-test significance
#             xticks = ax.get_xticks()
#             unique_genotypes = sorted(plot_data["genotype"].unique())
#             for i_tick, genotype in enumerate(unique_genotypes):
#                 for j, age_group in enumerate(["Young", "Old"]):
#                     count = len(
#                         plot_data[
#                             (plot_data["genotype"] == genotype)
#                             & (plot_data["epiAge_binary"] == age_group)
#                         ]
#                     )
#                     if count > 0:
#                         # Place count above each box
#                         x_pos = i_tick - 0.2 + 0.4 * j  # Young left, Old right
#                         ax.text(
#                             x_pos,
#                             plot_data["expr"].max() + 0.1,
#                             str(count),
#                             ha="center",
#                             va="bottom",
#                             fontsize=11,
#                         )
#                 # T-test between Young and Old for this genotype
#                 expr_young = plot_data[
#                     (plot_data["genotype"] == genotype) & (plot_data["epiAge_binary"] == "Young")
#                 ]["expr"]
#                 expr_old = plot_data[
#                     (plot_data["genotype"] == genotype) & (plot_data["epiAge_binary"] == "Old")
#                 ]["expr"]
#                 # Ensure numeric dtype for t-test
#                 expr_young = pd.to_numeric(expr_young, errors="coerce")
#                 expr_old = pd.to_numeric(expr_old, errors="coerce")
#                 expr_young = expr_young.dropna()
#                 expr_old = expr_old.dropna()
#                 if len(expr_young) > 1 and len(expr_old) > 1:
#                     t_stat, p_val = ttest_ind(expr_young, expr_old, equal_var=False)
#                     # Annotate significance above the boxes
#                     x_center = i_tick
#                     y_max = max(expr_young.max() if not expr_young.empty else 0,
#                                 expr_old.max() if not expr_old.empty else 0)
#                     y_annot = y_max + 0.3
#                     if p_val < 0.001:
#                         sig = "***"
#                     elif p_val < 0.01:
#                         sig = "**"
#                     elif p_val < 0.05:
#                         sig = "*"
#                     else:
#                         sig = "n.s."
#                     ax.text(
#                         x_center,
#                         y_annot,
#                         sig,
#                         ha="center",
#                         va="bottom",
#                         fontsize=15,
#                         color="black",
#                         fontweight="bold",
#                     )

#             # Axis labels and title with subtitle
#             plt.xlabel("Genotype", fontsize=15)
#             plt.ylabel("Chromatin Accessibility", fontsize=15)
            
#             # Add title with SNP: chr#:bp
#             # Extract chromosome and position from SNP name (e.g., '1_200166672_C_G' -> 'chr 1:200166672')
#             snp_parts = snp.split(":")
#             if len(snp_parts) >= 2:
#                 snp_chr_bp = f"chr{snp_parts[0]}:{snp_parts[1]}"
#             else:
#                 snp_chr_bp = snp  # fallback
#             plt.title(f"SNP: {snp_chr_bp}", fontsize=16)

#             plt.xticks(fontsize=14)
#             plt.yticks(fontsize=14)

#             # Remove top and right spines for classic look
#             ax.spines["top"].set_visible(False)
#             ax.spines["right"].set_visible(False)

#             # Place legend outside the plot
#             handles, labels = ax.get_legend_handles_labels()
#             ax.legend(
#                 handles,
#                 labels,
#                 title="EpiAge Group",
#                 title_fontsize=14,
#                 fontsize=13,
#                 loc="center left",
#                 bbox_to_anchor=(1.02, 0.5),
#                 borderaxespad=0,
#                 frameon=False,
#             )

#             plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space for legend

#             # Add pos_ or neg_ prefix to filename based on Estimate
#             prefix = "pos_" if estimate > 0 else "neg_"
#             filename = f"{i+1}_{cell_type}_{peak}_caSNPxPeak_boxplot.png"
#             plot_path = f"figures/3-regression/gwas_boxplots/{filename}"
#             print(f"Saving plot to {plot_path}")
#             plt.savefig(plot_path, dpi=300, bbox_inches="tight")
#             plt.close()

#             # Explicitly delete large objects and collect garbage to free memory
#             del expr_df, cov_df, epiage_df, geno_df, peak_expr, snp_geno, plot_data
#             import gc
#             gc.collect()
#             print("Memory cleaned up for this variant.")

#         except Exception as e:
#             import traceback
#             print(f"Error plotting {cell_type} {snp}: {e}")
#             traceback.print_exc()
#             continue

# # Plot strongest interactions
# # print("\nPlotting strongest positive interactions...")
# # load_and_plot_interactions(strongest_positive)

# # Plot strongest interactions
# print("\nPlotting GWAS hits...")
# load_and_plot_interactions(gwas_hits)

# # print("Plotting strongest negative interactions...")  
# # load_and_plot_interactions(strongest_negative)
