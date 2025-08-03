#!/usr/bin/env python3
"""
NK Cell UMAP Analysis with Genotype and Epigenetic Age Interactions

This script creates comprehensive UMAP visualizations for NK cells, showing:
1. Overall epigenetic age distribution
2. Genotype-stratified plots for significant SNP-epigenetic age interactions
3. Multi-panel visualizations showing interaction effects

Author: Peter C Allen
"""

import os
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import re
import glob
import numpy as np

print("Starting NK cell UMAP analysis with genotype interactions...")

# Configuration
fig_dir = "figures/5-epigeneticAge-Genotype-UMAPs"
os.makedirs(fig_dir, exist_ok=True)

# Load NK-specific h5ad file directly
nk_file_path = "output/20250630_celltype_subsets/repeat1_merged_dataset.h5ad"
print(f"Loading NK data from: {nk_file_path}")

adata_merged = sc.read_h5ad(nk_file_path)
print(f"Loaded NK data: {adata_merged.shape} (cells x features)")
print(f"Data appears to be raw counts - max value: {adata_merged.X.max()}, mean: {adata_merged.X.mean():.3f}")

# Load and attach epigenetic age data
base_dir = "output/20250630_celltype_subsets"
age_path = os.path.join(base_dir, "NK_epitrace_age.csv")

if os.path.exists(age_path):
    print(f"Loading epigenetic age data from: {age_path}")
    age = pd.read_csv(age_path, index_col=0)
    if "EpiTraceAge_iterative" in age.columns:
        # Match barcodes between the datasets
        common_barcodes = adata_merged.obs.index.intersection(age.index)
        print(f"Found {len(common_barcodes)} common barcodes between NK data and EpiTrace age data")
        
        if len(common_barcodes) > 0:
            adata_merged = adata_merged[common_barcodes].copy()
            adata_merged.obs = adata_merged.obs.join(age.loc[common_barcodes, ["EpiTraceAge_iterative"]], how="left")
            print(f"Successfully attached epigenetic age data to {adata_merged.shape[0]} cells")
        else:
            print("Warning: No common barcodes found between datasets!")
            # You might want to exit here or handle this case
    else:
        print(f"Column 'EpiTraceAge_iterative' not found in age file.")
else:
    print(f"Epigenetic age file not found at: {age_path}")
    print("Proceeding without epigenetic age data...")

print(f"\nFinal NK data shape: {adata_merged.shape}")
print(f"Available metadata columns: {adata_merged.obs.columns.tolist()}")

# --- 2. Process NK Data for UMAP ---

print("\nRunning analysis on NK object...")
print(f"Data summary: {adata_merged.shape} (cells x features)")
print(f"Data range before normalization: {adata_merged.X.min():.3f} to {adata_merged.X.max():.3f}")

# Apply normalization and log transformation (handles raw counts)
sc.pp.normalize_total(adata_merged)
sc.pp.log1p(adata_merged)
print(f"Data range after normalization/log1p: {adata_merged.X.min():.3f} to {adata_merged.X.max():.3f}")

sc.pp.pca(adata_merged)
sc.pp.neighbors(adata_merged)
sc.tl.umap(adata_merged)

adata_merged.write_h5ad("output/merged_NK_cells.h5ad")
print("Saved processed NK data to: output/merged_NK_cells.h5ad")

# --- 3. Plot UMAP by Epigenetic Age ---
if "EpiTraceAge_iterative" in adata_merged.obs.columns:
    print("\nGenerating UMAP plot colored by EpiTrace Age...")
    plt.figure(figsize=(7, 6))
    sc.pl.umap(
        adata_merged,
        color="EpiTraceAge_iterative",
        show=False,
        color_map="viridis",
        title="NK Cells Epigenetic Age",
        frameon=False,
        vmax=1.0
    )

    epiage_plot_path = os.path.join(fig_dir, "NK_merged_EpiTraceAge.png")
    plt.savefig(epiage_plot_path, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Successfully saved EpiTrace Age UMAP to: {epiage_plot_path}")
    
    # Generate hexagonal binned version
    print("Generating hexagonal binned UMAP plot colored by EpiTrace Age...")
    plt.figure(figsize=(7, 6))
    
    # Extract UMAP coordinates and EpiTrace Age values
    umap_coords = adata_merged.obsm['X_umap']
    epitrace_values = adata_merged.obs['EpiTraceAge_iterative'].values
    
    # Remove any NaN values
    valid_mask = ~np.isnan(epitrace_values)
    umap_coords_clean = umap_coords[valid_mask]
    epitrace_values_clean = epitrace_values[valid_mask]
    
    # Create hexbin plot
    plt.hexbin(umap_coords_clean[:, 0], umap_coords_clean[:, 1], 
               C=epitrace_values_clean, 
               gridsize=50, 
               cmap='viridis', 
               vmax=1.0,
               reduce_C_function=np.mean)
    plt.colorbar(label='EpiTrace Age')
    plt.title("NK Cells Epigenetic Age (Hexagonal Binned)")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")

    epiage_hex_plot_path = os.path.join(fig_dir, "NK_merged_EpiTraceAge_hexbin.png")
    plt.savefig(epiage_hex_plot_path, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Successfully saved hexagonal binned EpiTrace Age UMAP to: {epiage_hex_plot_path}")
else:
    print("\nSkipping EpiTrace Age plots - EpiTraceAge_iterative column not found")

# --- 4. Identify Significant Peaks and Plot UMAPs by Genotype/Epigenetic Age ---

print("\nIdentifying significant SNP-peak interactions for NK cells...")
interaction_results_file = "output/4-regression-results/combined_results_NK.csv"

if not os.path.exists(interaction_results_file):
    print(f"Interaction results file not found at: {interaction_results_file}")
    print("Skipping genotype analysis.")
else:
    nk_results = pd.read_csv(interaction_results_file)

    # Filter for significant results and select top 50
    significant_peaks = nk_results[
        (nk_results['peak'].isin(['chr4:55993465-55994658', 'chr6:31243369-31243947', 'chr7:5897477-5897880'])) &
        (nk_results['adj_pval'] < 0.05) &
        (nk_results['Estimate'] > 1)
    ].nsmallest(50, 'adj_pval')

    if significant_peaks.empty:
        print("No significant peaks found for NK cells with adj_pval < 0.05.")
    else:
        print(f"Found {len(significant_peaks)} significant peaks to plot by genotype and epigenetic age.")

        donor_id_column = 'donor_id'

        if donor_id_column not in adata_merged.obs.columns:
            print(f"Error: Donor ID column '{donor_id_column}' not found in AnnData object.")
            print(f"Available columns: {adata_merged.obs.columns.tolist()}")
        elif "EpiTraceAge_iterative" not in adata_merged.obs.columns:
            print("Warning: EpiTraceAge_iterative not found - skipping age-based analysis")
        else:
            for idx, row in significant_peaks.iterrows():
                peak = row["peak"]
                snp = row["snp"]
                print(f"\nProcessing SNP: {snp} for Peak: {peak}")

                # Ensure accessibility for this peak is in adata.obs
                # Check if peak is already in obs.columns, if not, extract from var_names
                if peak not in adata_merged.obs.columns:
                    if peak in adata_merged.var_names:
                        try:
                            peak_idx = list(adata_merged.var_names).index(peak)
                            adata_merged.obs[peak] = adata_merged.X[:, peak_idx].toarray().flatten()
                            print(f"Extracted accessibility for {peak} from var_names to obs.")
                        except Exception as e:
                            print(f"Could not extract accessibility for {peak}: {e}")
                            continue
                    else:
                        print(f"Peak {peak} not found in var_names or obs.")
                        continue
                else:
                    print(f"Peak {peak} already exists in obs.columns.")
                    
                # If peak exists in both var_names and obs, we need to ensure we use the obs version
                # Remove from var_names to avoid ambiguity if it exists in both
                if peak in adata_merged.var_names and peak in adata_merged.obs.columns:
                    print(f"Peak {peak} exists in both var_names and obs. Removing from var_names to avoid ambiguity.")
                    # Create a mask to keep all genes except this peak
                    keep_genes = adata_merged.var_names != peak
                    adata_merged = adata_merged[:, keep_genes].copy()

                match = re.match(r"^chr(\d+)", peak)
                if not match:
                    print(f"Could not extract chromosome from peak: {peak}")
                    continue

                chr_num = match.group(1)
                geno_file = f"data/genotype/wgs/TenK10K_TOB_ATAC_renamed_chr{chr_num}_common_variants.raw"

                if not os.path.exists(geno_file):
                    print(f"Genotype file not found: {geno_file}")
                    continue

                try:
                        with open(geno_file, 'r') as f:
                            header = f.readline().strip().split()

                        snp_cols = [i for i, col in enumerate(header) if col.startswith(snp)]
                        if not snp_cols:
                            print(f"SNP {snp} not found in genotype file.")
                            continue

                        iid_col_idx = header.index('IID')
                        snp_col_idx = snp_cols[0]

                        geno_df = pd.read_csv(geno_file, sep=r'\s+', usecols=[iid_col_idx, snp_col_idx], index_col='IID')
                        geno_df.columns = ['genotype']

                        # Parse alleles from SNP name
                        snp_parts = snp.split(':')
                        if len(snp_parts) == 4:
                            ref_allele = snp_parts[2]
                            alt_allele = snp_parts[3]
                        else:
                            ref_allele = 'REF'
                            alt_allele = 'ALT'

                        genotype_numeric_to_str = {
                            '0': ref_allele * 2,
                            '1': ref_allele + alt_allele,
                            '2': alt_allele * 2
                        }

                        # Merge genotype into AnnData object
                        adata_merged.obs[snp] = adata_merged.obs[donor_id_column].map(geno_df['genotype']).astype(pd.Int64Dtype()).astype(str)

                        mask = adata_merged.obs[snp].notna() & (adata_merged.obs[snp] != 'nan')
                        adata_plot = adata_merged[mask].copy()

                        # Bin epigenetic age
                        adata_plot.obs['age_group'] = pd.cut(
                            adata_plot.obs['EpiTraceAge_iterative'],
                            bins=[-np.inf, 0.5, np.inf],
                            labels=['young', 'old']
                        )

                        genotype_order = ['0', '1', '2']
                        adata_plot.obs[snp] = pd.Categorical(adata_plot.obs[snp], categories=genotype_order, ordered=True)
                        adata_plot.obs[f"{snp}_genotype_str"] = adata_plot.obs[snp].map(genotype_numeric_to_str)

                        # Create 6-panel plot: 2 rows (young, old) × 3 columns (genotypes)
                        fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)

                        for row_idx, age_group in enumerate(['young', 'old']):
                            for col_idx, gt in enumerate(genotype_order):
                                ax = axes[row_idx, col_idx]
                                cells = (
                                    (adata_plot.obs[snp] == gt) &
                                    (adata_plot.obs['age_group'] == age_group)
                                )
                                gt_label = genotype_numeric_to_str.get(gt, gt)
                                if cells.sum() == 0:
                                    ax.set_title(f"{gt_label} | {age_group}\nNo cells")
                                    ax.axis('off')
                                    continue
                                adata_sub = adata_plot[cells].copy()
                                # Ensure we're using the peak accessibility from obs
                                if peak in adata_sub.obs.columns:
                                    sc.pl.umap(
                                        adata_sub,
                                        color=peak,  # color by accessibility at the peak
                                        show=False,
                                        ax=ax,
                                        title=f"{gt_label} | {age_group}",
                                        frameon=False,
                                        color_map="viridis",
                                        size=20,
                                        vmax=1.0
                                    )
                                else:
                                    # If peak is not in obs, just plot without coloring
                                    sc.pl.umap(
                                        adata_sub,
                                        show=False,
                                        ax=ax,
                                        title=f"{gt_label} | {age_group}",
                                        frameon=False,
                                        size=20
                                    )

                        fig.suptitle(
                            f"NK Cells UMAP by Genotype and Epigenetic Age at {snp}\n"
                            f"Colored by Accessibility at {peak}",
                            fontsize=16
                        )
                        plot_path = os.path.join(fig_dir, f"NK_UMAP_{snp}_genotype_age_peak_{peak.replace(':', '_')}.png")
                        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                        plt.savefig(plot_path, bbox_inches="tight", dpi=300)
                        plt.close()
                        print(f"Successfully saved 6-panel genotype/age UMAP to: {plot_path}")

                        # Create hexagonal binned version of the 6-panel plot
                        fig_hex, axes_hex = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)

                        for row_idx, age_group in enumerate(['young', 'old']):
                            for col_idx, gt in enumerate(genotype_order):
                                ax_hex = axes_hex[row_idx, col_idx]
                                cells = (
                                    (adata_plot.obs[snp] == gt) &
                                    (adata_plot.obs['age_group'] == age_group)
                                )
                                gt_label = genotype_numeric_to_str.get(gt, gt)
                                if cells.sum() == 0:
                                    ax_hex.set_title(f"{gt_label} | {age_group}\nNo cells")
                                    ax_hex.axis('off')
                                    continue
                                    
                                adata_sub = adata_plot[cells].copy()
                                
                                # Extract UMAP coordinates for this subset
                                umap_coords_sub = adata_sub.obsm['X_umap']
                                
                                # Ensure we're using the peak accessibility from obs
                                if peak in adata_sub.obs.columns:
                                    peak_values = adata_sub.obs[peak].values
                                    # Remove any NaN values
                                    valid_mask = ~np.isnan(peak_values)
                                    if valid_mask.sum() > 0:  # Check if there are valid values
                                        umap_coords_clean = umap_coords_sub[valid_mask]
                                        peak_values_clean = peak_values[valid_mask]
                                        
                                        # Create hexbin plot
                                        hb = ax_hex.hexbin(umap_coords_clean[:, 0], umap_coords_clean[:, 1], 
                                                          C=peak_values_clean, 
                                                          gridsize=20, 
                                                          cmap='viridis', 
                                                          vmax=1.0,
                                                          reduce_C_function=np.mean)
                                    else:
                                        ax_hex.text(0.5, 0.5, 'No valid data', 
                                                   horizontalalignment='center',
                                                   verticalalignment='center',
                                                   transform=ax_hex.transAxes)
                                else:
                                    # If peak is not in obs, just plot density
                                    ax_hex.hexbin(umap_coords_sub[:, 0], umap_coords_sub[:, 1], 
                                                 gridsize=20, 
                                                 cmap='Blues')
                                
                                ax_hex.set_title(f"{gt_label} | {age_group}")
                                ax_hex.set_xlabel("UMAP1")
                                ax_hex.set_ylabel("UMAP2")

                        fig_hex.suptitle(
                            f"NK Cells UMAP by Genotype and Epigenetic Age at {snp} (Hexagonal Binned)\n"
                            f"Colored by Accessibility at {peak}",
                            fontsize=16
                        )
                        plot_hex_path = os.path.join(fig_dir, f"NK_UMAP_{snp}_genotype_age_peak_{peak.replace(':', '_')}_hexbin.png")
                        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                        plt.savefig(plot_hex_path, bbox_inches="tight", dpi=300)
                        plt.close()
                        print(f"Successfully saved 6-panel hexagonal binned genotype/age UMAP to: {plot_hex_path}")

                except Exception as e:
                    print(f"An error occurred while processing SNP {snp}: {e}")

print("\nScript finished.")
