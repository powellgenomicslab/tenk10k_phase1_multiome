#!/usr/bin/env python3
"""
NK Cell UMAP Plotting Functions

This module contains functions for creating UMAP visualizations of NK cells
with genotype and epigenetic age interactions.

Usage:
    from src.nk_umap_plotting import plot_nk_genotype_interactions
    
    # Run the plotting function
    plot_nk_genotype_interactions(
        adata_merged, 
        interaction_results_file="output/4-regression-results/combined_results_NK.csv",
        fig_dir="figures/5-epigeneticAge-Genotype-UMAPs"
    )

Author: Peter C Allen
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import scanpy as sc


def plot_nk_epigenetic_age(adata_merged, 
                          fig_dir="figures/5-epigeneticAge-Genotype-UMAPs"):
    """
    Plot NK cell UMAP visualizations colored by epigenetic age.
    
    Creates both regular scatter plots and hexagonal binned versions.
    
    Parameters:
    -----------
    adata_merged : AnnData
        Processed AnnData object with NK cells, UMAP coordinates, and EpiTraceAge_iterative
    fig_dir : str
        Output directory for figures
        
    Returns:
    --------
    None
        Saves plots to fig_dir
    """
    
    # Create output directory
    os.makedirs(fig_dir, exist_ok=True)
    
    if "EpiTraceAge_iterative" not in adata_merged.obs.columns:
        print("Warning: EpiTraceAge_iterative not found - skipping epigenetic age plots")
        return
    
    print("\nGenerating UMAP plot colored by EpiTrace Age...")
    plt.rcParams.update({'font.size': 20})  # Set global font size
    plt.figure(figsize=(7, 6), dpi=300)
    sc.pl.umap(
        adata_merged,
        color="EpiTraceAge_iterative",
        show=False,
        color_map="Spectral_r",
        title="NK Cells Epigenetic Age",
        frameon=False,
        vmax=1.0,
        size=20  # Increased cell size
    )

    epiage_plot_path = os.path.join(fig_dir, "NK_merged_EpiTraceAge.png")
    plt.savefig(epiage_plot_path, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Successfully saved EpiTrace Age UMAP to: {epiage_plot_path}")
    
    # Generate hexagonal binned version
    print("Generating hexagonal binned UMAP plot colored by EpiTrace Age...")
    plt.rcParams.update({'font.size': 20})  # Ensure font size for hexbin plot
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
               gridsize=20, 
               cmap='Spectral_r', 
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


def plot_nk_genotype_interactions(adata_merged, 
                                 interaction_results_file="output/4-regression-results/combined_results_NK.csv",
                                 fig_dir="figures/5-epigeneticAge-Genotype-UMAPs",
                                 target_peaks=['chr4:55993465-55994658', 'chr6:31243369-31243947', 'chr7:5897477-5897880'],
                                 adj_pval_threshold=0.05,
                                 estimate_threshold=1,
                                 max_peaks=50):
    """
    Plot NK cell UMAP visualizations by genotype and epigenetic age interactions.
    
    Parameters:
    -----------
    adata_merged : AnnData
        Processed AnnData object with NK cells, UMAP coordinates, and metadata
    interaction_results_file : str
        Path to CSV file containing interaction results
    fig_dir : str
        Output directory for figures
    target_peaks : list
        List of specific peaks to filter for
    adj_pval_threshold : float
        Adjusted p-value threshold for significance
    estimate_threshold : float
        Minimum estimate threshold
    max_peaks : int
        Maximum number of peaks to process
        
    Returns:
    --------
    None
        Saves plots to fig_dir
    """
    
    # Create output directory
    os.makedirs(fig_dir, exist_ok=True)
    
    print(f"\nIdentifying significant SNP-peak interactions for NK cells...")
    
    if not os.path.exists(interaction_results_file):
        print(f"Interaction results file not found at: {interaction_results_file}")
        print("Skipping genotype analysis.")
        return
    
    nk_results = pd.read_csv(interaction_results_file)

    # Filter for significant results and select top peaks
    significant_peaks = nk_results[
        (nk_results['peak'].isin(target_peaks)) &
        (nk_results['adj_pval'] < adj_pval_threshold) &
        (nk_results['Estimate'] > estimate_threshold)
    ].nsmallest(max_peaks, 'adj_pval')

    if significant_peaks.empty:
        print(f"No significant peaks found for NK cells with adj_pval < {adj_pval_threshold}.")
        return
    else:
        print(f"Found {len(significant_peaks)} significant peaks to plot by genotype and epigenetic age.")

    donor_id_column = 'donor_id'

    if donor_id_column not in adata_merged.obs.columns:
        print(f"Error: Donor ID column '{donor_id_column}' not found in AnnData object.")
        print(f"Available columns: {adata_merged.obs.columns.tolist()}")
        return
    elif "EpiTraceAge_iterative" not in adata_merged.obs.columns:
        print("Warning: EpiTraceAge_iterative not found - skipping age-based analysis")
        return
    
    for idx, row in significant_peaks.iterrows():
        peak = row["peak"]
        snp = row["snp"]
        print(f"\nProcessing SNP: {snp} for Peak: {peak}")

        # Ensure accessibility for this peak is in adata.obs
        if peak not in adata_merged.obs.columns:
            if peak in adata_merged.var_names:
                try:
                    peak_idx = list(adata_merged.var_names).index(peak)
                    adata_merged.obs[peak] = adata_merged.X[:, peak_idx].flatten()
                    print(f"Extracted accessibility for {peak} from var_names to obs.")
                except Exception as e:
                    print(f"Could not extract accessibility for {peak}: {e}")
                    continue
            else:
                print(f"Peak {peak} not found in var_names or obs.")
                continue
        else:
            print(f"Peak {peak} already exists in obs.columns.")
            
        if peak in adata_merged.var_names and peak in adata_merged.obs.columns:
            print(f"Peak {peak} exists in both var_names and obs. Will use obs version for plotting.")

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

            # Merge genotype into AnnData object metadata
            adata_merged.obs[snp] = adata_merged.obs[donor_id_column].map(geno_df['genotype']).astype(pd.Int64Dtype()).astype(str)

            # Work with indices instead of copying data to avoid memory issues
            mask = adata_merged.obs[snp].notna() & (adata_merged.obs[snp] != 'nan')
            
            # Create a temporary DataFrame with just the metadata we need
            plot_obs = adata_merged.obs[mask].copy()
            umap_coords = adata_merged.obsm['X_umap'][mask]

            # Bin epigenetic age
            plot_obs['age_group'] = pd.cut(
                plot_obs['EpiTraceAge_iterative'],
                bins=[-np.inf, 0.5, np.inf],
                labels=['young', 'old']
            )

            genotype_order = ['0', '1', '2']
            plot_obs[snp] = pd.Categorical(plot_obs[snp], categories=genotype_order, ordered=True)
            plot_obs[f"{snp}_genotype_str"] = plot_obs[snp].map(genotype_numeric_to_str)

            # Create 6-panel plot: 2 rows (young, old) × 3 columns (genotypes)
            fig, axes = plt.subplots(2, 3, figsize=(18, 14), sharex=True, sharey=True)

            # Calculate global min/max values for consistent color scaling across all scatter plots
            if peak in plot_obs.columns:
                scatter_global_vmin = plot_obs[peak].min()
                scatter_global_vmax = plot_obs[peak].max()
            else:
                scatter_global_vmin, scatter_global_vmax = None, None

            # Set global font size
            plt.rcParams.update({'font.size': 20})

            for row_idx, age_group in enumerate(['young', 'old']):
                for col_idx, gt in enumerate(genotype_order):
                    ax = axes[row_idx, col_idx]
                    # Use DataFrame indices instead of copying AnnData objects
                    cells_mask = (
                        (plot_obs[snp] == gt) &
                        (plot_obs['age_group'] == age_group)
                    )
                    gt_label = genotype_numeric_to_str.get(gt, gt)
                    if cells_mask.sum() == 0:
                        ax.set_title(f"{gt_label} | {age_group}\nNo cells", fontsize=20)
                        ax.axis('off')
                        continue
                    
                    # Get UMAP coordinates and peak values for this subset
                    sub_umap = umap_coords[cells_mask]
                    
                    # Plot using matplotlib directly to avoid AnnData copying
                    if peak in plot_obs.columns:
                        peak_values = plot_obs.loc[cells_mask, peak].values
                        
                        # Sort cells by chromatin accessibility (low to high) so high values are plotted on top
                        sort_idx = np.argsort(peak_values)
                        sub_umap_sorted = sub_umap[sort_idx]
                        peak_values_sorted = peak_values[sort_idx]
                        
                        scatter = ax.scatter(sub_umap_sorted[:, 0], sub_umap_sorted[:, 1], 
                                           c=peak_values_sorted, cmap='viridis', s=20,
                                           vmin=scatter_global_vmin, vmax=scatter_global_vmax)
                        # Add colorbar horizontally at the bottom for the first plot only
                        if row_idx == 0 and col_idx == 0:
                            cbar = plt.colorbar(scatter, ax=axes, orientation='horizontal', 
                                              shrink=0.8, aspect=30, pad=0.25)
                            cbar.set_label('Chromatin Accessibility', fontsize=20)
                            cbar.ax.tick_params(labelsize=16)
                    else:
                        # If peak is not available, just plot points
                        ax.scatter(sub_umap[:, 0], sub_umap[:, 1], s=20, alpha=0.6)
                    
                    ax.set_title(f"{gt_label} | {age_group}", fontsize=20)
                    ax.set_xlabel("UMAP1", fontsize=20)
                    ax.set_ylabel("UMAP2", fontsize=20)

            fig.suptitle(
                f"NK Cells UMAP by Genotype and Epigenetic Age at {snp}\n"
                f"Colored by Accessibility at {peak}",
                fontsize=20
            )
            plot_path = os.path.join(fig_dir, f"NK_UMAP_{snp}_genotype_age_peak_{peak.replace(':', '_')}.png")
            # Use subplots_adjust instead of tight_layout to avoid colorbar warnings
            plt.subplots_adjust(left=0.05, right=0.95, top=0.85, bottom=0.35, wspace=0.1, hspace=0.35)
            plt.savefig(plot_path, bbox_inches="tight", dpi=300)
            plt.close()
            print(f"Successfully saved 6-panel genotype/age UMAP to: {plot_path}")

            # Create hexagonal binned version of the 6-panel plot
            fig_hex, axes_hex = plt.subplots(2, 3, figsize=(18, 14), sharex=True, sharey=True)

            # Calculate global min/max values for consistent color scaling across all hexbins
            if peak in plot_obs.columns:
                global_vmin = plot_obs[peak].min()
                global_vmax = plot_obs[peak].max()
            else:
                global_vmin, global_vmax = None, None

            # Store hexbin objects for colorbar
            hexbin_objects = []

            for row_idx, age_group in enumerate(['young', 'old']):
                for col_idx, gt in enumerate(genotype_order):
                    ax_hex = axes_hex[row_idx, col_idx]
                    # Use DataFrame indices instead of copying AnnData objects
                    cells_mask = (
                        (plot_obs[snp] == gt) &
                        (plot_obs['age_group'] == age_group)
                    )
                    gt_label = genotype_numeric_to_str.get(gt, gt)
                    if cells_mask.sum() == 0:
                        ax_hex.set_title(f"{gt_label} | {age_group}\nNo cells")
                        ax_hex.axis('off')
                        continue
                        
                    # Get UMAP coordinates for this subset
                    sub_umap = umap_coords[cells_mask]
                    
                    # Create hexbin plot
                    if peak in plot_obs.columns:
                        peak_values = plot_obs.loc[cells_mask, peak].values
                        # Remove any NaN values
                        valid_mask = ~np.isnan(peak_values)
                        if valid_mask.sum() > 0:  # Check if there are valid values
                            umap_coords_clean = sub_umap[valid_mask]
                            peak_values_clean = peak_values[valid_mask]
                            
                            # Create hexbin plot with consistent color scaling
                            hb = ax_hex.hexbin(umap_coords_clean[:, 0], umap_coords_clean[:, 1], 
                                              C=peak_values_clean, 
                                              gridsize=20, 
                                              cmap='viridis', 
                                              reduce_C_function=np.mean,
                                              vmin=global_vmin, 
                                              vmax=global_vmax)
                            hexbin_objects.append(hb)
                        else:
                            ax_hex.text(0.5, 0.5, 'No valid data', 
                                       horizontalalignment='center',
                                       verticalalignment='center',
                                       transform=ax_hex.transAxes)
                    else:
                        # If peak is not in obs, just plot density
                        hb = ax_hex.hexbin(sub_umap[:, 0], sub_umap[:, 1], 
                                     gridsize=20, 
                                     cmap='Blues')
                        hexbin_objects.append(hb)
                    
                    ax_hex.set_title(f"{gt_label} | {age_group}")
                    ax_hex.set_xlabel("UMAP1")
                    ax_hex.set_ylabel("UMAP2")

            # Add horizontal colorbar at bottom for hexbin plots
            if hexbin_objects and peak in plot_obs.columns:
                cbar_hex = fig_hex.colorbar(hexbin_objects[0], ax=axes_hex, orientation='horizontal', 
                                          shrink=0.8, aspect=30, pad=0.25)
                cbar_hex.set_label('Chromatin Accessibility', fontsize=20)
                cbar_hex.ax.tick_params(labelsize=16)

            fig_hex.suptitle(
                f"NK Cells UMAP by Genotype and Epigenetic Age at {snp} (Hexagonal Binned)\n"
                f"Colored by Accessibility at {peak}",
                fontsize=16
            )
            plot_hex_path = os.path.join(fig_dir, f"NK_UMAP_{snp}_genotype_age_peak_{peak.replace(':', '_')}_hexbin.png")
            # Use subplots_adjust instead of tight_layout to avoid colorbar warnings
            plt.subplots_adjust(left=0.05, right=0.95, top=0.85, bottom=0.35, wspace=0.1, hspace=0.35)
            plt.savefig(plot_hex_path, bbox_inches="tight", dpi=300)
            plt.close()
            print(f"Successfully saved 6-panel hexagonal binned genotype/age UMAP to: {plot_hex_path}")

        except Exception as e:
            print(f"An error occurred while processing SNP {snp}: {e}")

    print("\nPlotting finished.")


def load_nk_data(file_path="output/merged_standardFull_NK_cells.h5ad"):
    """
    Convenience function to load NK cell data.
    
    Parameters:
    -----------
    file_path : str
        Path to the NK cell h5ad file
        
    Returns:
    --------
    adata_merged : AnnData
        Loaded AnnData object
    """
    
    print(f"Loading NK cell data from: {file_path}")
    adata_merged = sc.read_h5ad(file_path)
    print(f"Loaded data: {adata_merged.shape} (cells x features)")
    print(f"Available metadata columns: {adata_merged.obs.columns.tolist()}")
    
    return adata_merged
