#!/usr/bin/env python3
"""
Interaction Boxplot Functions

This module contains functions for creating boxplots of SNP x EpiAge interactions
on chromatin accessibility, including statistical testing and visualization.

Usage:
    from src.interaction_boxplots import load_and_plot_interactions
    
    # Create boxplots for variants
    load_and_plot_interactions(
        variants_df, 
        output_dir="figures/3-regression/caSNPxPeak-boxplots"
    )

Author: Peter C Allen
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import gc
from scipy import stats


def load_and_plot_interactions(
    variants_df,
    output_dir="figures/3-regression/caSNPxPeak-boxplots",
    expr_threshold=10,
    figure_size=(12,6), 
    palette={"Young": "#7FB2FF", "Old": "#D6C7A1"}
):
    """
    Create boxplots for SNP x EpiAge interactions on chromatin accessibility, with t-test significance.
    
    Parameters:
    -----------
    variants_df : pd.DataFrame
        DataFrame containing interaction results with columns: cell_type_x, peak, snp, Estimate
    output_dir : str
        Output directory for boxplot figures
    expr_threshold : float
        Maximum chromatin accessibility value to include in plots
    figure_size : tuple
        Figure size as (width, height)
    palette : dict
        Color palette for Young/Old age groups
        
    Returns:
    --------
    None
        Saves boxplot figures to output_dir
    """
    
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    for i, (idx, row) in enumerate(variants_df.iterrows()):
        print(f"\nProcessing variant {i+1}/{len(variants_df)}")
        cell_type = row["cell_type_x"]
        peak = row["peak"]
        snp = row["snp"]
        estimate = row["Estimate"]
        print(f"Cell type: {cell_type}, Peak: {peak}, SNP: {snp}")

        match = re.match(r"^chr(\d+)", peak)
        if match:
            chr_num = match.group(1)
        else:
            print(f"Could not extract chromosome number from peak: {peak}")
            continue

        try:
            # Load expression data
            expr_file = (
                f"data/expressionBeds/{cell_type}/ExpressionBeds/chr{chr_num}.bed.gz"
            )
            print(f"Loading expression data from {expr_file}")
            expr_df = pd.read_csv(expr_file, sep="\t", compression="gzip", index_col=3)

            # Load covariates and epiAge
            cov_file = f"data/expressionBeds/{cell_type}/covariates.txt"
            print(f"Loading covariates from {cov_file}")
            cov_df = pd.read_csv(cov_file, index_col=0).T
            cov_df.index = cov_df.index.str.replace("X", "")

            epiage_file = f"data/20250630_epiages/{cell_type}_epiage.txt"
            print(f"Loading epiAge data from {epiage_file}")
            epiage_df = pd.read_csv(epiage_file, sep="\t").dropna(
                subset=["median_epiage"]
            )
            epiage_df["epiAge_binary"] = np.where(
                epiage_df["median_epiage"] > 0.5, "Old", "Young"
            )

            # Load genotype data efficiently - only read header first to find SNP column
            geno_file = f"data/genotype/wgs/TenK10K_TOB_ATAC_renamed_chr{chr_num}_common_variants.raw"
            print(f"Loading genotype data from {geno_file}")
            
            # Read header to find which column contains our SNP
            with open(geno_file, 'r') as f:
                header = f.readline().strip().split()
            
            # Find SNP column (may have _A, _G, _T, _C suffix)
            snp_cols = [i for i, col in enumerate(header) if col.startswith(snp)]
            if not snp_cols:
                print(f"SNP {snp} not found in genotype file headers")
                continue
            
            # Load only IID column and the SNP column
            iid_col = header.index('IID')
            snp_col = snp_cols[0]  # Take first match
            usecols = [iid_col, snp_col]
            
            geno_df = pd.read_csv(geno_file, sep=r"\s+", usecols=usecols, index_col='IID')
            geno_df.columns = [header[snp_col].split('_')[0]]  # Clean column name

            if peak not in expr_df.index:
                print(f"Peak {peak} not found in expression data for {cell_type}.")
                continue
            
            # Use the cleaned column name (should be the SNP name)
            snp_col_name = geno_df.columns[0]
            if snp_col_name != snp:
                print(f"Expected SNP {snp} but got column {snp_col_name}")

            print("Preparing data for plotting...")
            peak_expr = expr_df.loc[peak].iloc[
                3:
            ]  # Skip first 3 columns (chr, start, end)
            snp_geno = geno_df.iloc[:, 0]  # Get the first (and only) column

            # Ensure indices match between peak_expr and snp_geno
            common_ids = peak_expr.index.intersection(snp_geno.index)
            if len(common_ids) == 0:
                print(
                    f"No overlapping individuals between expression and genotype for {cell_type}, {peak}, {snp}."
                )
                continue

            peak_expr = peak_expr.loc[common_ids]
            snp_geno = snp_geno.loc[common_ids]

            # Merge data
            plot_data = pd.DataFrame(
                {
                    "ID": peak_expr.index,
                    "expr": peak_expr.values,
                    "genotype": snp_geno.astype(str),
                }
            )

            if "individual" not in epiage_df.columns:
                print(f"'individual' column missing in epiage file for {cell_type}.")
                continue

            print("Merging with epiAge data...")
            plot_data = plot_data.merge(
                epiage_df[["individual", "epiAge_binary"]],
                left_on="ID",
                right_on="individual",
                how="left",
            )
            plot_data = plot_data.dropna()

            # Remove rows where genotype is NA before plotting
            plot_data = plot_data[plot_data["genotype"].notna()]
            plot_data = plot_data[plot_data["genotype"] != "nan"]

            # Filter out chromatin accessibility values above threshold
            plot_data = plot_data[plot_data["expr"] <= expr_threshold]

            # Order epiAge_binary: Young first, then Old
            plot_data["epiAge_binary"] = pd.Categorical(
                plot_data["epiAge_binary"], categories=["Young", "Old"], ordered=True
            )

            # Create the boxplot
            _create_interaction_boxplot(
                plot_data, snp, peak, cell_type, estimate, 
                output_dir, i+1, figure_size, palette
            )

            # Explicitly delete large objects and collect garbage to free memory
            del expr_df, cov_df, epiage_df, geno_df, peak_expr, snp_geno, plot_data
            gc.collect()
            print("Memory cleaned up for this variant.")

        except Exception as e:
            import traceback
            print(f"Error plotting {cell_type} {snp}: {e}")
            traceback.print_exc()
            continue


def _create_interaction_boxplot(
    plot_data, snp, peak, cell_type, estimate, 
    output_dir, variant_num, figure_size, palette
):
    """
    Internal function to create and save a single interaction boxplot.
    
    Parameters:
    -----------
    plot_data : pd.DataFrame
        Prepared data for plotting with columns: expr, genotype, epiAge_binary
    snp : str
        SNP identifier
    peak : str
        Peak identifier
    cell_type : str
        Cell type name
    estimate : float
        Interaction estimate value
    output_dir : str
        Output directory for figure
    variant_num : int
        Variant number for filename
    figure_size : tuple
        Figure size as (width, height)
    palette : dict
        Color palette for age groups
    """
    
    print("Plotting boxplot...")
    plt.figure(figsize=figure_size)
    # Set classic theme: white background, no grid, black axes
    plt.style.use("default")
    ax = plt.gca()
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.grid(False)

    # Create custom genotype labels based on SNP
    def get_genotype_labels(snp_id):
        """Get custom genotype labels based on SNP alleles"""
        if "6:31245967:T:A" in snp_id:
            return {'0': 'TT', '1': 'TA', '2': 'AA'}
        else:
            # Parse alleles from SNP name for other SNPs
            snp_parts = snp_id.split(":")
            if len(snp_parts) >= 4:
                ref_allele = snp_parts[2]
                alt_allele = snp_parts[3]
                return {
                    '0': ref_allele * 2,
                    '1': ref_allele + alt_allele,
                    '2': alt_allele * 2
                }
            else:
                return {'0': '0', '1': '1', '2': '2'}
    
    genotype_labels = get_genotype_labels(snp)

        # Draw violin plot by age group (Young/Old) with genotype subgroups
    ax = sns.violinplot(
        data=plot_data,
        x="genotype",
        y="expr",
        hue="epiAge_binary",
        palette=palette,
        order=sorted(plot_data["genotype"].unique()),
        hue_order=["Young", "Old"],
        dodge=True,
        inner=None,  # No inner markings
        linewidth=0,
        saturation=0.7
    )
    
    # Get the positions of the violin plots to place boxplots correctly
    violin_positions = []
    for i, genotype in enumerate(sorted(plot_data["genotype"].unique())):
        # For dodged plots, Young is left (-0.2) and Old is right (+0.2) relative to center
        violin_positions.extend([i - 0.2, i + 0.2])
    
    # Create separate datasets for Young and Old to overlay boxplots at correct positions
    young_data = plot_data[plot_data["epiAge_binary"] == "Young"]
    old_data = plot_data[plot_data["epiAge_binary"] == "Old"]
    
    # Overlay boxplots for Young age group
    if not young_data.empty:
        for i, genotype in enumerate(sorted(plot_data["genotype"].unique())):
            genotype_data = young_data[young_data["genotype"] == genotype]
            if not genotype_data.empty:
                bp = ax.boxplot(
                    genotype_data["expr"].values,
                    positions=[i - 0.2],  # Young position
                    widths=0.1,
                    patch_artist=True,
                    boxprops=dict(facecolor="white", edgecolor="black", linewidth=1.5, alpha=0.3),
                    medianprops=dict(color="black", linewidth=2),
                    whiskerprops=dict(color="black", linewidth=1.5),
                    capprops=dict(color="black", linewidth=1.5),
                    flierprops=dict(marker='o', markersize=3, markerfacecolor='black', markeredgecolor='black')
                )
    
    # Overlay boxplots for Old age group  
    if not old_data.empty:
        for i, genotype in enumerate(sorted(plot_data["genotype"].unique())):
            genotype_data = old_data[old_data["genotype"] == genotype]
            if not genotype_data.empty:
                bp = ax.boxplot(
                    genotype_data["expr"].values,
                    positions=[i + 0.2],  # Old position
                    widths=0.1,
                    patch_artist=True,
                    boxprops=dict(facecolor="white", edgecolor="black", linewidth=1.5, alpha=0.3),
                    medianprops=dict(color="black", linewidth=2),
                    whiskerprops=dict(color="black", linewidth=1.5),
                    capprops=dict(color="black", linewidth=1.5),
                    flierprops=dict(marker='o', markersize=3, markerfacecolor='black', markeredgecolor='black')
                )

    # Add regression lines for Young and Old age groups
    _add_regression_lines(ax, plot_data, palette)

    # Add sample counts only (no t-test significance)
    _add_sample_counts_only(ax, plot_data, genotype_labels)

    # Axis labels and title with larger fonts
    plt.xlabel("Genotype", fontsize=24)
    plt.ylabel("Chromatin Accessibility", fontsize=24)
    
    # Add title with SNP: chr#:bp and extra spacing
    snp_parts = snp.split(":")
    if len(snp_parts) >= 2:
        snp_chr_bp = f"chr{snp_parts[0]}:{snp_parts[1]}"
    else:
        snp_chr_bp = snp  # fallback
    plt.title(f"SNP: {snp_chr_bp}", fontsize=26, pad=20)  # Added pad=20 for spacing

    # Increase font sizes for ticks
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    
    # Update x-tick labels with custom genotype labels
    # The violin plot with hue creates multiple ticks per genotype
    # We need to get the current tick positions and labels
    current_ticks = ax.get_xticks()
    unique_genotypes = sorted(plot_data["genotype"].unique())
    
    # For dodged violins, we have 2 ticks per genotype (Young, Old)
    # We want to place the genotype label at the center of each pair
    if len(current_ticks) == len(unique_genotypes) * 2:
        # Set tick positions to be at the center of each genotype group
        new_tick_positions = [i for i in range(len(unique_genotypes))]
        custom_labels = [genotype_labels.get(gt, gt) for gt in unique_genotypes]
        ax.set_xticks(new_tick_positions)
        ax.set_xticklabels(custom_labels)
    else:
        # Fallback: use original approach if tick structure is unexpected
        custom_labels = [genotype_labels.get(gt, gt) for gt in unique_genotypes]
        ax.set_xticks(range(len(unique_genotypes)))
        ax.set_xticklabels(custom_labels)

    # Set y-axis maximum to accommodate sample count labels
    ymax = plot_data["expr"].max() + 0.1
    ax.set_ylim(bottom=ax.get_ylim()[0], top=ymax)

    # Remove top and right spines for classic look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Place legend outside the plot with larger fonts
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(
        handles,
        labels,
        title="EpiAge Group",
        title_fontsize=22,
        fontsize=20,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0,
        frameon=False,
    )

    # Add slope annotations below the legend
    _add_slope_annotations(ax, legend)

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space for legend

    # Add pos_ or neg_ prefix to filename based on Estimate
    prefix = "pos_" if estimate > 0 else "neg_"
    filename = f"{variant_num}_{cell_type}_{peak}_caSNPxPeak_boxplot.png"
    plot_path = os.path.join(output_dir, filename)
    print(f"Saving plot to {plot_path}")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()


def _add_regression_lines(ax, plot_data, palette):
    """
    Add regression lines for Young and Old age groups showing trend of chromatin 
    accessibility across genotypes.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object to add regression lines to
    plot_data : pd.DataFrame
        Plot data with columns: genotype, epiAge_binary, expr
    palette : dict
        Color palette for age groups (Young/Old)
    """
    
    # Convert genotype to numeric for regression
    plot_data_numeric = plot_data.copy()
    try:
        plot_data_numeric['genotype_numeric'] = pd.to_numeric(plot_data_numeric['genotype'])
    except (ValueError, TypeError):
        print("Warning: Could not convert genotype to numeric values for regression")
        return
    
    # Get unique genotypes for positioning
    unique_genotypes = sorted(plot_data["genotype"].unique())
    
    # Add regression lines for each age group
    for age_group in ["Young", "Old"]:
        age_data = plot_data_numeric[plot_data_numeric["epiAge_binary"] == age_group].copy()
        
        # Check if we have sufficient data for regression
        if len(age_data) < 2:
            print(f"Warning: Insufficient data points ({len(age_data)}) for {age_group} regression line")
            continue
            
        # Check if we have variation in genotype values
        unique_genotypes_in_age = age_data['genotype_numeric'].nunique()
        if unique_genotypes_in_age < 2:
            print(f"Warning: No genotype variation for {age_group} group - cannot draw regression line")
            continue
        
        # Remove any NaN values
        age_data = age_data.dropna(subset=['genotype_numeric', 'expr'])
        
        if len(age_data) < 2:
            print(f"Warning: Insufficient valid data points for {age_group} after removing NaN values")
            continue
        
        try:
            # Get x and y values as numpy arrays
            x_vals = age_data['genotype_numeric'].values.astype(float)
            y_vals = age_data['expr'].values.astype(float)
            
            # Ensure we have valid data
            if len(x_vals) != len(y_vals) or len(x_vals) < 2:
                print(f"Warning: Invalid data dimensions for {age_group}")
                continue
                
            # Perform linear regression using numpy polyfit as backup
            try:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
            except:
                # Fallback to numpy polyfit
                coeffs = np.polyfit(x_vals, y_vals, 1)
                slope, intercept = coeffs[0], coeffs[1]
            
            # Calculate regression line points positioned at violin plot locations
            x_positions = []
            y_positions = []
            
            for i, genotype in enumerate(unique_genotypes):
                genotype_numeric = float(genotype)
                y_pred = slope * genotype_numeric + intercept
                
                # Position based on age group (same as violin plot positions)
                if age_group == "Young":
                    x_pos = i - 0.2
                else:  # Old
                    x_pos = i + 0.2
                    
                x_positions.append(x_pos)
                y_positions.append(y_pred)
            
            # Draw regression line
            ax.plot(
                x_positions, 
                y_positions, 
                color=palette[age_group], 
                linewidth=2.5, 
                linestyle='-',
                alpha=0.8,
                label=f'{age_group} trend'
            )
            
            # Store slope for later annotation (will be added after legend)
            if not hasattr(ax, '_slope_annotations'):
                ax._slope_annotations = []
            ax._slope_annotations.append((age_group, slope, palette[age_group]))
            
        except Exception as e:
            print(f"Warning: Could not calculate regression line for {age_group}: {e}")
            continue


def _add_slope_annotations(ax, legend):
    """
    Add slope annotations below the legend.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object containing the plot
    legend : matplotlib.legend.Legend
        The legend object to position slopes below
    """
    
    if not hasattr(ax, '_slope_annotations') or not ax._slope_annotations:
        return
    
    # Get legend position
    legend_bbox = legend.get_window_extent()
    fig = ax.get_figure()
    
    # Convert legend position to axes coordinates
    legend_bbox_axes = legend_bbox.transformed(ax.transAxes.inverted())
    
    # Position slopes below the legend
    y_start = legend_bbox_axes.y0 - 0.05  # Start slightly below legend
    
    for i, (age_group, slope, color) in enumerate(ax._slope_annotations):
        y_pos = y_start - (i * 0.08)  # Space out the slope annotations
        
        ax.text(
            legend_bbox_axes.x0,  # Align with left edge of legend
            y_pos,
            f'{age_group} slope: {slope:.3f}',
            color=color,
            fontsize=16,
            fontweight='bold',
            ha='left',
            va='top',
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.3", facecolor='white', edgecolor=color, alpha=0.8)
        )


def _add_sample_counts_only(ax, plot_data, genotype_labels):
    """
    Add only sample counts to boxplot (no statistical tests).
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object to annotate
    plot_data : pd.DataFrame
        Plot data with columns: genotype, epiAge_binary, expr
    genotype_labels : dict
        Mapping of numeric genotypes to string labels
    """
    
    unique_genotypes = sorted(plot_data["genotype"].unique())
    
    for i_tick, genotype in enumerate(unique_genotypes):
        for j, age_group in enumerate(["Young", "Old"]):
            count = len(
                plot_data[
                    (plot_data["genotype"] == genotype)
                    & (plot_data["epiAge_binary"] == age_group)
                ]
            )
            if count > 0:
                # Place count above each box
                x_pos = i_tick - 0.2 + 0.4 * j  # Young left, Old right
                ax.text(
                    x_pos,
                    plot_data["expr"].max() + 0.1,
                    str(count),
                    ha="center",
                    va="bottom",
                    fontsize=16,  # Increased font size
                    fontweight="normal"  # Changed from "bold" to "normal"
                )


def filter_gwas_hits(filtered_df, gwas_hits_table):
    """
    Filter interaction results to match GWAS hits.
    
    Parameters:
    -----------
    filtered_df : pd.DataFrame
        Filtered interaction results with significant hits
    gwas_hits_table : pd.DataFrame
        Table with GWAS hits containing columns: cell_type, Chr, Peak
        
    Returns:
    --------
    pd.DataFrame
        Filtered dataframe containing only GWAS hits
    """
    
    # Standardize cell type column names to match filtered_df
    gwas_hits_table = gwas_hits_table.copy()
    gwas_hits_table["cell_type"] = gwas_hits_table["cell_type"].str.replace("_", " ")

    # Filter filtered_df to match cell_type and peak in gwas_hits_table
    gwas_hits = filtered_df[
        filtered_df.apply(
            lambda row: ((row["cell_type"].replace("_", " ") in gwas_hits_table["cell_type"].values) and
                         (row["peak"] in gwas_hits_table[gwas_hits_table["cell_type"] == row["cell_type"].replace("_", " ")]["Peak"].values)),
            axis=1
        )
    ]
    
    return gwas_hits


def filter_strongest_interactions(filtered_df, maf_threshold=0.12, n_positive=50, n_negative=5):
    """
    Filter for strongest positive and negative interaction estimates.
    
    Parameters:
    -----------
    filtered_df : pd.DataFrame
        Filtered interaction results
    maf_threshold : float
        Minimum minor allele frequency threshold
    n_positive : int
        Number of strongest positive interactions to return
    n_negative : int
        Number of strongest negative interactions to return
        
    Returns:
    --------
    tuple
        (strongest_positive, strongest_negative) DataFrames
    """
    
    # Filter by MAF and get strongest interactions
    high_maf_df = filtered_df[filtered_df['MAF'] > maf_threshold]
    
    strongest_positive = high_maf_df[high_maf_df['Estimate'] > 0].nlargest(n_positive, 'Estimate')
    strongest_negative = high_maf_df[high_maf_df['Estimate'] < 0].nsmallest(n_negative, 'Estimate')
    
    return strongest_positive, strongest_negative
