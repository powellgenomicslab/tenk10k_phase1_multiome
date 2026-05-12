#!/usr/bin/env python3
"""
Functions for Genotype × Epigenetic Age Interaction Analysis

This module contains all the functions used in 5-interactionPlotting.py
for processing, analyzing, and visualizing interaction results.

Author: Peter C Allen
"""

import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binomtest
import re


def load_and_combine_results(input_dir="output/4-regression-results"):
    """
    Load and combine interaction results from all chromosomes and cell types.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing interaction results files
        
    Returns:
    --------
    tuple
        (combined_df, significant_df, filtered_df) - raw results, significant results, and filtered results
    """
    
    print("Processing genotype × epigenetic age interaction results...")
    
    # Find all interaction results files
    file_paths = glob.glob(os.path.join(input_dir, "interaction_results_*_chr[0-9]*.csv"))
    cell_types = list({re.sub(r".*_results_(.+)_chr\d+\.csv", r"\1", os.path.basename(fp)) for fp in file_paths})

    print(f"Found results for {len(cell_types)} cell types")

    # Combine results by cell type
    for cell_type in cell_types:
        print(f"Processing {cell_type}...")
        
        # Find all chromosome files for this cell type
        cell_type_files = [fp for fp in file_paths if re.search(rf"interaction_results_{cell_type}_chr\d+\.csv", os.path.basename(fp))]
        
        if not cell_type_files:
            continue
        
        # Combine results across chromosomes
        combined_results = pd.concat([pd.read_csv(f) for f in cell_type_files], ignore_index=True)
        combined_results = combined_results.sort_values('adj_pval')
        
        # Save combined results
        output_path = f"{input_dir}/combined_results_{cell_type}.csv"
        combined_results.to_csv(output_path, index=False)

    # Combine all cell types for overall analysis
    print("Combining results across all cell types...")
    combined_files = glob.glob(f'{input_dir}/combined_results_*.csv')
    all_results = []

    for filepath in combined_files:
        cell_type = os.path.basename(filepath).replace('combined_results_', '').replace('.csv', '')
        df = pd.read_csv(filepath)
        df['cell_type'] = cell_type
        all_results.append(df)

    # Create master dataframe
    combined_df = pd.concat(all_results, ignore_index=True)
    print(f"Total tests performed: {len(combined_df):,}")

    # Filter for significant results
    significant_df = combined_df[combined_df['adj_pval'] < 0.05]
    print(f"Significant interactions (FDR < 0.05): {len(significant_df):,}")

    # Load color palette and merge
    palette_df = pd.read_csv("data/colour_palette_table.tsv", sep="\t")
    palette_df.columns = palette_df.columns.str.strip()
    palette_df = palette_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    # Map cell_type_plot from palette_df to filtered_df
    filtered_df = significant_df.merge(
        palette_df[["wg2_scpred_prediction", "cell_type"]],
        left_on="cell_type",
        right_on="wg2_scpred_prediction",
        how="left",
    )
    filtered_df.rename(columns={"cell_type_y": "cell_type_plot"}, inplace=True)
    
    return combined_df, significant_df, filtered_df


def annotate_snp_sharing(filtered_df):
    """
    Annotate SNPs as unique (cell type-specific) or shared across cell types.
    
    Parameters:
    -----------
    filtered_df : pd.DataFrame
        Filtered interaction results
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with snp_type column added (unique/shared)
    """
    
    # Calculate unique and shared SNPs per cell type
    snp_celltype = filtered_df.groupby("snp")["cell_type_plot"].nunique().reset_index()
    snp_celltype.columns = ["snp", "n_celltypes"]

    # Merge with filtered_df to annotate SNPs as unique/shared
    filtered_df = filtered_df.merge(snp_celltype, on="snp", how="left")
    filtered_df["snp_type"] = np.where(filtered_df["n_celltypes"] == 1, "unique", "shared")
    
    return filtered_df


def get_celltype_labels():
    """
    Get formatted cell type labels with subscripts for plotting.
    
    Returns:
    --------
    dict
        Dictionary mapping cell type names to formatted labels
    """
    
    return {
        "CD14 Mono": r"CD14$_{\mathrm{Mono}}$",
        "CD16 Mono": r"CD16$_{\mathrm{Mono}}$",
        "CD4 Naive": r"CD4$_{\mathrm{Naive}}$",
        "CD8 Naive": r"CD8$_{\mathrm{Naive}}$",
        "B naive": r"B$_{\mathrm{naive}}$",
        "B Naive": r"B$_{\mathrm{naive}}$",  # Handle both variants
        "B intermediate": r"B$_{\mathrm{intermediate}}$",
        "B memory": r"B$_{\mathrm{memory}}$",
        "CD4 TCM": r"CD4$_{\mathrm{TCM}}$",
        "CD8 TCM": r"CD8$_{\mathrm{TCM}}$",
        "CD8 Proliferating": r"CD8$_{\mathrm{Proliferating}}$",
        "NK CD56bright": r"NK$_{\mathrm{CD56bright}}$",
        "CD4 TEM": r"CD4$_{\mathrm{TEM}}$",
        "CD8 TEM": r"CD8$_{\mathrm{TEM}}$",
        "CD4 CTL": r"CD4$_{\mathrm{CTL}}$",
        "Treg": "Treg",
        "CD4 Proliferating": r"CD4$_{\mathrm{Proliferating}}$",
        "gdT": "gdT",
        "MAIT": "MAIT",
        "dnT": "dnT",
        "ILC": "ILC",
        "NK": "NK",
        "NK Proliferating": r"NK$_{\mathrm{Proliferating}}$",
        "Plasmablast": "Plasmablast",
        "cDC2": "cDC2",
        "pDC": "pDC",
        "cDC1": "cDC1",
        "ASDC": "ASDC",
        "HSPC": "HSPC"
    }


def get_celltype_order():
    """
    Get the desired order for cell types in plots.
    
    Returns:
    --------
    list
        List of cell type names in desired order
    """
    
    return [
        "CD4 TCM", "CD4 Naive", "CD4 TEM", "CD4 CTL", "Treg", "CD4 Proliferating",
        "gdT", "MAIT", "dnT", "ILC", "CD8 TEM", "CD8 Naive", "CD8 TCM", "CD8 Proliferating",
        "NK", "NK CD56bright", "NK Proliferating", "B naive", "B intermediate", "B memory",
        "Plasmablast", "CD14 Mono", "CD16 Mono", "cDC2", "pDC", "cDC1", "ASDC", "HSPC"
    ]


def plot_unique_shared_snps(filtered_df, output_path="figures/3-regression/numberSig-caSNPs-unique-shared-by-celltype.png"):
    """
    Create stacked bar plot showing unique vs shared SNPs by cell type.
    
    Parameters:
    -----------
    filtered_df : pd.DataFrame
        Filtered interaction results with snp_type annotation
    output_path : str
        Path to save the plot
        
    Returns:
    --------
    None
        Saves plot to output_path
    """
    
    print("\nGenerating Plot: Unique vs shared SNPs by cell type...")
    
    celltype_labels = get_celltype_labels()
    desired_order = get_celltype_order()

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

    # Calculate total SNPs per cell type
    snp_counts_ordered = pd.DataFrame({
        "cell_type_plot": desired_order,
        "snp": [snp_type_counts.loc[ct, ["unique", "shared"]].sum() if ct in snp_type_counts.index else 0 for ct in desired_order]
    })

    # Create the plot
    plt.figure(figsize=(16, 8))
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
            fontsize=18,
            fontweight="bold"
        )

    plt.xlabel("", fontsize=18)
    plt.ylabel("Number of Significant Interactions", fontsize=20)

    # Format x-tick labels with subscripts
    formatted_labels = [celltype_labels.get(ct, ct) for ct in snp_type_counts.index]
    plt.xticks(range(len(snp_type_counts.index)), formatted_labels, rotation=90, fontsize=20)
    plt.yticks(fontsize=20)

    plt.xlim(-0.5, len(snp_type_counts.index) - 0.5)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    plt.legend(title="SNP Type", fontsize=18, title_fontsize=20, loc="upper right")
    plt.tight_layout()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Plot saved: {output_path}")


def analyze_estimate_directions(filtered_df):
    """
    Perform statistical analysis of estimate directions (positive vs negative).
    
    Parameters:
    -----------
    filtered_df : pd.DataFrame
        Filtered interaction results
        
    Returns:
    --------
    dict
        Dictionary with analysis results
    """
    
    print("\nPerforming statistical analysis of estimate directions...")

    # Count positive and negative estimates
    positive_count = (filtered_df['Estimate'] > 0).sum()
    negative_count = (filtered_df['Estimate'] < 0).sum()
    total_count = positive_count + negative_count

    if total_count > 0:
        positive_percentage = (positive_count / total_count) * 100
        negative_percentage = (negative_count / total_count) * 100
    else:
        positive_percentage = 0
        negative_percentage = 0

    # Perform binomial test
    p_value = binomtest(positive_count, n=total_count, p=0.5, alternative='greater')
    
    results = {
        'positive_count': positive_count,
        'negative_count': negative_count,
        'total_count': total_count,
        'positive_percentage': positive_percentage,
        'negative_percentage': negative_percentage,
        'binomial_pvalue': p_value.pvalue
    }
    
    print(f"Positive estimates: {positive_count} ({positive_percentage:.2f}%)")
    print(f"Negative estimates: {negative_count} ({negative_percentage:.2f}%)")
    print(f"Binomial test p-value (positive > negative): {p_value.pvalue:.4g}")
    
    return results


def prepare_gwas_hits_data(filtered_df):
    """
    Prepare GWAS hits data for boxplot analysis.
    
    Parameters:
    -----------
    filtered_df : pd.DataFrame
        Filtered interaction results
        
    Returns:
    --------
    pd.DataFrame
        GWAS hits data ready for plotting
    """
    
    print("\nPreparing data for interaction boxplots...")

    # Add cell_type column for compatibility
    filtered_df['cell_type'] = filtered_df['cell_type_x']

    # Define GWAS hits table
    gwas_hits_table = pd.DataFrame({
        "cell_type": ["NK", "NK", "NK"],
        "Chr": [4, 6, 7],
        "Peak": ["chr4:55993465-55994658", "chr6:31243369-31243947", "chr7:5897477-5897880"]
    })

    # Import and use the boxplot filtering function
    import sys
    import os
    # Add parent directory to path to access interaction_boxplots
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, parent_dir)
    
    from interaction_boxplots import filter_gwas_hits
    gwas_hits = filter_gwas_hits(filtered_df, gwas_hits_table)
    print(f"Found {len(gwas_hits)} GWAS hits to plot")
    
    return gwas_hits


def create_gwas_boxplots(gwas_hits, output_dir="figures/3-regression/gwas_boxplots"):
    """
    Create boxplots for GWAS hits.
    
    Parameters:
    -----------
    gwas_hits : pd.DataFrame
        GWAS hits data
    output_dir : str
        Output directory for boxplots
        
    Returns:
    --------
    None
        Saves boxplots to output_dir
    """
    
    print("\nGenerating GWAS hit boxplots...")
    
    # Import the boxplot function
    import sys
    import os
    # Add parent directory to path to access interaction_boxplots
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, parent_dir)
    
    from interaction_boxplots import load_and_plot_interactions
    
    # Create boxplots for GWAS hits
    load_and_plot_interactions(
        gwas_hits, 
        output_dir=output_dir
    )


def create_strongest_interaction_boxplots(filtered_df, 
                                        maf_threshold=0.12, 
                                        n_positive=50, 
                                        n_negative=5,
                                        output_dir_positive="figures/3-regression/strongest_positive_boxplots",
                                        output_dir_negative="figures/3-regression/strongest_negative_boxplots"):
    """
    Create boxplots for strongest positive and negative interactions.
    
    Parameters:
    -----------
    filtered_df : pd.DataFrame
        Filtered interaction results
    maf_threshold : float
        Minimum minor allele frequency threshold
    n_positive : int
        Number of strongest positive interactions
    n_negative : int
        Number of strongest negative interactions
    output_dir_positive : str
        Output directory for positive interaction boxplots
    output_dir_negative : str
        Output directory for negative interaction boxplots
        
    Returns:
    --------
    None
        Saves boxplots to output directories
    """
    
    print("\nFiltering strongest interactions...")
    
    # Import the filtering and plotting functions
    import sys
    import os
    # Add parent directory to path to access interaction_boxplots
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, parent_dir)
    
    from interaction_boxplots import filter_strongest_interactions, load_and_plot_interactions
    
    strongest_positive, strongest_negative = filter_strongest_interactions(
        filtered_df, 
        maf_threshold=maf_threshold, 
        n_positive=n_positive, 
        n_negative=n_negative
    )

    print(f"\nPlotting {len(strongest_positive)} strongest positive interactions...")
    load_and_plot_interactions(
        strongest_positive, 
        output_dir=output_dir_positive
    )

    print(f"Plotting {len(strongest_negative)} strongest negative interactions...")
    load_and_plot_interactions(
        strongest_negative, 
        output_dir=output_dir_negative
    )


def run_full_analysis(input_dir="output/4-regression-results", 
                     output_dir="figures/",
                     create_gwas_plots=True,
                     create_strongest_plots=False):
    """
    Run the complete interaction analysis workflow.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing interaction results
    output_dir : str
        Base output directory for figures
    create_gwas_plots : bool
        Whether to create GWAS hit boxplots
    create_strongest_plots : bool
        Whether to create strongest interaction boxplots
        
    Returns:
    --------
    dict
        Dictionary containing all analysis results and dataframes
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load and combine data
    combined_df, significant_df, filtered_df = load_and_combine_results(input_dir)
    
    # Annotate SNP sharing
    filtered_df = annotate_snp_sharing(filtered_df)
    
    # Create unique/shared SNPs plot
    plot_unique_shared_snps(filtered_df)
    
    # Analyze estimate directions
    direction_analysis = analyze_estimate_directions(filtered_df)
    
    # Prepare GWAS hits
    gwas_hits = prepare_gwas_hits_data(filtered_df)
    
    # Create plots
    if create_gwas_plots:
        create_gwas_boxplots(gwas_hits)
    
    if create_strongest_plots:
        create_strongest_interaction_boxplots(filtered_df)
    
    # Return all results
    return {
        'combined_df': combined_df,
        'significant_df': significant_df,
        'filtered_df': filtered_df,
        'gwas_hits': gwas_hits,
        'direction_analysis': direction_analysis
    }
