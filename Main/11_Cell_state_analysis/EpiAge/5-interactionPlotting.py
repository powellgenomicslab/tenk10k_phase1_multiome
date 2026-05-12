#!/usr/bin/env python3
"""
Genotype × Epigenetic Age Interaction Results Visualization

This script processes and visualizes the results from genotype-epigenetic age 
interaction analyses across all cell types. It creates summary plots showing
the distribution of significant interactions and their characteristics.

All functions have been modularized into 5-functions.py for better code organization.

Author: Peter C Allen
"""

# Import the analysis functions from the same directory
from functions_5 import (
    load_and_combine_results,
    annotate_snp_sharing,
    plot_unique_shared_snps,
    analyze_estimate_directions,
    prepare_gwas_hits_data,
    create_gwas_boxplots,
    create_strongest_interaction_boxplots,
    run_full_analysis
)

# ============================================================================
# MAIN ANALYSIS WORKFLOW
# ============================================================================

def main():
    """
    Main function to run the complete interaction analysis.
    """
    
    # Configuration
    input_dir = "output/4-regression-results"
    output_dir = "figures/"
    
    print("Starting genotype × epigenetic age interaction analysis...")
    
    # Option 1: Run the complete analysis in one function call
    print("\n" + "="*60)
    print("RUNNING COMPLETE ANALYSIS")
    print("="*60)
    
    results = run_full_analysis(
        input_dir=input_dir,
        output_dir=output_dir,
        create_gwas_plots=True,
        create_strongest_plots=False  # Set to True if you want strongest interaction plots
    )
    
    print(f"\nAnalysis complete! Results summary:")
    print(f"- Total tests: {len(results['combined_df']):,}")
    print(f"- Significant interactions: {len(results['significant_df']):,}")
    print(f"- GWAS hits for plotting: {len(results['gwas_hits'])}")
    print(f"- Positive estimates: {results['direction_analysis']['positive_count']} ({results['direction_analysis']['positive_percentage']:.2f}%)")
    print(f"- Negative estimates: {results['direction_analysis']['negative_count']} ({results['direction_analysis']['negative_percentage']:.2f}%)")


def main_step_by_step():
    """
    Alternative main function that runs analysis step by step for more control.
    Uncomment this and comment out main() if you want to run steps individually.
    """
    
    # Configuration
    input_dir = "output/4-regression-results"
    output_dir = "figures/"
    
    print("Starting step-by-step genotype × epigenetic age interaction analysis...")
    
    # Step 1: Load and combine data
    print("\n" + "="*60)
    print("STEP 1: DATA LOADING AND PREPROCESSING")
    print("="*60)
    
    combined_df, significant_df, filtered_df = load_and_combine_results(input_dir)
    
    # Step 2: Annotate SNP sharing
    print("\n" + "="*60)
    print("STEP 2: SNP SHARING ANNOTATION")
    print("="*60)
    
    filtered_df = annotate_snp_sharing(filtered_df)
    
    # Step 3: Create unique/shared SNPs plot
    print("\n" + "="*60)
    print("STEP 3: UNIQUE VS SHARED SNPS PLOT")
    print("="*60)
    
    plot_unique_shared_snps(filtered_df)
    
    # Step 4: Analyze estimate directions
    print("\n" + "="*60)
    print("STEP 4: STATISTICAL ANALYSIS")
    print("="*60)
    
    direction_analysis = analyze_estimate_directions(filtered_df)
    
    # Step 5: Prepare and create GWAS boxplots
    print("\n" + "="*60)
    print("STEP 5: GWAS HIT BOXPLOTS")
    print("="*60)
    
    gwas_hits = prepare_gwas_hits_data(filtered_df)
    create_gwas_boxplots(gwas_hits)
    
    # Step 6: Optional - Create strongest interaction boxplots
    # Uncomment the following lines if you want these plots:
    
    # print("\n" + "="*60)
    # print("STEP 6: STRONGEST INTERACTION BOXPLOTS")
    # print("="*60)
    # 
    # create_strongest_interaction_boxplots(
    #     filtered_df,
    #     maf_threshold=0.12,
    #     n_positive=50,
    #     n_negative=5
    # )
    
    print("\nStep-by-step analysis complete!")


# ============================================================================
# CUSTOM ANALYSIS EXAMPLES
# ============================================================================

def custom_analysis_example():
    """
    Example of how to run custom analysis with specific parameters.
    """
    
    print("Running custom analysis example...")
    
    # Load data
    combined_df, significant_df, filtered_df = load_and_combine_results()
    filtered_df = annotate_snp_sharing(filtered_df)
    
    # Custom plot with different output path
    plot_unique_shared_snps(
        filtered_df, 
        output_path="figures/custom/my_custom_snp_plot.png"
    )
    
    # Create boxplots with custom parameters
    gwas_hits = prepare_gwas_hits_data(filtered_df)
    create_gwas_boxplots(
        gwas_hits, 
        output_dir="figures/custom/my_gwas_boxplots"
    )
    
    # Custom strongest interactions analysis
    create_strongest_interaction_boxplots(
        filtered_df,
        maf_threshold=0.15,  # Higher MAF threshold
        n_positive=20,       # Fewer positive interactions
        n_negative=10,       # More negative interactions
        output_dir_positive="figures/custom/top20_positive",
        output_dir_negative="figures/custom/top10_negative"
    )


# ============================================================================
# SCRIPT EXECUTION
# ============================================================================

if __name__ == "__main__":
    # Run the main analysis
    main()
    
    # Alternatively, run step-by-step analysis (uncomment the line below)
    # main_step_by_step()
    
    # Or run custom analysis (uncomment the line below)
    # custom_analysis_example()
    
    print("\nScript completed successfully!")

# ============================================================================
# USAGE NOTES
# ============================================================================
"""
This script has been completely modularized. All functions are now in 5-functions.py.

Main functions available:
- run_full_analysis(): Complete workflow in one call
- load_and_combine_results(): Load and process data
- plot_unique_shared_snps(): Create stacked bar plot
- analyze_estimate_directions(): Statistical analysis
- create_gwas_boxplots(): GWAS hit boxplots
- create_strongest_interaction_boxplots(): Strongest interaction boxplots

You can also import individual functions:
from src.2_subset_fullTenk10k.5_functions import plot_unique_shared_snps, analyze_estimate_directions

For boxplot functions, they're still available from:
from src.interaction_boxplots import load_and_plot_interactions, filter_gwas_hits
"""
