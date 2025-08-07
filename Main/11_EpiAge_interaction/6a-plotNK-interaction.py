# Import the functions
import importlib
import src.nk_umap_plotting
importlib.reload(src.nk_umap_plotting)
from src.nk_umap_plotting import plot_nk_genotype_interactions, plot_nk_epigenetic_age, load_nk_data


# Load your data (if you haven't already)
adata_merged = load_nk_data("output/merged_standardFull_NK_cells.h5ad")

# Plot epigenetic age UMAPs (both regular and hexbin versions)
plot_nk_epigenetic_age(adata_merged)

# Run the genotype interaction plotting function with default parameters
plot_nk_genotype_interactions(adata_merged)

# Or with custom parameters
plot_nk_genotype_interactions(
    adata_merged, 
    interaction_results_file="output/4-regression-results/combined_results_NK.csv",
    fig_dir="figures/5-epigeneticAge-Genotype-UMAPs",
    target_peaks=['chr4:55993465-55994658', 'chr6:31243369-31243947', 'chr7:5897477-5897880'],
    adj_pval_threshold=0.05,
    estimate_threshold=1,
    max_peaks=50
)