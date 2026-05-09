import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import sparse
from typing import Optional, Tuple, Dict
import logging
import glob
import gc

pool_files = glob.glob("/g/data/ei56/od8037/TenK10K/CreateObjects/Objects/*.h5ad") 
celltype_col = "predicted.id"

global_sums = {}
global_counts = {}
peak_names = None

print(f"Found {len(pool_files)} pool files to process.")

for i, file in enumerate(pool_files):
    print(f"Processing pool {i+1}/{len(pool_files)}: {file}")
    
    # Load just this pool into memory
    adata = sc.read_h5ad(file)
    
    # Capture peak names from the first file to use later
    if peak_names is None:
        peak_names = adata.var_names.copy()
    
    # Tally up sums and cell counts for each cell type in this pool
    for ct in adata.obs[celltype_col].unique():
        mask = adata.obs[celltype_col] == ct
        
        # Subset the AnnData object first, then grab X
        X_ct = adata[mask].X
        
        n_cells_in_pool = X_ct.shape[0]
        # Sum across cells (axis=0). Converts sparse matrix to dense 1D array
        ct_sum = np.asarray(X_ct.sum(axis=0)).flatten()
        
        if ct not in global_sums:
            global_sums[ct] = ct_sum
            global_counts[ct] = n_cells_in_pool
        else:
            global_sums[ct] += ct_sum
            global_counts[ct] += n_cells_in_pool
            
    # Clean up memory before loading the next pool
    del adata
    gc.collect()

print("Finished processing all pools. Calculating global means...")

# 3. Calculate Global Means (Accessibility Frequency)
means_dict = {}
for ct in global_sums.keys():
    print(f"Calculating mean for cell type: {ct}")
    total_cells = global_counts[ct]
    means_dict[ct] = global_sums[ct] / total_cells

# 4. Create the final Pseudobulk AnnData object
df_means = pd.DataFrame(means_dict, index=peak_names)

cell_counts = [global_counts[ct] for ct in df_means.columns]
obs_df = pd.DataFrame({"n_cells": cell_counts}, index=df_means.columns)

adata_pseudo = sc.AnnData(
    X=df_means.T.values,
    obs=obs_df,
    var=pd.DataFrame(index=df_means.index)
)

print(f"Final pseudobulk shape: {adata_pseudo.shape}")

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# =============================================================================
# Step 1: Compute mean accessibility per peak per cell type
# =============================================================================

def compute_celltype_mean_accessibility(
    adata: sc.AnnData,
    celltype_col: str = "cell_type",
    is_pseudobulk: bool = False,
    min_cells_per_type: int = 50,
) -> pd.DataFrame:
    """
    Compute mean chromatin accessibility per peak per cell type.
    
    Parameters
    ----------
    adata : AnnData
        scATAC-seq AnnData object.
        - If cell-level: obs are cells, var are peaks, X is counts/accessibility.
        - If pseudobulk: obs are cell types (or donor x cell type), var are peaks.
    celltype_col : str
        Column in adata.obs containing cell type labels (cell-level mode).
    is_pseudobulk : bool
        If True, treat adata as already pseudobulked (obs = cell types, X = mean accessibility).
    min_cells_per_type : int
        Minimum number of cells required per cell type (cell-level mode).
    
    Returns
    -------
    pd.DataFrame
        Shape (n_peaks, n_celltypes), mean accessibility per peak per cell type.
    """
    if is_pseudobulk:
        # If already pseudobulk, obs index should be cell types
        logger.info(f"Using pseudobulk mode. Shape: {adata.shape}")
        X = adata.X
        if sparse.issparse(X):
            X = X.toarray()
        mean_df = pd.DataFrame(
            X, index=adata.obs_names, columns=adata.var_names
        ).T  # peaks x cell_types
        return mean_df

    # Cell-level mode: compute mean per cell type
    logger.info(f"Computing mean accessibility from cell-level data. Shape: {adata.shape}")
    assert celltype_col in adata.obs.columns, (
        f"Column '{celltype_col}' not found in adata.obs. "
        f"Available columns: {list(adata.obs.columns)}"
    )
    
    cell_types = adata.obs[celltype_col].unique()
    logger.info(f"Found {len(cell_types)} cell types")
    
    # Filter cell types with too few cells
    ct_counts = adata.obs[celltype_col].value_counts()
    valid_cts = ct_counts[ct_counts >= min_cells_per_type].index.tolist()
    logger.info(
        f"Keeping {len(valid_cts)} cell types with >= {min_cells_per_type} cells "
        f"(dropped: {set(cell_types) - set(valid_cts)})"
    )
    
    mean_dict = {}
    for ct in valid_cts:
        mask = adata.obs[celltype_col] == ct
        X_ct = adata[mask].X
        if sparse.issparse(X_ct):
            X_ct = X_ct.toarray()
        mean_dict[ct] = np.mean(X_ct, axis=0).flatten()
    
    mean_df = pd.DataFrame(mean_dict, index=adata.var_names)  # peaks x cell_types
    logger.info(f"Mean accessibility matrix shape: {mean_df.shape}")
    return mean_df


# =============================================================================
# Step 2: Compute Shannon entropy and specificity scores
# =============================================================================

def compute_specificity_scores(
    mean_accessibility: pd.DataFrame,
    pseudocount: float = 1e-10,
) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Compute Shannon entropy-based specificity scores Q(p|t) for each peak.
    
    Following Schug et al. (2005):
        p_i = q_i / sum(q_i)          # convert to probabilities
        H(p) = -sum(p_t * log2(p_t))  # Shannon entropy
        Q(p|t) = H(p) - log2(p_t)     # specificity score (lower = more specific to t)
    
    Parameters
    ----------
    mean_accessibility : pd.DataFrame
        Shape (n_peaks, n_celltypes). Mean accessibility values.
    pseudocount : float
        Small value added to avoid log(0).
    
    Returns
    -------
    Q_scores : pd.DataFrame
        Shape (n_peaks, n_celltypes). Specificity scores Q(p|t).
    H_scores : pd.Series
        Shape (n_peaks,). Entropy H(p) for each peak.
    """
    logger.info("Computing Shannon entropy-based specificity scores...")
    
    # q_i = relative accessibility scores (ensure non-negative)
    q = mean_accessibility.values.copy()
    q = np.maximum(q, 0)  # ensure non-negative
    q = q + pseudocount   # avoid division by zero
    
    # Convert to probabilities: p_i = q_i / sum(q_i)
    row_sums = q.sum(axis=1, keepdims=True)
    p = q / row_sums  # (n_peaks, n_celltypes)
    
    # Shannon entropy: H(p) = -sum(p_t * log2(p_t))
    log2_p = np.log2(p)
    H = -np.sum(p * log2_p, axis=1)  # (n_peaks,)
    
    # Specificity score: Q(p|t) = H(p) - log2(p_t)
    # Lower Q means more specific to cell type t
    Q = H[:, np.newaxis] - log2_p  # (n_peaks, n_celltypes)
    
    Q_scores = pd.DataFrame(Q, index=mean_accessibility.index, columns=mean_accessibility.columns)
    H_scores = pd.Series(H, index=mean_accessibility.index, name="entropy")
    
    logger.info(
        f"Entropy range: [{H.min():.3f}, {H.max():.3f}], "
        f"max possible: {np.log2(mean_accessibility.shape[1]):.3f}"
    )
    
    return Q_scores, H_scores


# =============================================================================
# Step 3: Estimate null distribution and compute p-values
# =============================================================================

def estimate_null_distribution(
    mean_accessibility: pd.DataFrame,
    n_celltypes: int,
    n_simulations: int = 500_000,
    pseudocount: float = 1e-10,
    random_seed: int = 42,
) -> np.ndarray:
    """
    Estimate the null distribution of minimum Q scores via simulation.
    
    Under the null hypothesis:
    - Each peak has an average accessibility level across all cell types
    - log2(fold changes from average) ~ N(0, s)
    - s is estimated from the top 50% least variable peaks
    
    Parameters
    ----------
    mean_accessibility : pd.DataFrame
        Shape (n_peaks, n_celltypes). Used to estimate s.
    n_celltypes : int
        Number of cell types.
    n_simulations : int
        Number of samples to draw for the empirical null distribution.
    pseudocount : float
        Small value to avoid log(0).
    random_seed : int
        Random seed for reproducibility.
    
    Returns
    -------
    null_Q_min : np.ndarray
        Shape (n_simulations,). Empirical null distribution of minimum Q(p|t) scores.
    """
    logger.info("Estimating null distribution of specificity scores...")
    
    # Estimate s from the top 50% least variable peaks
    # Compute log2 fold changes from the mean for each peak
    q = mean_accessibility.values.copy()
    q = np.maximum(q, pseudocount)
    
    row_means = q.mean(axis=1, keepdims=True)
    log2_fc = np.log2(q / row_means)  # log2 fold changes from the average
    
    # Compute variance of log2 fold changes per peak
    peak_var = np.var(log2_fc, axis=1)
    
    # Select top 50% least variable peaks
    var_threshold = np.median(peak_var)
    least_variable_mask = peak_var <= var_threshold
    n_least_var = least_variable_mask.sum()
    logger.info(f"Using {n_least_var} least variable peaks (50%) to estimate s")
    
    # Estimate s as the std of log2 fold changes from least variable peaks
    s = np.std(log2_fc[least_variable_mask])
    logger.info(f"Estimated s (std of log2 fold changes under null): {s:.4f}")
    
    # Simulate null distribution
    rng = np.random.default_rng(random_seed)
    
    # Draw log2 fold changes from N(0, s)
    null_log2_fc = rng.normal(0, s, size=(n_simulations, n_celltypes))  # (n_sim, n_ct)
    
    # Convert to relative accessibility: q_i = 2^(X_i)
    null_q = np.power(2, null_log2_fc)
    null_q = null_q + pseudocount
    
    # Convert to probabilities
    null_p = null_q / null_q.sum(axis=1, keepdims=True)
    
    # Compute entropy
    null_log2_p = np.log2(null_p)
    null_H = -np.sum(null_p * null_log2_p, axis=1)  # (n_sim,)
    
    # Compute Q(p|t) for each simulation
    null_Q = null_H[:, np.newaxis] - null_log2_p  # (n_sim, n_ct)
    
    # The minimum Q across cell types is the "most specific" score
    null_Q_min = null_Q.min(axis=1)  # (n_sim,)
    
    logger.info(
        f"Null Q_min distribution: mean={null_Q_min.mean():.3f}, "
        f"std={null_Q_min.std():.3f}, "
        f"5th percentile={np.percentile(null_Q_min, 5):.3f}"
    )
    
    return null_Q_min


def compute_pvalues(
    Q_scores: pd.DataFrame,
    null_Q_min: np.ndarray,
) -> pd.Series:
    """
    Compute p-values for each peak's specificity score against the null distribution.
    
    Parameters
    ----------
    Q_scores : pd.DataFrame
        Shape (n_peaks, n_celltypes). Observed specificity scores.
    null_Q_min : np.ndarray
        Shape (n_simulations,). Empirical null distribution of minimum Q scores.
    
    Returns
    -------
    pvalues : pd.Series
        Shape (n_peaks,). P-values for each peak.
    """
    logger.info("Computing p-values against null distribution...")
    
    # For each peak, the observed minimum Q across cell types
    obs_Q_min = Q_scores.min(axis=1).values  # (n_peaks,)
    
    # Sort null distribution for efficient p-value computation
    null_sorted = np.sort(null_Q_min)
    n_null = len(null_sorted)
    
    # P-value = fraction of null samples with Q_min <= observed Q_min
    # (one-sided test: lower Q means more specific)
    pvals = np.searchsorted(null_sorted, obs_Q_min) / n_null
    
    pvalues = pd.Series(pvals, index=Q_scores.index, name="pvalue")
    
    n_sig_001 = (pvalues < 0.01).sum()
    n_sig_025 = (pvalues < 0.025).sum()
    n_sig_05 = (pvalues < 0.05).sum()
    logger.info(
        f"Significant peaks: p<0.01: {n_sig_001}, p<0.025: {n_sig_025}, p<0.05: {n_sig_05}"
    )
    
    return pvalues


# =============================================================================
# Step 4: Identify cell-type-restricted peaks
# =============================================================================

def identify_restricted_peaks(
    Q_scores: pd.DataFrame,
    pvalues: pd.Series,
    p_cutoff: float = 0.025,
) -> pd.DataFrame:
    """
    Identify cell-type-restricted peaks and assign each to its most-specific cell type.
    
    Parameters
    ----------
    Q_scores : pd.DataFrame
        Shape (n_peaks, n_celltypes). Specificity scores.
    pvalues : pd.Series
        Shape (n_peaks,). P-values.
    p_cutoff : float
        P-value cutoff for significance.
    
    Returns
    -------
    restricted_peaks : pd.DataFrame
        DataFrame with columns: peak, assigned_celltype, Q_min, pvalue
    """
    sig_mask = pvalues < p_cutoff
    sig_peaks = Q_scores.loc[sig_mask]
    sig_pvals = pvalues.loc[sig_mask]
    
    # Assign each peak to the cell type with the minimum Q score
    assigned_ct = sig_peaks.idxmin(axis=1)
    Q_min = sig_peaks.min(axis=1)
    
    restricted_peaks = pd.DataFrame({
        "peak": sig_peaks.index,
        "assigned_celltype": assigned_ct.values,
        "Q_min": Q_min.values,
        "pvalue": sig_pvals.values,
    })
    
    logger.info(f"Identified {len(restricted_peaks)} cell-type-restricted peaks (p < {p_cutoff})")
    logger.info("Peaks per cell type:")
    ct_counts = restricted_peaks["assigned_celltype"].value_counts()
    for ct, count in ct_counts.items():
        logger.info(f"  {ct}: {count}")
    
    return restricted_peaks


# =============================================================================
# Step 5: Generate heatmap
# =============================================================================

def plot_heatmap(
    mean_accessibility: pd.DataFrame,
    restricted_peaks: pd.DataFrame,
    celltype_order: Optional[list] = None,
    figsize: Tuple[int, int] = (16, 20),
    cmap: str = "YlOrRd",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    log2_transform: bool = True,
    row_normalize: bool = False,
    title: str = "Cell-type-restricted chromatin accessibility peaks",
    save_path: Optional[str] = None,
    show_celltype_boundaries: bool = True,
    dpi: int = 150,
) -> plt.Figure:
    """
    Generate a heatmap of cell-type-restricted peaks.
    
    Rows = peaks (grouped by assigned cell type, sorted by specificity within each group)
    Columns = cell types
    Color = log2-transformed chromatin accessibility
    
    Parameters
    ----------
    mean_accessibility : pd.DataFrame
        Shape (n_peaks, n_celltypes). Mean accessibility values.
    restricted_peaks : pd.DataFrame
        Output from identify_restricted_peaks().
    celltype_order : list, optional
        Custom ordering of cell types along the x-axis. If None, uses hierarchical
        clustering of cell types based on their accessibility profiles.
    figsize : tuple
        Figure size.
    cmap : str
        Colormap.
    vmin, vmax : float, optional
        Color scale limits. If None, auto-determined.
    log2_transform : bool
        Whether to log2-transform accessibility values.
    row_normalize : bool
        Whether to z-score normalize each row (peak) across cell types.
    title : str
        Plot title.
    save_path : str, optional
        Path to save the figure.
    show_celltype_boundaries : bool
        Whether to draw horizontal lines between cell type groups.
    dpi : int
        Figure resolution.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    logger.info("Generating Figure 2E-style heatmap...")
    
    # Get the accessibility matrix for restricted peaks only
    restricted_peak_names = restricted_peaks["peak"].values
    mat = mean_accessibility.loc[restricted_peak_names].copy()
    
    # Determine cell type order
    if celltype_order is None:
        # Use hierarchical clustering of cell types
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist
        
        ct_mat = mat.T.values  # (n_celltypes, n_peaks)
        if log2_transform:
            ct_mat_for_clust = np.log2(np.maximum(ct_mat, 1e-10) + 1)
        else:
            ct_mat_for_clust = ct_mat.copy()
        
        dist = pdist(ct_mat_for_clust, metric="correlation")
        Z = linkage(dist, method="ward")
        order_idx = leaves_list(Z)
        celltype_order = [mat.columns[i] for i in order_idx]
        logger.info(f"Cell type order (hierarchical clustering): {celltype_order}")
    
    # Reorder columns
    mat = mat[celltype_order]
    
    # Order peaks: group by assigned cell type (in column order), then sort by Q_min within group
    peak_order = []
    group_boundaries = []  # for drawing horizontal lines
    group_labels = []
    
    for ct in celltype_order:
        ct_peaks = restricted_peaks[restricted_peaks["assigned_celltype"] == ct].copy()
        if len(ct_peaks) == 0:
            continue
        ct_peaks = ct_peaks.sort_values("Q_min", ascending=True)
        peak_order.extend(ct_peaks["peak"].tolist())
        group_boundaries.append(len(peak_order))
        group_labels.append(ct)
    
    mat = mat.loc[peak_order]
    
    # Transform values
    plot_mat = mat.values.copy()
    if log2_transform:
        plot_mat = np.log2(np.maximum(plot_mat, 1e-10) + 1)
    
    if row_normalize:
        row_means = plot_mat.mean(axis=1, keepdims=True)
        row_stds = plot_mat.std(axis=1, keepdims=True)
        row_stds[row_stds == 0] = 1
        plot_mat = (plot_mat - row_means) / row_stds
    
    # Auto-determine color limits
    if vmin is None:
        vmin = np.percentile(plot_mat, 1)
    if vmax is None:
        vmax = np.percentile(plot_mat, 99)
    
    # Create figure
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 2, width_ratios=[30, 1], wspace=0.02)
    
    ax_heat = fig.add_subplot(gs[0])
    ax_cbar = fig.add_subplot(gs[1])
    
    # Plot heatmap using imshow for efficiency with large matrices
    im = ax_heat.imshow(
        plot_mat,
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="none",
    )
    
    # Draw cell type group boundaries
    if show_celltype_boundaries:
        for boundary in group_boundaries[:-1]:
            ax_heat.axhline(y=boundary - 0.5, color="white", linewidth=0.5, alpha=0.8)
    
    # # Add cell type group labels on the right
    # prev_boundary = 0
    # for i, (boundary, label) in enumerate(zip(group_boundaries, group_labels)):
    #     mid = (prev_boundary + boundary) / 2
    #     ax_heat.text(
    #         len(celltype_order) + 0.5, mid, label,
    #         fontsize=5, va="center", ha="left",
    #         clip_on=False,
    #     )
    #     prev_boundary = boundary
    
    # X-axis: cell type labels
    ax_heat.set_xticks(range(len(celltype_order)))
    ax_heat.set_xticklabels(celltype_order, rotation=90, fontsize=9, ha="center")
    ax_heat.tick_params(axis="x", bottom=True, top=False, labelbottom=True)
    
    # Y-axis: hide individual peak labels (too many)
    ax_heat.set_yticks([])
    ax_heat.set_ylabel(f"Cell-type-restricted peaks (n={len(peak_order):,})", fontsize=12)
    
    # Colorbar
    cbar_label = "log₂(accessibility + 1)" if log2_transform else "Accessibility"
    if row_normalize:
        cbar_label = "Z-score"
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label(cbar_label, fontsize=12) 
    cbar.ax.tick_params(labelsize=10)
    
    ax_heat.set_title(title, fontsize=16, pad=15)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        logger.info(f"Figure saved to {save_path}")
    
    return fig


# =============================================================================
# Step 6: Summary statistics and additional plots
# =============================================================================

def plot_summary_stats(
    restricted_peaks: pd.DataFrame,
    H_scores: pd.Series,
    pvalues: pd.Series,
    save_path: Optional[str] = None,
    dpi: int = 150,
) -> plt.Figure:
    """
    Plot summary statistics of the entropy analysis.
    
    Parameters
    ----------
    restricted_peaks : pd.DataFrame
        Output from identify_restricted_peaks().
    H_scores : pd.Series
        Entropy scores for all peaks.
    pvalues : pd.Series
        P-values for all peaks.
    save_path : str, optional
        Path to save the figure.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Histogram of entropy scores
    ax = axes[0, 0]
    ax.hist(H_scores, bins=100, color="steelblue", edgecolor="none", alpha=0.8)
    ax.set_xlabel("Shannon entropy H(p)")
    ax.set_ylabel("Number of peaks")
    ax.set_title("Distribution of peak entropy scores")
    ax.axvline(H_scores.median(), color="red", linestyle="--", label=f"Median: {H_scores.median():.2f}")
    ax.legend()
    
    # 2. Histogram of p-values
    ax = axes[0, 1]
    ax.hist(pvalues, bins=100, color="coral", edgecolor="none", alpha=0.8)
    ax.set_xlabel("P-value")
    ax.set_ylabel("Number of peaks")
    ax.set_title("Distribution of specificity p-values")
    ax.axvline(0.025, color="red", linestyle="--", label="p = 0.025")
    ax.legend()
    
    # 3. Bar plot of peaks per cell type
    ax = axes[1, 0]
    ct_counts = restricted_peaks["assigned_celltype"].value_counts().sort_values(ascending=True)
    ax.barh(range(len(ct_counts)), ct_counts.values, color="teal", edgecolor="none")
    ax.set_yticks(range(len(ct_counts)))
    ax.set_yticklabels(ct_counts.index, fontsize=7)
    ax.set_xlabel("Number of restricted peaks")
    ax.set_title("Cell-type-restricted peaks per cell type")
    
    # 4. Distribution of Q_min for significant peaks
    ax = axes[1, 1]
    ax.hist(restricted_peaks["Q_min"], bins=50, color="mediumpurple", edgecolor="none", alpha=0.8)
    ax.set_xlabel("Minimum Q(p|t) score")
    ax.set_ylabel("Number of peaks")
    ax.set_title("Specificity scores of restricted peaks")
    
    plt.suptitle("Entropy-based cell-type specificity analysis summary", fontsize=14, y=1.02)
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        logger.info(f"Summary figure saved to {save_path}")
    
    return fig


# =============================================================================
# Main pipeline function
# =============================================================================

def run_entropy_pipeline(
    adata: sc.AnnData,
    celltype_col: str = "cell_type",
    is_pseudobulk: bool = False,
    min_cells_per_type: int = 50,
    p_cutoff: float = 0.025,
    n_simulations: int = 500_000,
    pseudocount: float = 1e-10,
    random_seed: int = 42,
    celltype_order: Optional[list] = None,
    log2_transform: bool = True,
    row_normalize: bool = False,
    cmap: str = "YlOrRd",
    figsize: Tuple[int, int] = (16, 20),
    save_prefix: Optional[str] = None,
    dpi: int = 150,
    show_celltype_boundaries: bool = True
) -> Dict:
    """
    Full pipeline: compute entropy-based specificity, identify restricted peaks,
    and generate Figure 2E-style heatmap.
    
    Parameters
    ----------
    adata : AnnData
        scATAC-seq AnnData object.
    celltype_col : str
        Column in adata.obs containing cell type labels.
    is_pseudobulk : bool
        If True, treat adata as already pseudobulked.
    min_cells_per_type : int
        Minimum cells per cell type (cell-level mode).
    p_cutoff : float
        P-value cutoff for identifying restricted peaks.
    n_simulations : int
        Number of simulations for null distribution.
    pseudocount : float
        Pseudocount to avoid log(0).
    random_seed : int
        Random seed.
    celltype_order : list, optional
        Custom ordering of cell types for the heatmap.
    log2_transform : bool
        Whether to log2-transform for heatmap visualization.
    row_normalize : bool
        Whether to z-score normalize rows in the heatmap.
    cmap : str
        Colormap for heatmap.
    figsize : tuple
        Figure size for heatmap.
    save_prefix : str, optional
        Prefix for saving output files. If None, no files saved.
    dpi : int
        Figure resolution.
    
    Returns
    -------
    result : dict
        Dictionary containing:
        - 'mean_accessibility': pd.DataFrame (peaks x cell types)
        - 'Q_scores': pd.DataFrame (peaks x cell types)
        - 'H_scores': pd.Series (entropy per peak)
        - 'pvalues': pd.Series (p-value per peak)
        - 'restricted_peaks': pd.DataFrame (significant peaks + assignments)
        - 'null_Q_min': np.ndarray (null distribution)
        - 'fig_heatmap': matplotlib.figure.Figure
        - 'fig_summary': matplotlib.figure.Figure
    """
    logger.info("=" * 70)
    logger.info("SHANNON ENTROPY-BASED CELL-TYPE SPECIFICITY ANALYSIS")
    logger.info("=" * 70)
    
    # Step 1: Compute mean accessibility
    mean_acc = compute_celltype_mean_accessibility(
        adata, celltype_col=celltype_col,
        is_pseudobulk=is_pseudobulk,
        min_cells_per_type=min_cells_per_type,
    )
    
    # Step 2: Compute specificity scores
    Q_scores, H_scores = compute_specificity_scores(mean_acc, pseudocount=pseudocount)
    
    # Step 3: Estimate null distribution
    n_celltypes = mean_acc.shape[1]
    null_Q_min = estimate_null_distribution(
        mean_acc, n_celltypes=n_celltypes,
        n_simulations=n_simulations,
        pseudocount=pseudocount,
        random_seed=random_seed,
    )
    
    # Step 4: Compute p-values
    pvalues = compute_pvalues(Q_scores, null_Q_min)
    
    # Step 5: Identify restricted peaks
    restricted_peaks = identify_restricted_peaks(Q_scores, pvalues, p_cutoff=p_cutoff)
    
    # Step 6: Generate heatmap
    heatmap_save = f"{save_prefix}_heatmap.png" if save_prefix else None
    fig_heatmap = plot_heatmap(
        mean_acc, restricted_peaks,
        celltype_order=celltype_order,
        figsize=figsize, cmap=cmap,
        log2_transform=log2_transform,
        row_normalize=row_normalize,
        save_path=heatmap_save,
        dpi=dpi,
        show_celltype_boundaries=show_celltype_boundaries
    )
    
    # Step 7: Summary plots
    summary_save = f"{save_prefix}_summary.png" if save_prefix else None
    fig_summary = plot_summary_stats(
        restricted_peaks, H_scores, pvalues,
        save_path=summary_save,
        dpi=dpi,
    )
    
    # Save tables if requested
    if save_prefix:
        restricted_peaks.to_csv(f"{save_prefix}_restricted_peaks.csv", index=False)
        Q_scores.to_csv(f"{save_prefix}_Q_scores.csv")
        pvalues.to_csv(f"{save_prefix}_pvalues.csv")
        logger.info(f"Tables saved with prefix: {save_prefix}")
    
    result = {
        "mean_accessibility": mean_acc,
        "Q_scores": Q_scores,
        "H_scores": H_scores,
        "pvalues": pvalues,
        "restricted_peaks": restricted_peaks,
        "null_Q_min": null_Q_min,
        "fig_heatmap": fig_heatmap,
        "fig_summary": fig_summary,
    }
    
    logger.info("=" * 70)
    logger.info(f"ANALYSIS COMPLETE: {len(restricted_peaks)} restricted peaks identified")
    logger.info("=" * 70)
    
    return result

# Process the pseudobulk object
output_prefix = "/g/data/fy54/od8037/TenK10K/HeatmapNew/Results/tenk10k"
adata_pseudo = adata_pseudo[~adata_pseudo.obs_names.isin(["ILC", "CD8_Proliferating"])].copy()
print(f"Loaded object with {adata_pseudo.shape[0]} cell types and {adata_pseudo.shape[1]} peaks.")
max_accessibility_per_peak = np.asarray(adata_pseudo.X.max(axis=0)).flatten()

# Keep peaks that are open in at least 5% (0.05) of cells in at least one cell type
keep_peaks_mask = max_accessibility_per_peak > 0.05
adata_filtered = adata_pseudo[:, keep_peaks_mask].copy()

print(f"Original peaks: {adata_pseudo.shape[1]}")
print(f"Filtered peaks: {adata_filtered.shape[1]}")

celltype_order = [
    "CD4_TCM", "CD4_Naive", "CD4_TEM", "CD4_CTL", "Treg", "CD4_Proliferating",
    "gdT", "MAIT", "dnT",
    "CD8_TEM", "CD8_Naive", "CD8_TCM",
    "NK", "NK_CD56bright", "NK_Proliferating",
    "B_naive", "B_intermediate", "B_memory", "Plasmablast",
    "CD14_Mono", "CD16_Mono",
    "cDC2", "pDC", "cDC1", "ASDC", "HSPC"
]

print("Starting entropy pipeline...")
results = run_entropy_pipeline(
    adata=adata_filtered,
    is_pseudobulk=True,
    p_cutoff=0.025,
    n_simulations=500_000,
    pseudocount=1e-10,
    celltype_order=celltype_order, 
    log2_transform=True,           
    cmap="viridis",     
    figsize=(14, 8),
    save_prefix=output_prefix, 
    dpi=300,
    show_celltype_boundaries=False
)

# =============================================================================
# Broad lineage version
# =============================================================================
coarse_mapping = {
    # CD4 T cells
    "CD4_TCM": "CD4_T", 
    "CD4_Naive": "CD4_T", 
    "CD4_TEM": "CD4_T", 
    "CD4_CTL": "CD4_T", 
    "Treg": "CD4_T", 
    "CD4_Proliferating": "CD4_T",
    
    # Unconventional T cells
    "gdT": "Unconven_T", 
    "MAIT": "Unconven_T", 
    "dnT": "Unconven_T",
    
    # CD8 T cells
    "CD8_TEM": "CD8_T", 
    "CD8_Naive": "CD8_T", 
    "CD8_TCM": "CD8_T",
    
    # Natural Killer cells
    "NK": "NK", 
    "NK_CD56bright": "NK", 
    "NK_Proliferating": "NK",
    
    # B cells
    "B_naive": "B", 
    "B_intermediate": "B", 
    "B_memory": "B", 
    "Plasmablast": "B",
    
    # Monocytes
    "CD14_Mono": "Mono", 
    "CD16_Mono": "Mono",
    
    # Dendritic cells
    "cDC2": "DC", 
    "pDC": "DC", 
    "cDC1": "DC", 
    "ASDC": "DC",
    
    # Hematopoietic Stem/Progenitor Cells
    "HSPC": "HSPC"
}
new_groups = [coarse_mapping.get(ct, ct) for ct in adata_pseudo.obs_names]

df_pseudo = pd.DataFrame(adata_pseudo.X, index=adata_pseudo.obs_names, columns=adata_pseudo.var_names)
cell_counts = adata_pseudo.obs['n_cells']

# Step A: Multiply mean frequencies by cell counts to get raw total hits
df_raw_sums = df_pseudo.mul(cell_counts, axis=0)
df_raw_sums["broad_lineage"] = new_groups

# Step B: Sum the raw hits by lineage
lineage_sums = df_raw_sums.groupby("broad_lineage").sum()

# Step C: Sum the total cells by lineage
obs_df = adata_pseudo.obs.copy()
obs_df["broad_lineage"] = new_groups
lineage_total_cells = obs_df.groupby("broad_lineage")['n_cells'].sum()

# Step D: Divide the summed hits by the summed cells to get the weighted mean
df_coarse = lineage_sums.div(lineage_total_cells, axis=0)

adata_coarse = sc.AnnData(
    X=df_coarse.values,
    obs=pd.DataFrame({"n_cells": lineage_total_cells}),
    var=adata_pseudo.var.copy()
)
max_accessibility_per_peak = np.asarray(adata_coarse.X.max(axis=0)).flatten()

# Keep peaks that are open in at least 5% (0.05) of cells in at least one cell type
keep_peaks_mask = max_accessibility_per_peak > 0.05
adata_filtered = adata_coarse[:, keep_peaks_mask].copy()

results = run_entropy_pipeline(
    adata=adata_filtered,
    is_pseudobulk=True,             
    p_cutoff=0.025,               
    n_simulations=500_000,        
    pseudocount=1e-10,
    log2_transform=True,         
    cmap="viridis",  
    figsize=(6, 4),
    save_prefix="/g/data/fy54/od8037/TenK10K/HeatmapNew/Results/broad",      
    dpi=300,
    show_celltype_boundaries=False
)
