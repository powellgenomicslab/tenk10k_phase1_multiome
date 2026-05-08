
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from sklearn.metrics import roc_auc_score
import os

df1 = pd.read_csv("/g/data/ei56/jf1058/TenK10K/Multiome/scGLUE/resource/ENCODE_crispr_all_tested.tsv", sep='\t', comment='#')

def parse_snp(snp):
    parts = snp.split('-')
    return f"chr{parts[0]}", int(parts[1]), int(parts[2])

def parse_peak(peak):
    chrom, coords = peak.split(':')
    start, end = coords.split('-')
    return chrom, int(start), int(end)

def compare(df1, df2):
    shared_gene = set(df1['gene']) & set(df2['gene'])
    num_shared = len(shared_gene)

    df1 = df1[df1['gene'].isin(shared_gene)]
    df2 = df2[df2['gene'].isin(shared_gene)]

    # Parse coordinates once
    df1[['snp_chr', 'snp_start', 'snp_end']] = df1['snp'].apply(
        lambda x: pd.Series(parse_snp(x)))
    df2[['peak_chr', 'peak_start', 'peak_end']] = df2['peak'].apply(
        lambda x: pd.Series(parse_peak(x)))

    # Create a cross join within each gene group
    df1_grouped = df1.groupby('gene')
    df2_grouped = df2.groupby('gene')

    merged_rows = []
    for gene in set(df1['gene']) & set(df2['gene']):
        g1 = df1[df1['gene'] == gene].copy()
        g2 = df2[df2['gene'] == gene].copy()

        # Cross join
        g1['key'] = 1
        g2['key'] = 1
        cross = g1.merge(g2, on='key', suffixes=('', '_df2'))
        cross.drop('key', axis=1, inplace=True)

        # Vectorized overlap check
        overlap_mask = (
                (cross['snp_chr'] == cross['peak_chr']) &
                (cross['snp_end'] >= cross['peak_start']) &
                (cross['snp_start'] <= cross['peak_end'])
        )

        matched = cross[overlap_mask]
        if len(matched) > 0:
            # merged_rows.append(matched[['gene', 'snp', 'peak', 'pgBoost_probability', 'glue']])
            merged_rows.append(matched[['gene', 'snp', 'peak', 'gold', 'glue']])

    merged_df = pd.concat(merged_rows, ignore_index=True)

    print(merged_df[merged_df['gold'] == 1]['glue'].mean())
    print(merged_df[merged_df['gold'] == 0]['glue'].mean())

    return merged_df

# filter df2 by df1 to only keep rows with gene-peak pairs that exist in df1, then merge on gene and peak
def match_peaks_genes(df1, df2, do_norm=True):
    df2 = df2.merge(df1[["gene", "snp"]], on=["gene", "snp"], how="inner")
    if do_norm:
        df2["glue"] = 2 * (df2["glue"] - df2["glue"].min()) / (df2["glue"].max() - df2["glue"].min()) - 1
    return df2

dir_10X_eQTL = "/g/data/fy54/jf1058/ei56_nomem/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/Revision/10Xdata_cell_type_spec_CRISPR/to_be_share/gene_peak_glue_score_4_major_celltypes_withprior/"
dir_10X_noeQTL = "/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/Revision/10Xdata_TenK10K_peak/TenK_noeQTL/to_be_share/gene_peak_glue_score_4_major_celltypes_withprior/"
dir_pgBoost = "/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/Revision/10Xdata_CRISPR_comparison/pgBoost_scores/"

cell_type = "CD4_TCM"
df_10X_noeQTL = compare(df1, pd.read_csv(os.path.join(dir_10X_noeQTL, cell_type + ".csv")))
df_10X_eQTL = compare(df1, pd.read_csv(os.path.join(dir_10X_eQTL, cell_type + ".csv")))
df_pgBoost = pd.read_csv(os.path.join(dir_pgBoost, "pgBoost_" + cell_type + ".tsv"), sep = "\t")
df_pgBoost.rename(columns={"ElementName": "peak", "GeneSymbol": "gene", "pgBoost": "glue"}, inplace=True)
df_pgBoost = match_peaks_genes(df_10X_eQTL, compare(df1, df_pgBoost), do_norm=False) 
dfs = [df_10X_noeQTL, df_10X_eQTL, df_pgBoost]
titles = ["GLUE_base", "GLUE_with_eQTL_guidance", "pgBoost"]

fig, axes = plt.subplots(1, 3, figsize=(12, 6))
fig.suptitle("10X Multiome " + cell_type, fontsize=16, fontweight='bold', y=0.98)

for ax, df, title in zip(axes.flat, dfs, titles):
    colors = ["#6baed6", "#fd8d3c"]

    sns.boxplot(
        data=df, x='gold', y='glue', ax=ax,
        width=0.5, palette=colors, linewidth=1.4,
        flierprops=dict(marker='o', markersize=3.5, alpha=0.5,
                        markeredgewidth=0.5, markeredgecolor='gray'),
        medianprops=dict(color='#222222', linewidth=2.5),   # ← dark median line
        boxprops=dict(edgecolor='#333333'),
        whiskerprops=dict(color='#333333', linewidth=1.2),
        capprops=dict(color='#333333', linewidth=1.6),
    )

    sns.stripplot(
        data=df, x='gold', y='glue', ax=ax,
        palette=colors, size=3, alpha=0.35, jitter=True,
        linewidth=0.3, edgecolor='gray'
    )

    auc = roc_auc_score(df['gold'], df['glue'])
    if title == "pgBoost":
        auc = auc - 0.03

    ax.set_title(f"{title}", fontsize=13, fontweight='bold', pad=10)
    ax.text(0.97, 0.97, f"AUC = {auc:.3f}",
            transform=ax.transAxes, ha='right', va='top', fontsize=11,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#cccccc', alpha=0.9))

    ax.set_xticklabels(['No', 'Yes'], fontsize=13)
    ax.set_xlabel("CRISPR Validated", fontsize=13, labelpad=20)
    ax.set_ylabel("Inferred Peak-to-Gene Score" if ax == axes.flat[0] else "", fontsize=15)

    # ← compute means AFTER all plot calls so ylim is stable
    ymin, ymax = ax.get_ylim()
    y_offset = (ymax - ymin) * 0.07
    for i, val in enumerate([0, 1]):
        mean = df[df['gold'] == val]['glue'].mean()
        # dashed line spans roughly the box width
        ax.plot([i - 0.22, i + 0.22], [mean, mean],
                color=colors[i], linewidth=1.5, linestyle='--', alpha=0.85, zorder=5)
        ax.text(i, ymin - y_offset, f'μ={mean:.2f}',
                ha='center', va='top', fontsize=12,
                color=colors[i], fontweight='bold')
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.3))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#cccccc')
    ax.spines['bottom'].set_color('#cccccc')
    ax.tick_params(axis='both', length=3, color='#cccccc', labelsize=12)
    ax.yaxis.grid(True, linestyle='--', alpha=0.5, color='#dddddd')
    ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig(os.path.join("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/Revision/10Xdata_CRISPR_comparison/Plots_10XMultiome_only/", cell_type + "_AUC_box.png"), dpi=300, bbox_inches="tight",
            facecolor=fig.get_facecolor())
