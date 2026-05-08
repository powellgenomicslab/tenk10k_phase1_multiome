#!/usr/bin/env python3

import argparse
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy import stats

# ─── Genotype extraction ─────────────────────────────────────────────────────

def extract_genotype_dosage(plink_prefix: str, variant_id: str, plink_bin: str) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = os.path.join(tmpdir, "geno")

        # Write variant ID to a filter file
        snp_file = os.path.join(tmpdir, "snp.txt")
        with open(snp_file, "w") as f:
            f.write(variant_id + "\n")

        cmd = [
            plink_bin,  # pass plink 1.9 binary here (e.g. "plink")
            "--bfile", plink_prefix,
            "--extract", snp_file,
            "--recode", "A",  # additive coding: 0/1/2 copies of A2
            "--out", out_prefix,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"plink failed for variant '{variant_id}'.\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

        raw_file = out_prefix + ".raw"
        if not os.path.exists(raw_file):
            raise FileNotFoundError(
                f"plink did not produce {raw_file}. "
                f"Check that variant ID '{variant_id}' exists in {plink_prefix}.bim"
            )

        geno_df = pd.read_csv(raw_file, sep=r"\s+")

    # Plink 1.9 names the dosage column  <SNPID>_<A2>  e.g. "15:19776288:T:TTC_TTC"
    meta_cols = {"FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"}
    dosage_cols = [c for c in geno_df.columns if c not in meta_cols]
    if len(dosage_cols) == 0:
        raise ValueError(f"No dosage column found in plink output. Columns: {list(geno_df.columns)}")
    if len(dosage_cols) > 1:
        print(f"[WARNING] Multiple dosage columns found: {dosage_cols}. Using first: {dosage_cols[0]}")

    dosage_col = dosage_cols[0]
    result_df = geno_df[["IID", dosage_col]].copy()
    result_df = result_df.rename(columns={dosage_col: "dosage", "IID": "donor_id"})

    # Round to integer genotype (0/1/2), drop missing
    result_df["dosage"] = pd.to_numeric(result_df["dosage"], errors="coerce")
    result_df = result_df.dropna(subset=["dosage"])
    result_df["genotype"] = result_df["dosage"].round().astype(int).clip(0, 2)

    # Human-readable label
    allele2 = variant_id.split(":")[-1] if ":" in variant_id else "ALT"
    allele1 = variant_id.split(":")[-2] if ":" in variant_id else "REF"
    label_map = {
        0: f"{allele1}/{allele1}",
        1: f"{allele1}/{allele2}",
        2: f"{allele2}/{allele2}",
    }
    result_df["genotype_label"] = result_df["genotype"].map(label_map)

    return result_df[["donor_id", "genotype", "genotype_label"]]


def prepare_data(tsv_path: str, peak: str, geno_df: pd.DataFrame,
                 n_bins: int) -> pd.DataFrame:

    needed_cols = ["donor_id", "barcode", "pseudotime", peak]
    print(f"[INFO] Loading TSV: {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t", usecols=needed_cols, low_memory=False)

    if peak not in df.columns:
        raise ValueError(
            f"Peak '{peak}' not found in TSV columns.\n"
            f"Available peaks (first 5): {[c for c in df.columns if c.startswith('chr')][:5]}"
        )

    df = df.dropna(subset=["pseudotime", peak])
    df["pseudotime"] = pd.to_numeric(df["pseudotime"], errors="coerce")
    df[peak] = pd.to_numeric(df[peak], errors="coerce")
    df = df.dropna(subset=["pseudotime", peak])

    # Bin pseudotime
    if n_bins == 3:
        bin_labels = ["Early", "Mid", "Late"]
    else:
        bin_labels = [f"Bin{i + 1}" for i in range(n_bins)]

    pt_min = df["pseudotime"].min()
    pt_max = df["pseudotime"].max()
    bin_edges = np.linspace(pt_min, pt_max, n_bins + 1)
    # Print bin edges so user knows the thresholds used
    print(f"[INFO] Pseudotime range: [{pt_min:.4f}, {pt_max:.4f}]")
    print(f"[INFO] Bin edges: {[f'{e:.4f}' for e in bin_edges]}")

    df["pseudotime_bin"] = pd.cut(
        df["pseudotime"], bins=bin_edges, labels=bin_labels,
        include_lowest=True
    )

    # Merge genotype (donor-level) into cell-level data
    df = df.merge(geno_df, on="donor_id", how="inner")
    print(f"[INFO] Cells after merging with genotype: {len(df)}")
    return df

def inverse_normal_transform(vals):

    from scipy.stats import norm, rankdata
    n    = len(vals)
    rank = rankdata(vals, method="average")
    return norm.ppf((rank - 0.375) / (n + 0.25))


def normalise_and_aggregate(df, peak):

    grp      = ["donor_id", "genotype", "genotype_label", "pseudotime_bin"]
    bin_cats = df["pseudotime_bin"].cat.categories

    agg = df.groupby(grp, observed=True)[peak].mean().reset_index()

    # Step 3: INT within each bin
    for bin_label in bin_cats:
        mask = agg["pseudotime_bin"] == bin_label
        pos_mask = mask & (agg[peak] > 0)
        if pos_mask.sum() > 1:
            agg.loc[pos_mask, peak] = inverse_normal_transform(agg.loc[pos_mask, peak].values)

    agg["pseudotime_bin"] = pd.Categorical(agg["pseudotime_bin"],
                                           categories=bin_cats, ordered=True)
    return {"mean": agg}

def plot(df, peak, variant_id, pval_ge, outdir, return_df=False):
    os.makedirs(outdir, exist_ok=True)
    safe_var = variant_id.replace(":", "_")

    # Aggregate raw counts first, then normalise to log-CPM
    agg_results = normalise_and_aggregate(df, peak)

    if return_df:
        return agg_results

    df_mean = agg_results["mean"]
    fig = _make_plot(df_mean, peak, variant_id, pval_ge,
                     "mean", "Donor mean log-CPM accessibility")
    out = os.path.join(outdir, f"{peak}__{safe_var}__donor_mean.pdf")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out}")

plink = '/g/data/ei56/ax3061/proj/tenk10k/caQTL/data/genotype/TenK10K_TOB_ATAC_renamed_chr16_common_variants_qced'
variant = '16:50336792:G:C'
plink_bin = '/g/data/ei56/jf1058/software/PLINK/plink'
geno_df = extract_genotype_dosage(plink, variant, plink_bin)

tsv = '/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Phenotypes/FinalTSVs_merged/chr16.tsv'
peak = 'chr16_50326863_50327752'
n_pseudotime_bins = 3
df = prepare_data(tsv, peak, geno_df, n_pseudotime_bins)

outdir = '/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/plots/'
df = plot(df, peak, variant, 6.76e-28, outdir, return_df=True)
df['mean'].to_csv('/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Plotting/plots/' + peak + '_' + variant.replace(':', '_') + '_df_' + str(n_pseudotime_bins) + 'bin.csv', index=False)
