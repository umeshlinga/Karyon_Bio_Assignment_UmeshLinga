#!/usr/bin/env python3
"""
01_download_and_qc.py
Per-sample loading and quality control for GSE136103 (Ramachandran 2019).

What this script does
---------------------
1. Discovers per-sample 10x folders under <raw_dir> (each folder must contain
   matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz). Sample folders are named
   "GSM<digits>_<donor>_<fraction>"; see utils.parse_sample_dir for parsing.
2. Loads each sample into AnnData, annotates donor/disease/fraction.
3. Computes QC metrics: n_genes_by_counts, total_counts, pct_counts_mt,
   pct_counts_ribo, pct_counts_hb.
4. Applies *liver-fibrosis-aware* QC thresholds:
       - min_genes >= 200
       - min_counts >= 500
       - pct_counts_mt <= 20%        (relaxed for stressed/diseased cells)
       - pct_counts_hb <= 50%        (drop heavy blood-only contamination)
5. Detects doublets per sample via scanpy.pp.scrublet (fast, no extra dep).
6. Writes:
       results/qc/qc_per_sample.csv       (pre/post counts, retention)
       results/qc/qc_violin.png           (per-sample violin plots)
       results/qc/qc_summary.md           (human-readable QC report)
       results/processed/qc_filtered.h5ad (concatenated filtered AnnData)

Usage
-----
python scripts/01_download_and_qc.py \
    --raw_dir /path/to/extracted_GSE136103/per_sample \
    --out_dir results \
    [--max_cells_per_sample 800]   # optional subsampling for fast demos
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils as U                                       # noqa: E402

sc.settings.verbosity = 1
sc.settings.figdir = "."
sc.settings.set_figure_params(dpi=120, facecolor="white")


MITO_MAX = 20.0      # %; relaxed for fibrotic liver
HB_MAX = 50.0        # %; allow some blood contam, drop only extreme
RIBO_MAX = 60.0      # %; drop only extreme ribo contamination
MIN_GENES = 200
MIN_COUNTS = 500


def add_qc_metrics(ad_obj):
    ad_obj.var["mt"] = ad_obj.var_names.str.startswith(("MT-", "MT."))
    ad_obj.var["ribo"] = ad_obj.var_names.str.startswith(("RPS", "RPL"))
    ad_obj.var["hb"] = ad_obj.var_names.str.match(r"^HB[ABDEGMQZ]\d?$")
    sc.pp.calculate_qc_metrics(
        ad_obj,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    return ad_obj


def run_doublet_detection(ad_obj, sample_key="sample_id"):
    """Per-sample doublet scoring with scrublet (built into scanpy)."""
    try:
        sc.pp.scrublet(ad_obj, batch_key=sample_key, verbose=False)
    except Exception as exc:                              # pragma: no cover
        print(f"[warn] scrublet failed ({exc}); marking doublet column False.")
        ad_obj.obs["predicted_doublet"] = False
        ad_obj.obs["doublet_score"] = 0.0
    return ad_obj


def apply_filters(ad_obj):
    pre = ad_obj.n_obs
    keep = (
        (ad_obj.obs["n_genes_by_counts"] >= MIN_GENES)
        & (ad_obj.obs["total_counts"] >= MIN_COUNTS)
        & (ad_obj.obs["pct_counts_mt"] <= MITO_MAX)
        & (ad_obj.obs["pct_counts_hb"] <= HB_MAX)
        & (ad_obj.obs["pct_counts_ribo"] <= RIBO_MAX)
        & (~ad_obj.obs["predicted_doublet"].astype(bool))
    )
    ad_obj = ad_obj[keep].copy()
    sc.pp.filter_genes(ad_obj, min_cells=3)
    return ad_obj, pre, int(keep.sum())


def make_qc_violin(ad_obj, out_path: Path) -> None:
    obs = ad_obj.obs
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, col, title in zip(
        axes,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        ["# genes / cell", "total UMI / cell", "% mito / cell"],
    ):
        groups = sorted(obs["sample_id"].unique())
        data = [obs.loc[obs["sample_id"] == g, col].values for g in groups]
        parts = ax.violinplot(data, showmeans=False, showmedians=True)
        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels([g.split("_", 1)[-1] for g in groups],
                           rotation=90, fontsize=6)
        ax.set_title(title)
        ax.set_yscale("log" if "counts" in col else "linear")
    plt.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw_dir", required=True, type=Path,
                    help="Directory of per-sample 10x folders.")
    ap.add_argument("--out_dir", default=Path("results"), type=Path)
    ap.add_argument("--max_cells_per_sample", type=int, default=None,
                    help="Optional cap to speed up demo runs.")
    ap.add_argument("--include_pbmc", action="store_true",
                    help="Include blood / PBMC samples (default: liver only).")
    args = ap.parse_args()

    qc_dir = args.out_dir / "qc"
    proc_dir = args.out_dir / "processed"
    qc_dir.mkdir(parents=True, exist_ok=True)
    proc_dir.mkdir(parents=True, exist_ok=True)

    samples = U.discover_samples(args.raw_dir, include_pbmc=args.include_pbmc)
    if not samples:
        raise SystemExit(f"No GSE136103-style samples found under {args.raw_dir}")
    sheet = U.build_sample_sheet(samples)
    sheet.to_csv(qc_dir / "sample_sheet.csv", index=False)
    print(f"Discovered {len(samples)} samples "
          f"({(sheet['disease'] == 'cirrhotic').sum()} cirrhotic, "
          f"{(sheet['disease'] == 'healthy').sum()} healthy)")

    adata = U.load_concat(samples, max_cells_per_sample=args.max_cells_per_sample)
    print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    adata = add_qc_metrics(adata)

    # Per-sample summary BEFORE filtering
    pre_summary = (
        adata.obs.groupby("sample_id", observed=True)
        .agg(n_cells_pre=("sample_id", "size"),
             med_genes=("n_genes_by_counts", "median"),
             med_counts=("total_counts", "median"),
             med_mt=("pct_counts_mt", "median"))
    )

    adata = run_doublet_detection(adata)
    pre_doublets = int(adata.obs["predicted_doublet"].sum())
    print(f"Predicted doublets (Scrublet): {pre_doublets:,}")

    adata, n_pre, n_post = apply_filters(adata)
    print(f"Filtered: {n_pre:,} -> {n_post:,} cells "
          f"({100.0 * n_post / max(n_pre, 1):.1f}% retained)")

    post_summary = (
        adata.obs.groupby("sample_id", observed=True)
        .agg(n_cells_post=("sample_id", "size"))
    )
    summary = pre_summary.join(post_summary, how="left").fillna(0)
    summary["retention_pct"] = (
        100.0 * summary["n_cells_post"] / summary["n_cells_pre"]
    )
    summary.to_csv(qc_dir / "qc_per_sample.csv")

    make_qc_violin(adata, qc_dir / "qc_violin.png")

    # Save raw counts as a layer before any normalization downstream.
    adata.layers["counts"] = adata.X.copy()
    adata.write_h5ad(proc_dir / "qc_filtered.h5ad")

    md = qc_dir / "qc_summary.md"
    with md.open("w", encoding="utf-8") as fh:
        fh.write("# GSE136103 — QC Summary\n\n")
        fh.write(f"- Samples discovered: **{len(samples)}**\n")
        fh.write(f"- Cells loaded (post-subsample): **{n_pre:,}**\n")
        fh.write(f"- Cells retained after QC: **{n_post:,}** "
                 f"({100.0 * n_post / max(n_pre, 1):.1f}%)\n")
        fh.write(f"- Predicted doublets (Scrublet, pre-filter): **{pre_doublets:,}**\n")
        fh.write("\n## QC thresholds (liver-fibrosis-aware)\n")
        fh.write(f"- min_genes >= {MIN_GENES}; min_counts >= {MIN_COUNTS}\n")
        fh.write(f"- pct_counts_mt <= {MITO_MAX}%  *(relaxed: activated HSCs and "
                 "SAMs in fibrotic livers legitimately have elevated mito)*\n")
        fh.write(f"- pct_counts_hb <= {HB_MAX}%; pct_counts_ribo <= {RIBO_MAX}%\n")
        fh.write("\n## Per-sample retention\n\n")
        fh.write(summary.round(1).to_markdown())
        fh.write("\n")
    print(f"Wrote QC report -> {md}")
    print(f"Wrote filtered AnnData -> {proc_dir / 'qc_filtered.h5ad'}")


if __name__ == "__main__":
    main()
