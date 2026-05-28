#!/usr/bin/env python3
"""
03_annotate_validate.py
Marker-gene based cell-type annotation + validation of the three disease-
relevant compartments required by the Karyon assignment:
  - Hepatic stellate / mesenchymal / myofibroblast-like cells
  - Macrophage / monocyte populations
  - Endothelial cells

Approach
--------
1. Score each cluster against canonical liver marker panels (utils.LIVER_MARKERS)
   using sc.tl.score_genes. The cluster is labeled by argmax of mean cluster
   score across panels, with a small penalty when multiple panels score high
   (mixed clusters get the "Mixed" suffix).
2. Compute fibrosis-signature scores per cell (SAM, Scar-EC, HSC-activation,
   ECM, TGFB) and verify that they are higher in cirrhotic vs healthy donors
   within each compartment (Mann-Whitney U on donor-level means).
3. Save UMAPs colored by cell_type + dotplot of canonical markers.

Inputs
------
results/processed/integrated.h5ad

Outputs
-------
results/processed/annotated.h5ad
results/figures/umap_celltype.png
results/figures/dotplot_markers.png
results/figures/violin_fibrosis_scores.png
results/annotation/cluster_marker_scores.csv
results/annotation/compartment_validation.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils as U                                       # noqa: E402

sc.settings.verbosity = 1


# Map marker-panel name -> coarse cell_type label written into obs
PANEL_TO_TYPE = {
    "HSC_quiescent":      "HSC_quiescent",
    "HSC_activated":      "HSC_activated/Myofibroblast",
    "Kupffer":            "Kupffer",
    "Monocyte":           "Monocyte",
    "SAM":                "Scar-Associated Macrophage",
    "LSEC":               "LSEC",
    "Endothelial_vasc":   "Endothelial_vasc",
    "Scar_endothelial":   "Scar-Endothelial",
    "Cholangiocyte":      "Cholangiocyte",
    "Hepatocyte":         "Hepatocyte",
    "T_cell":             "T_cell",
    "B_cell":             "B_cell",
    "Plasma":             "Plasma",
    "NK":                 "NK",
    "pDC":                "pDC",
    "Mast":               "Mast",
}

# Coarse compartment label used for downstream DE and prioritization.
TYPE_TO_COMPARTMENT = {
    "HSC_quiescent":              "Mesenchyme",
    "HSC_activated/Myofibroblast":"Mesenchyme",
    "Kupffer":                    "Macrophage",
    "Monocyte":                   "Macrophage",
    "Scar-Associated Macrophage": "Macrophage",
    "LSEC":                       "Endothelial",
    "Endothelial_vasc":           "Endothelial",
    "Scar-Endothelial":           "Endothelial",
    "Cholangiocyte":              "Cholangiocyte",
    "Hepatocyte":                 "Hepatocyte",
    "T_cell":                     "T/NK",
    "B_cell":                     "B/Plasma",
    "Plasma":                     "B/Plasma",
    "NK":                         "T/NK",
    "pDC":                        "DC",
    "Mast":                       "Mast",
}


def score_panels(adata):
    """Per-cell score of each marker panel (mean z-scored expression of panel)."""
    for name, genes in U.LIVER_MARKERS.items():
        present = [g for g in genes if g in adata.raw.var_names]
        if not present:
            adata.obs[f"score__{name}"] = 0.0
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score__{name}",
                          use_raw=True, random_state=0)
    return adata


def score_fibrosis_signatures(adata):
    for name, genes in U.FIBROSIS_SIGNATURES.items():
        present = [g for g in genes if g in adata.raw.var_names]
        if not present:
            adata.obs[f"sig__{name}"] = 0.0
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=f"sig__{name}",
                          use_raw=True, random_state=0)
    return adata


def annotate_clusters(adata, cluster_key="leiden"):
    score_cols = [c for c in adata.obs.columns if c.startswith("score__")]
    cluster_means = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
    cluster_means.columns = [c.replace("score__", "") for c in cluster_means.columns]

    # Assign each cluster the panel with the highest mean score.
    # If the top panel is < 0 OR within 0.1 of the second panel, flag as Mixed.
    labels: dict = {}
    rationale_rows = []
    for cl, row in cluster_means.iterrows():
        sorted_panels = row.sort_values(ascending=False)
        top_panel = sorted_panels.index[0]
        top_score = sorted_panels.iloc[0]
        gap = top_score - sorted_panels.iloc[1] if len(sorted_panels) > 1 else 1.0
        cell_type = PANEL_TO_TYPE.get(top_panel, top_panel)
        if top_score < 0.05:
            cell_type = "Unknown"
        elif gap < 0.05:
            cell_type = f"{cell_type} (mixed/{PANEL_TO_TYPE.get(sorted_panels.index[1], sorted_panels.index[1])})"
        labels[cl] = cell_type
        rationale_rows.append({
            "cluster": cl,
            "cell_type": cell_type,
            "top_panel": top_panel,
            "top_score": float(top_score),
            "gap_to_second": float(gap),
        })

    adata.obs["cell_type"] = adata.obs[cluster_key].map(labels).astype("category")
    compartment_map = {ct: TYPE_TO_COMPARTMENT.get(ct.split(" (mixed")[0], "Other")
                       for ct in adata.obs["cell_type"].cat.categories}
    adata.obs["compartment"] = (
        adata.obs["cell_type"].map(compartment_map).astype("category")
    )
    return adata, cluster_means, pd.DataFrame(rationale_rows)


def validate_compartments(adata, sig_key_prefix="sig__"):
    """Donor-level Mann-Whitney U of fibrosis signature scores: cirrhotic vs healthy
    within each compartment of interest."""
    rows = []
    sig_cols = [c for c in adata.obs.columns if c.startswith(sig_key_prefix)]
    for comp in ["Mesenchyme", "Macrophage", "Endothelial"]:
        sub = adata.obs[adata.obs["compartment"] == comp]
        if sub.empty:
            continue
        for sig in sig_cols:
            donor_means = (sub.groupby(["donor", "disease"], observed=True)[sig]
                              .mean().reset_index())
            healthy = donor_means.loc[donor_means["disease"] == "healthy", sig]
            cirrh = donor_means.loc[donor_means["disease"] == "cirrhotic", sig]
            if len(healthy) < 2 or len(cirrh) < 2:
                u, p = float("nan"), float("nan")
            else:
                u, p = stats.mannwhitneyu(cirrh, healthy, alternative="greater")
            rows.append({
                "compartment": comp,
                "signature": sig.replace(sig_key_prefix, ""),
                "n_donor_healthy": int(len(healthy)),
                "n_donor_cirrh": int(len(cirrh)),
                "mean_healthy": float(healthy.mean()) if len(healthy) else float("nan"),
                "mean_cirrh": float(cirrh.mean()) if len(cirrh) else float("nan"),
                "delta": float(cirrh.mean() - healthy.mean())
                         if len(healthy) and len(cirrh) else float("nan"),
                "mwu_p_one_sided": float(p),
            })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_h5ad", required=True, type=Path)
    ap.add_argument("--out_dir", default=Path("results"), type=Path)
    args = ap.parse_args()

    fig_dir = args.out_dir / "figures"
    proc_dir = args.out_dir / "processed"
    ann_dir = args.out_dir / "annotation"
    for d in (fig_dir, proc_dir, ann_dir):
        d.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.in_h5ad)
    if adata.raw is None:
        adata.raw = adata
    print(f"Loaded {adata.n_obs:,} cells from {args.in_h5ad}")

    adata = score_panels(adata)
    adata = score_fibrosis_signatures(adata)

    adata, cluster_means, rationale = annotate_clusters(adata)
    cluster_means.to_csv(ann_dir / "cluster_marker_scores.csv")
    rationale.to_csv(ann_dir / "cluster_label_rationale.csv", index=False)
    print("Cell-type labels:")
    print(adata.obs["cell_type"].value_counts().to_string())

    validation = validate_compartments(adata)
    validation.to_csv(ann_dir / "compartment_validation.csv", index=False)
    print("\nCompartment validation (signature in cirrh vs healthy):")
    print(validation.round(4).to_string(index=False))

    # UMAP by cell_type
    ax = sc.pl.umap(adata, color="cell_type", show=False, legend_fontsize=6,
                    size=8, frameon=False)
    fig = ax.figure if hasattr(ax, "figure") else plt.gcf()
    fig.set_size_inches(8.5, 6.5)
    fig.savefig(fig_dir / "umap_celltype.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # Dotplot of canonical markers per cell_type
    marker_dot = {k: [g for g in v if g in adata.raw.var_names][:4]
                  for k, v in U.LIVER_MARKERS.items()}
    marker_dot = {k: v for k, v in marker_dot.items() if v}
    try:
        sc.pl.dotplot(adata, marker_dot, groupby="cell_type",
                      standard_scale="var", show=False,
                      use_raw=True, figsize=(10, 6))
        plt.savefig(fig_dir / "dotplot_markers.png", dpi=140, bbox_inches="tight")
        plt.close()
    except Exception as exc:
        print(f"[warn] dotplot failed ({exc}); skipping.")

    # Violin of fibrosis signatures by disease, faceted on compartment
    sig_cols = ["sig__SAM_signature", "sig__ScarEC_signature",
                "sig__HSC_activation", "sig__ECM_organization"]
    present = [c for c in sig_cols if c in adata.obs.columns]
    if present:
        fig, axes = plt.subplots(1, len(present), figsize=(4 * len(present), 4))
        if len(present) == 1:
            axes = [axes]
        for ax, col in zip(axes, present):
            data_by_disease = []
            labels = []
            for d in ("healthy", "cirrhotic"):
                vals = adata.obs.loc[adata.obs["disease"] == d, col].values
                data_by_disease.append(vals)
                labels.append(d)
            ax.violinplot(data_by_disease, showmedians=True)
            ax.set_xticks([1, 2]); ax.set_xticklabels(labels)
            ax.set_title(col.replace("sig__", ""))
        plt.tight_layout()
        plt.savefig(fig_dir / "violin_fibrosis_scores.png", dpi=140,
                    bbox_inches="tight")
        plt.close()

    adata.write_h5ad(proc_dir / "annotated.h5ad")
    print(f"Wrote annotated AnnData -> {proc_dir / 'annotated.h5ad'}")


if __name__ == "__main__":
    main()
