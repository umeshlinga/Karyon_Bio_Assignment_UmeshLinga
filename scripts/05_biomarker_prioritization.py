#!/usr/bin/env python3
"""
05_biomarker_prioritization.py
Hybrid rule-based + ML biomarker / target prioritization for liver fibrosis.

Inputs
------
results/de/pseudobulk_all_compartments.csv     (from 04_de_pathway_cci.py)
results/processed/annotated.h5ad               (for cell-type-specific expression)
results/de/cci_evidence_panel.csv              (optional CCI score)

What it does
------------
1. Engineers a per-gene, per-compartment feature matrix:
       - signed log2FC (cirrhotic vs healthy)         |log2FC|
       - -log10(padj)
       - compartment specificity tau (across cell types)
       - detection rate in cirrhotic compartment
       - druggability flag (DGIdb/Open Targets curated subset)
       - surface_or_secreted flag (translational tractability)
       - pathway-membership boolean (fibrosis signature gene)
       - literature-fibrosis flag (curated)
       - CCI evidence (max interaction_score where this gene is ligand OR receptor)

2. Composite rule-based prioritization score
       0.20 * specificity_norm
     + 0.20 * abs_log2FC_norm
     + 0.10 * negLogP_norm
     + 0.10 * detection_rate_cirrh
     + 0.10 * pathway_relevance
     + 0.15 * druggability
     + 0.10 * literature_flag
     + 0.05 * cci_evidence_norm

3. ML cross-check: train a RandomForestRegressor on the engineered features
   using literature_fibrosis as a noisy positive proxy. We treat it as a
   sanity-check estimator that gives feature importances + a model-ranked
   score; the final ranking blends rule-based and model-based scores 70/30.

4. Outputs
       results/biomarker_ranked_table.csv          (Top N ranked candidates)
       results/biomarker_feature_matrix.csv        (full feature matrix)
       results/figures/feature_importance.png      (RF importances)
       results/figures/top_biomarkers.png          (bar plot of top 15)
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
from scipy import sparse
from sklearn.ensemble import RandomForestRegressor

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils as U                                       # noqa: E402


PATHWAY_GENES = set().union(*U.FIBROSIS_SIGNATURES.values())


def build_compartment_expression(adata) -> pd.DataFrame:
    """Mean log-normalized expression per (gene x compartment)."""
    if adata.raw is not None:
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names
    if sparse.issparse(X):
        X = X.tocsc()
    out = {}
    for comp in adata.obs["compartment"].cat.categories:
        idx = np.where(adata.obs["compartment"].values == comp)[0]
        if len(idx) == 0:
            continue
        sub = X[idx]
        m = np.asarray(sub.mean(axis=0)).ravel()
        out[comp] = m
    return pd.DataFrame(out, index=var_names)


def build_compartment_detection_rate(adata) -> pd.DataFrame:
    """Fraction of cells expressing each gene per (gene x compartment), cirrhotic only."""
    if adata.raw is not None:
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names
    if sparse.issparse(X):
        X = X.tocsc()
    obs = adata.obs
    out = {}
    for comp in obs["compartment"].cat.categories:
        idx = np.where((obs["compartment"].values == comp)
                       & (obs["disease"].values == "cirrhotic"))[0]
        if len(idx) == 0:
            out[comp] = np.zeros(len(var_names))
            continue
        sub = X[idx]
        if sparse.issparse(sub):
            det = np.asarray((sub > 0).sum(axis=0)).ravel() / sub.shape[0]
        else:
            det = (sub > 0).mean(axis=0)
        out[comp] = det
    return pd.DataFrame(out, index=var_names)


def assemble_features(de_df: pd.DataFrame,
                      expr_df: pd.DataFrame,
                      det_df: pd.DataFrame,
                      cci_df) -> pd.DataFrame:
    rows = []
    tau = U.specificity_tau(expr_df)

    cci_score_by_gene = {}
    if cci_df is not None and not cci_df.empty:
        for col in ("ligand", "receptor"):
            g = cci_df.groupby(col)["interaction_score"].max()
            for k, v in g.items():
                cci_score_by_gene[k] = max(cci_score_by_gene.get(k, 0.0), float(v))

    for _, row in de_df.iterrows():
        gene = row["gene"]
        comp = row["compartment"]
        lfc = float(row.get("log2FoldChange", np.nan))
        padj = float(row.get("padj", np.nan))
        if not np.isfinite(lfc) or not np.isfinite(padj):
            continue
        if padj <= 0:
            padj = 1e-300
        detection = (float(det_df.loc[gene, comp])
                     if gene in det_df.index and comp in det_df.columns else 0.0)
        rows.append({
            "gene": gene,
            "compartment": comp,
            "log2FoldChange": lfc,
            "abs_log2FC": abs(lfc),
            "padj": padj,
            "negLogP": -np.log10(padj),
            "tau_specificity": float(tau.get(gene, 0.0)),
            "detection_rate_cirrh": detection,
            "pathway_relevance": int(gene in PATHWAY_GENES),
            "druggability": int(gene in U.DRUGGABLE_GENES),
            "surface_or_secreted": int(gene in U.SURFACE_OR_SECRETED),
            "literature_fibrosis": int(gene in U.LITERATURE_FIBROSIS),
            "cci_evidence": float(cci_score_by_gene.get(gene, 0.0)),
        })
    return pd.DataFrame(rows)


def rule_based_score(feat: pd.DataFrame) -> pd.Series:
    s = (
        0.20 * U.minmax(feat["tau_specificity"].values)
        + 0.20 * U.minmax(feat["abs_log2FC"].values)
        + 0.10 * U.minmax(feat["negLogP"].values)
        + 0.10 * feat["detection_rate_cirrh"].clip(0, 1).values
        + 0.10 * feat["pathway_relevance"].values
        + 0.15 * feat["druggability"].values
        + 0.10 * feat["literature_fibrosis"].values
        + 0.05 * U.minmax(feat["cci_evidence"].values)
    )
    return pd.Series(s, index=feat.index, name="rule_score")


def ml_score(feat: pd.DataFrame, out_dir: Path) -> pd.Series:
    """RandomForest regressor with literature_fibrosis as a noisy target."""
    X_cols = [
        "abs_log2FC", "negLogP", "tau_specificity", "detection_rate_cirrh",
        "pathway_relevance", "druggability", "surface_or_secreted",
        "cci_evidence",
    ]
    X = feat[X_cols].values.astype(float)
    y = feat["literature_fibrosis"].values.astype(float)
    if y.sum() < 5:
        return pd.Series(np.zeros(len(feat)), index=feat.index, name="ml_score")
    rf = RandomForestRegressor(n_estimators=300, max_depth=None,
                               min_samples_leaf=3, random_state=0, n_jobs=-1)
    rf.fit(X, y)
    pred = rf.predict(X)
    imp = pd.Series(rf.feature_importances_, index=X_cols).sort_values()
    fig, ax = plt.subplots(figsize=(6, 4))
    imp.plot.barh(ax=ax, color="steelblue")
    ax.set_xlabel("Random Forest feature importance")
    ax.set_title("Prioritization model — feature importances")
    fig.tight_layout()
    fig.savefig(out_dir / "feature_importance.png", dpi=140,
                bbox_inches="tight")
    plt.close(fig)
    return pd.Series(U.minmax(pred), index=feat.index, name="ml_score")


def translational_category(row) -> str:
    if row["surface_or_secreted"] and row["druggability"]:
        return "Therapeutic (druggable, tractable)"
    if row["surface_or_secreted"]:
        return "Diagnostic / surface biomarker"
    if row["pathway_relevance"] and row["druggability"]:
        return "Therapeutic (pathway target)"
    if row["pathway_relevance"]:
        return "Mechanistic biomarker"
    return "Exploratory"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de_csv", required=True, type=Path)
    ap.add_argument("--in_h5ad", required=True, type=Path)
    ap.add_argument("--cci_csv", type=Path, default=None)
    ap.add_argument("--out_dir", default=Path("results"), type=Path)
    ap.add_argument("--top_k", type=int, default=20)
    args = ap.parse_args()

    fig_dir = args.out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    de_df = pd.read_csv(args.de_csv)
    de_df = de_df[(de_df["padj"] < 0.05) & (de_df["log2FoldChange"].abs() > 0.5)]
    print(f"DE input rows (filtered): {len(de_df):,}; compartments: "
          f"{sorted(de_df['compartment'].unique())}")

    cci_df = (pd.read_csv(args.cci_csv)
              if args.cci_csv and args.cci_csv.exists() else None)

    adata = sc.read_h5ad(args.in_h5ad)
    expr_df = build_compartment_expression(adata)
    det_df = build_compartment_detection_rate(adata)

    feat = assemble_features(de_df, expr_df, det_df, cci_df)
    if feat.empty:
        raise SystemExit("No features assembled (no significant DE rows).")
    print(f"Feature matrix: {len(feat):,} (gene x compartment) rows")

    feat["rule_score"] = rule_based_score(feat)
    feat["ml_score"] = ml_score(feat, fig_dir)
    feat["final_score"] = 0.7 * feat["rule_score"] + 0.3 * feat["ml_score"]
    feat["translational_category"] = feat.apply(translational_category, axis=1)

    feat = feat.sort_values("final_score", ascending=False)
    feat_unique = feat.drop_duplicates(subset=["gene"], keep="first")

    feat.to_csv(args.out_dir / "biomarker_feature_matrix.csv", index=False)
    top = feat_unique.head(args.top_k).copy()
    top.insert(0, "Rank", np.arange(1, len(top) + 1))
    top = top[[
        "Rank", "gene", "compartment", "translational_category",
        "log2FoldChange", "padj", "tau_specificity", "detection_rate_cirrh",
        "pathway_relevance", "druggability", "surface_or_secreted",
        "literature_fibrosis", "cci_evidence",
        "rule_score", "ml_score", "final_score",
    ]]
    top.columns = [
        "Rank", "Gene", "Primary_Compartment", "Translational_Category",
        "log2FC_cirrh_vs_healthy", "padj", "Tau_Specificity",
        "Detection_Rate_Cirrh", "Fibrosis_Pathway_Member", "Druggable_Flag",
        "Surface_or_Secreted", "Literature_Fibrosis", "CCI_Evidence",
        "Rule_Score", "ML_Score", "Final_Prioritization_Score",
    ]
    out_csv = args.out_dir / "biomarker_ranked_table.csv"
    top.round(4).to_csv(out_csv, index=False)
    print(f"\nWrote top-{args.top_k} ranked biomarker table -> {out_csv}")
    print(top.head(15).to_string(index=False))

    show = top.head(15)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    colors = {"Mesenchyme": "#a06b3c", "Macrophage": "#3c7aa0",
              "Endothelial": "#c0392b"}
    bar_colors = [colors.get(c, "#888") for c in show["Primary_Compartment"]]
    ax.barh(show["Gene"][::-1], show["Final_Prioritization_Score"][::-1],
            color=bar_colors[::-1])
    ax.set_xlabel("Final prioritization score")
    ax.set_title("Top 15 prioritized fibrosis biomarker / target candidates")
    fig.tight_layout()
    fig.savefig(fig_dir / "top_biomarkers.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
