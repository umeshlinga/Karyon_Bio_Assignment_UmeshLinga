#!/usr/bin/env python3
"""
02_preprocess_integrate.py
Normalization, HVG selection, dimensionality reduction, batch integration
(Harmony if available, otherwise BBKNN fallback), neighbors, UMAP, Leiden.

Design notes
------------
- Donor is used as the integration covariate. Disease status is INTENTIONALLY
  preserved — never put it into the batch key, which would erase fibrosis biology.
- Validation: after integration we check that disease-status separation is
  preserved within mesenchymal/macrophage/endothelial compartments by inspecting
  silhouette w.r.t. cell-type proxy markers (deferred to script 03).
- For commodity hardware we use Harmony via harmonypy if installed; otherwise
  fall back to BBKNN. Both preserve per-cell gene expression while correcting
  graph/embedding structure across donors.

Inputs
------
results/processed/qc_filtered.h5ad   (from 01_download_and_qc.py)

Outputs
-------
results/processed/integrated.h5ad    (with X_pca, X_umap, clusters)
results/figures/umap_donor.png
results/figures/umap_disease.png
results/figures/umap_cluster.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils as U                                       # noqa: E402

sc.settings.verbosity = 1


def try_harmony(adata, batch_key):
    try:
        import harmonypy as hm
    except ImportError:
        return False
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_key)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != adata.n_obs:
        Z = Z.T
    if Z.shape != (adata.n_obs, adata.obsm["X_pca"].shape[1]):
        raise RuntimeError(
            f"Harmony output shape {Z.shape} incompatible with X_pca "
            f"{adata.obsm['X_pca'].shape}; expected (n_obs, n_pcs).")
    adata.obsm["X_pca_harmony"] = Z
    return True


def try_bbknn(adata, batch_key):
    try:
        import bbknn                                    # noqa: F401
    except ImportError:
        return False
    sc.external.pp.bbknn(adata, batch_key=batch_key,
                        neighbors_within_batch=3)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_h5ad", required=True, type=Path)
    ap.add_argument("--out_dir", default=Path("results"), type=Path)
    ap.add_argument("--n_hvg", type=int, default=3000)
    ap.add_argument("--n_pcs", type=int, default=40)
    ap.add_argument("--leiden_resolution", type=float, default=0.8)
    ap.add_argument("--batch_key", default="donor")
    args = ap.parse_args()

    proc_dir = args.out_dir / "processed"
    fig_dir = args.out_dir / "figures"
    proc_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.in_h5ad)
    print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # 1. Normalize + log
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 2. HVGs (donor-aware: pick genes variable in many donors)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg,
                                batch_key=args.batch_key, flavor="seurat")
    n_hvg_used = int(adata.var["highly_variable"].sum())
    print(f"Selected {n_hvg_used} HVGs (donor-aware)")

    # 3. Scale + PCA on HVGs
    adata.raw = adata
    adata_h = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_h, max_value=10)
    sc.tl.pca(adata_h, n_comps=args.n_pcs, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_h.obsm["X_pca"]

    # 4. Integration
    used = "none"
    if try_harmony(adata, args.batch_key):
        used = "harmony"
        sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
    elif try_bbknn(adata, args.batch_key):
        used = "bbknn"
    else:
        print("[warn] No batch-correction backend; falling back to plain PCA.")
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
    print(f"Integration backend: {used}")

    # 5. UMAP + Leiden
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=args.leiden_resolution, key_added="leiden")
    print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

    # 6. UMAP figures
    for color, fname in [("donor",   "umap_donor.png"),
                         ("disease", "umap_disease.png"),
                         ("leiden",  "umap_cluster.png")]:
        ax = sc.pl.umap(adata, color=color, show=False, legend_fontsize=6,
                        size=8, frameon=False)
        fig = ax.figure if hasattr(ax, "figure") else plt.gcf()
        fig.set_size_inches(7, 6)
        fig.savefig(fig_dir / fname, dpi=140, bbox_inches="tight")
        plt.close(fig)

    out_path = proc_dir / "integrated.h5ad"
    adata.write_h5ad(out_path)
    print(f"Wrote integrated AnnData -> {out_path}")


if __name__ == "__main__":
    main()
