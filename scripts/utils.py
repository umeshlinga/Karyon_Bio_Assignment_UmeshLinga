#!/usr/bin/env python3
"""
utils.py
Shared helpers for the Karyon Bio liver-fibrosis biomarker discovery pipeline.

Provides:
    - GSE136103 sample-sheet parsing from the canonical filename scheme
      "GSM<id>_<donor>_<fraction>" (e.g. healthy1_cd45+, cirrhotic2_cd45-A)
    - Per-sample 10x .mtx loader that concatenates into a single AnnData
    - Marker-gene dictionaries grounded in Ramachandran et al. 2019
    - Pseudobulk aggregation (donor x cell_type) for donor-aware DE
    - Min-max + safe-divide helpers used by the prioritization scorer
"""

from __future__ import annotations

import gzip
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


# ---------------------------------------------------------------------------
# Marker panels (Ramachandran 2019 Nature + standard liver atlases)
# ---------------------------------------------------------------------------

# Major NPC compartments expected in GSE136103
LIVER_MARKERS: Dict[str, List[str]] = {
    # Mesenchymal / stromal lineage
    "HSC_quiescent":      ["LRAT", "RBP1", "REEP6", "CYGB", "PPARG", "GFAP"],
    "HSC_activated":      ["ACTA2", "COL1A1", "COL3A1", "TAGLN", "PDGFRB",
                           "LUM", "DCN", "LOX", "LOXL2", "TIMP1", "POSTN",
                           "SPARC", "CTHRC1", "COL1A2"],
    # Macrophage lineage
    "Kupffer":            ["CD163", "VSIG4", "MARCO", "TIMD4", "MAF",
                           "MRC1", "LYVE1", "CD5L"],
    "Monocyte":           ["VCAN", "FCN1", "S100A8", "S100A9", "LYZ", "CD14"],
    "SAM":                ["TREM2", "CD9", "GPNMB", "SPP1", "CCR2", "FABP5",
                           "MMP9"],
    # Endothelial lineage
    "LSEC":               ["STAB2", "STAB1", "FCN3", "CLEC4G", "CLEC1B",
                           "OIT3"],
    "Endothelial_vasc":   ["VWF", "RAMP3", "PECAM1", "CDH5"],
    "Scar_endothelial":   ["ACKR1", "PLVAP", "RGCC"],
    # Epithelial
    "Cholangiocyte":      ["KRT19", "KRT7", "EPCAM", "SOX9", "CFTR"],
    "Hepatocyte":         ["ALB", "TF", "APOB", "HNF4A", "CYP3A4"],
    # Immune
    "T_cell":             ["CD3D", "CD3E", "CD8A", "CD4", "TRAC"],
    "B_cell":             ["CD79A", "MS4A1", "CD19"],
    "Plasma":             ["IGKC", "MZB1", "JCHAIN"],
    "NK":                 ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "pDC":                ["LILRA4", "IRF7", "CLEC4C"],
    "Mast":               ["TPSAB1", "TPSB2", "CPA3"],
}

# Disease-state signatures used for AddModuleScore-style scoring
FIBROSIS_SIGNATURES: Dict[str, List[str]] = {
    "SAM_signature":          ["TREM2", "CD9", "GPNMB", "SPP1", "FABP5", "MMP9"],
    "ScarEC_signature":       ["ACKR1", "PLVAP", "RGCC"],
    "HSC_activation":         ["ACTA2", "COL1A1", "COL3A1", "PDGFRB",
                               "LOX", "LOXL2", "TIMP1", "CTHRC1"],
    "ECM_organization":       ["COL1A1", "COL1A2", "COL3A1", "FN1", "LOX",
                               "LOXL1", "LOXL2", "DCN", "BGN"],
    "TGF_beta":               ["TGFB1", "TGFBR1", "TGFBR2", "SMAD3", "SMAD7",
                               "SERPINE1"],
}

# Druggable / tractable gene flags. Coarse first-pass mapping based on
# DGIdb + Open Targets snapshots (refresh in production via API).
DRUGGABLE_GENES = {
    # Receptors with clinical-stage inhibitors
    "PDGFRB", "PDGFRA", "TGFBR1", "TGFBR2", "NOTCH1", "NOTCH2", "NOTCH3",
    "ACKR1", "CCR2", "CXCR4", "ITGAV", "ITGB1", "ITGB6", "EGFR", "FGFR1",
    "MMP9", "MMP2", "LOXL2", "TIMP1", "TREM2", "CD9", "CSF1R",
    # Secreted ligands
    "TGFB1", "PDGFB", "PDGFA", "CTGF", "WNT5A", "IL6", "TNF", "TNFSF12",
    "TNFRSF12A", "CCL2", "CXCL12", "SPP1", "POSTN", "SERPINE1",
}

SURFACE_OR_SECRETED = {
    # Surface receptors and ligands; potential antibody / soluble biomarker
    "TREM2", "CD9", "PDGFRB", "PDGFRA", "ACKR1", "PLVAP", "CCR2", "CSF1R",
    "TGFBR1", "TGFBR2", "ITGA11", "TNFRSF12A", "NOTCH1", "NOTCH2", "NOTCH3",
    "VCAM1", "ICAM1", "CDH5", "PECAM1", "EPCAM", "KRT19",
    # Soluble / matricellular candidates suitable for serum readouts
    "SPP1", "POSTN", "TIMP1", "SERPINE1", "MMP9", "CTHRC1", "COL3A1",
    "PIIINP",
}

# Literature support flag — proteins repeatedly named in liver fibrosis CCI /
# scar-niche papers. Intentionally conservative: presence boosts score,
# absence does not penalize beyond zero.
LITERATURE_FIBROSIS = {
    "TREM2", "CD9", "GPNMB", "SPP1", "ACKR1", "PLVAP", "PDGFRB", "PDGFRA",
    "TNFRSF12A", "TNFSF12", "LOXL2", "TIMP1", "ACTA2", "COL1A1", "COL3A1",
    "CTHRC1", "POSTN", "CCR2", "CSF1R", "TGFB1",
}


# ---------------------------------------------------------------------------
# Sample-sheet parsing
# ---------------------------------------------------------------------------

# Per-sample folder name convention: GSM<digits>_<donor>_<fraction>
_SAMPLE_RE = re.compile(
    r"^(?P<gsm>GSM\d+)_(?P<donor>(?:healthy|cirrhotic|blood)\d+)_(?P<fraction>cd45[+\-][A-Za-z]?)$"
)


@dataclass(frozen=True)
class SampleInfo:
    gsm: str
    donor: str
    fraction: str
    disease: str            # "healthy" or "cirrhotic" (or "blood")
    tissue: str             # "liver" or "PBMC"
    path: Path


def parse_sample_dir(path: Path) -> Optional[SampleInfo]:
    """Parse a per-sample folder; return None if it does not match GSE136103."""
    m = _SAMPLE_RE.match(path.name)
    if m is None:
        return None
    donor = m.group("donor")
    if donor.startswith("healthy"):
        disease, tissue = "healthy", "liver"
    elif donor.startswith("cirrhotic"):
        disease, tissue = "cirrhotic", "liver"
    else:
        disease, tissue = "blood", "PBMC"
    return SampleInfo(
        gsm=m.group("gsm"),
        donor=donor,
        fraction=m.group("fraction"),
        disease=disease,
        tissue=tissue,
        path=path,
    )


def discover_samples(root: Path,
                     include_pbmc: bool = False) -> List[SampleInfo]:
    """Scan a directory of per-sample 10x folders and build a sample sheet."""
    out: List[SampleInfo] = []
    for sub in sorted(p for p in root.iterdir() if p.is_dir()):
        info = parse_sample_dir(sub)
        if info is None:
            continue
        if info.tissue == "PBMC" and not include_pbmc:
            continue
        out.append(info)
    return out


def build_sample_sheet(samples: Iterable[SampleInfo]) -> pd.DataFrame:
    rows = [{
        "sample_id": s.gsm + "_" + s.donor + "_" + s.fraction,
        "gsm": s.gsm,
        "donor": s.donor,
        "fraction": s.fraction,
        "disease": s.disease,
        "tissue": s.tissue,
        "path": str(s.path),
    } for s in samples]
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 10x loader
# ---------------------------------------------------------------------------

def _read_genes(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None,
                         names=["gene_id", "gene_symbol"], usecols=[0, 1])
    return df


def _read_barcodes(path: Path) -> List[str]:
    with gzip.open(path, "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def load_sample(info: SampleInfo) -> ad.AnnData:
    """Read one 10x mtx sample into AnnData (cells x genes)."""
    import scipy.io as sio
    mtx_path = info.path / "matrix.mtx.gz"
    with gzip.open(mtx_path, "rb") as fh:
        mat = sio.mmread(fh).tocsr()           # genes x cells
    barcodes = _read_barcodes(info.path / "barcodes.tsv.gz")
    genes = _read_genes(info.path / "features.tsv.gz")

    # Use gene symbols as var_names; resolve duplicates by ensembl-suffix
    var = genes.set_index("gene_symbol")
    if not var.index.is_unique:
        var.index = pd.Index(
            [f"{s}|{gid}" if c > 1 else s
             for s, c, gid in zip(var.index,
                                  var.groupby(level=0).cumcount() + 1,
                                  var["gene_id"])]
        )

    obs = pd.DataFrame({
        "barcode": barcodes,
        "sample_id": info.gsm + "_" + info.donor + "_" + info.fraction,
        "donor": info.donor,
        "fraction": info.fraction,
        "disease": info.disease,
        "tissue": info.tissue,
        "gsm": info.gsm,
    })
    obs.index = (obs["sample_id"] + "_" + obs["barcode"]).values

    X = mat.T.tocsr().astype(np.float32)        # cells x genes
    a = ad.AnnData(X=X, obs=obs, var=var)
    return a


def load_concat(samples: Iterable[SampleInfo],
                max_cells_per_sample: Optional[int] = None,
                random_state: int = 0) -> ad.AnnData:
    """Load several samples and concatenate. Optional per-sample subsampling
    for fast end-to-end demonstration runs on commodity hardware."""
    rng = np.random.default_rng(random_state)
    parts: List[ad.AnnData] = []
    for s in samples:
        a = load_sample(s)
        if max_cells_per_sample is not None and a.n_obs > max_cells_per_sample:
            idx = rng.choice(a.n_obs, size=max_cells_per_sample, replace=False)
            a = a[np.sort(idx)].copy()
        parts.append(a)
    if not parts:
        raise RuntimeError("No samples loaded.")
    merged = ad.concat(parts, join="outer", merge="first", index_unique=None)
    merged.var_names_make_unique()
    return merged


# ---------------------------------------------------------------------------
# Pseudobulk aggregation
# ---------------------------------------------------------------------------

def pseudobulk(adata: ad.AnnData,
               donor_key: str = "donor",
               group_key: str = "cell_type",
               layer: Optional[str] = "counts",
               min_cells: int = 10) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Sum counts per (donor x group). Returns (counts_df [genes x samples],
    sample_meta_df [samples x metadata])."""
    X = adata.layers[layer] if layer else adata.X
    if sparse.issparse(X):
        X = X.tocsr()
    obs = adata.obs[[donor_key, group_key, "disease"]].copy()
    obs["key"] = obs[donor_key].astype(str) + "__" + obs[group_key].astype(str)
    counts: Dict[str, np.ndarray] = {}
    meta_rows = []
    for key, idx in obs.groupby("key", observed=True).groups.items():
        rows = adata.obs_names.get_indexer(idx)
        if len(rows) < min_cells:
            continue
        sub = X[rows]
        col_sum = np.asarray(sub.sum(axis=0)).ravel() if sparse.issparse(sub) \
            else sub.sum(axis=0)
        counts[key] = col_sum.astype(np.int64)
        donor, group = key.split("__", 1)
        meta_rows.append({
            "pseudobulk_id": key,
            "donor": donor,
            "cell_type": group,
            "n_cells": int(len(rows)),
            "disease": adata.obs.loc[idx, "disease"].iloc[0],
        })
    cdf = pd.DataFrame(counts, index=adata.var_names).fillna(0).astype(np.int64)
    mdf = pd.DataFrame(meta_rows).set_index("pseudobulk_id").loc[cdf.columns]
    return cdf, mdf


# ---------------------------------------------------------------------------
# Scoring helpers
# ---------------------------------------------------------------------------

def minmax(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    lo, hi = np.nanmin(x), np.nanmax(x)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < 1e-12:
        return np.zeros_like(x)
    return (x - lo) / (hi - lo)


def safe_log2fc(a: np.ndarray, b: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    return np.log2((np.asarray(a) + eps) / (np.asarray(b) + eps))


def specificity_tau(expr_by_group: pd.DataFrame) -> pd.Series:
    """Tau specificity index across cell groups.
    expr_by_group: genes x groups mean expression (or detection rate)."""
    x = expr_by_group.values.astype(float)
    row_max = x.max(axis=1, keepdims=True)
    row_max[row_max == 0] = 1.0
    x_hat = x / row_max
    tau = (1.0 - x_hat).sum(axis=1) / max(x.shape[1] - 1, 1)
    return pd.Series(tau, index=expr_by_group.index, name="tau_specificity")
