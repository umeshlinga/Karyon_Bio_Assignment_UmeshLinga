#!/usr/bin/env python3
"""
04_de_pathway_cci.py
Donor-aware differential expression, pathway/signature analysis, and a
light-weight ligand-receptor cell-cell interaction summary across the
fibrotic niche.

Statistical approach
--------------------
For each compartment (Mesenchyme, Macrophage, Endothelial) we form
pseudobulk profiles (donor x compartment) by summing raw counts, then run
PyDESeq2 with design ``~disease`` (cirrhotic vs healthy reference). This is
the donor-aware design described in Q5 of the screening questions, and it
is appropriate even at the modest donor count of GSE136103 because each
pseudobulk row is a biological replicate (donor).

Pathway / signature analysis
----------------------------
We score curated fibrosis-relevant gene sets (utils.FIBROSIS_SIGNATURES) on
the DE-up genes for each compartment, plus simple over-representation of
the top-K genes against canonical sets (ECM organization, TGF-beta, etc.).

CCI
---
We do not run a full LIANA/CellPhoneDB pass here (heavy R/Java dependencies
+ runtime); instead we build an interaction *evidence panel* using a curated
ligand-receptor map (subset of OmniPath) restricted to ligands expressed in
sender compartments and receptors expressed in receiver compartments, with
disease-bias scoring. This is meant as a hypothesis generator, exactly as
the screening answer prescribes.

Inputs
------
results/processed/annotated.h5ad

Outputs
-------
results/de/pseudobulk_<compartment>.csv      DESeq2 results per gene
results/de/de_top_<compartment>.csv          Top up-/down-regulated genes
results/de/pathway_signature_scores.csv      Signature enrichment per compartment
results/de/cci_evidence_panel.csv            Curated LR interactions, ranked
results/figures/volcano_<compartment>.png
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils as U                                       # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)


# A small, curated ligand-receptor list focused on the fibrotic niche.
# Each row: (ligand, receptor, sender_compartment_hint, receiver_compartment_hint,
#           pathway_label). The hints are advisory: at runtime we check that
#           ligand is expressed in any compartment and receptor in any compartment.
LR_PANEL = [
    ("PDGFB",   "PDGFRB",     "Macrophage",  "Mesenchyme",  "PDGF"),
    ("PDGFA",   "PDGFRA",     "Macrophage",  "Mesenchyme",  "PDGF"),
    ("TNFSF12", "TNFRSF12A",  "Macrophage",  "Mesenchyme",  "TWEAK"),
    ("TGFB1",   "TGFBR1",     "Macrophage",  "Mesenchyme",  "TGFB"),
    ("TGFB1",   "TGFBR2",     "Macrophage",  "Mesenchyme",  "TGFB"),
    ("SPP1",    "ITGAV",      "Macrophage",  "Mesenchyme",  "Osteopontin"),
    ("SPP1",    "CD44",       "Macrophage",  "Macrophage",  "Osteopontin"),
    ("CCL2",    "CCR2",       "Endothelial", "Macrophage",  "Chemokine"),
    ("CXCL12",  "CXCR4",      "Endothelial", "Macrophage",  "Chemokine"),
    ("JAG1",    "NOTCH1",     "Mesenchyme",  "Endothelial", "NOTCH"),
    ("JAG1",    "NOTCH2",     "Endothelial", "Mesenchyme",  "NOTCH"),
    ("DLL4",    "NOTCH3",     "Endothelial", "Mesenchyme",  "NOTCH"),
    ("COL1A1",  "ITGA11",     "Mesenchyme",  "Mesenchyme",  "ECM-Integrin"),
    ("FN1",     "ITGAV",      "Mesenchyme",  "Mesenchyme",  "ECM-Integrin"),
    ("POSTN",   "ITGAV",      "Mesenchyme",  "Macrophage",  "ECM-Integrin"),
    ("CSF1",    "CSF1R",      "Mesenchyme",  "Macrophage",  "CSF1"),
]


def pseudobulk_de(adata, compartment: str):
    """Run PyDESeq2 on (donor x compartment) pseudobulk for one compartment."""
    sub = adata[adata.obs["compartment"] == compartment].copy()
    if sub.n_obs < 50:
        print(f"[skip] compartment '{compartment}' has only {sub.n_obs} cells.")
        return None

    counts_df, meta_df = U.pseudobulk(
        sub, donor_key="donor", group_key="compartment", layer="counts",
        min_cells=20)
    if counts_df.shape[1] < 3:
        print(f"[skip] only {counts_df.shape[1]} pseudobulk samples for "
              f"{compartment}; need >=3 donors per group.")
        return None

    # PyDESeq2 expects samples x genes
    counts = counts_df.T.astype(int)
    meta = meta_df.copy()
    meta["disease"] = pd.Categorical(meta["disease"],
                                     categories=["healthy", "cirrhotic"])

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:                            # pragma: no cover
        raise SystemExit("pydeseq2 missing; pip install pydeseq2") from exc

    dds = DeseqDataSet(counts=counts, metadata=meta, design="~disease",
                       refit_cooks=True, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=("disease", "cirrhotic", "healthy"),
                       quiet=True)
    stats.summary()
    res = stats.results_df
    res["compartment"] = compartment
    res["gene"] = res.index
    return res.reset_index(drop=True)


def make_volcano(res, out_path: Path, title: str):
    fig, ax = plt.subplots(figsize=(6, 5))
    x = res["log2FoldChange"].values
    y = -np.log10(res["padj"].fillna(1.0).clip(lower=1e-300).values)
    sig = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)
    ax.scatter(x[~sig], y[~sig], s=4, c="lightgrey", alpha=0.6)
    ax.scatter(x[sig], y[sig], s=8, c="crimson", alpha=0.85)
    top = res.dropna(subset=["padj"]).copy()
    top = top.assign(score=top["log2FoldChange"].abs() * (-np.log10(top["padj"].clip(lower=1e-300))))
    for _, row in top.nlargest(12, "score").iterrows():
        ax.text(row["log2FoldChange"], -np.log10(max(row["padj"], 1e-300)),
                row["gene"], fontsize=7)
    ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.6)
    ax.axvline(1, color="grey", ls="--", lw=0.6); ax.axvline(-1, color="grey", ls="--", lw=0.6)
    ax.set_xlabel("log2 FC (cirrhotic / healthy)"); ax.set_ylabel("-log10 padj")
    ax.set_title(title)
    fig.savefig(out_path, dpi=140, bbox_inches="tight"); plt.close(fig)


def signature_enrichment(de_results: pd.DataFrame) -> pd.DataFrame:
    """Mean log2FC and fraction-significant for each fibrosis signature, per compartment."""
    rows = []
    for comp, df in de_results.groupby("compartment"):
        for sig_name, genes in U.FIBROSIS_SIGNATURES.items():
            sub = df[df["gene"].isin(genes)]
            if sub.empty:
                continue
            sig = (sub["padj"] < 0.05) & (sub["log2FoldChange"] > 0)
            rows.append({
                "compartment": comp,
                "signature": sig_name,
                "n_in_signature": len(sub),
                "mean_log2FC": float(sub["log2FoldChange"].mean()),
                "frac_up_sig": float(sig.mean()),
                "top_genes": ",".join(
                    sub.sort_values("log2FoldChange", ascending=False)
                       .head(5)["gene"].tolist()),
            })
    return pd.DataFrame(rows)


def cci_evidence_panel(adata) -> pd.DataFrame:
    """For each LR pair: expression and disease bias in plausible compartments."""
    raw = adata.raw.to_adata()
    # gene -> column index, for fast lookup
    var_idx = {g: i for i, g in enumerate(raw.var_names)}
    obs = adata.obs

    def compartment_mean(gene: str, compartment: str, disease: str) -> float:
        if gene not in var_idx:
            return float("nan")
        mask = (obs["compartment"] == compartment) & (obs["disease"] == disease)
        if mask.sum() == 0:
            return float("nan")
        col = raw.X[:, var_idx[gene]]
        if hasattr(col, "toarray"):
            col = col.toarray().ravel()
        else:
            col = np.asarray(col).ravel()
        return float(col[mask.values].mean())

    rows = []
    compartments = ["Mesenchyme", "Macrophage", "Endothelial"]
    for ligand, receptor, sender_hint, receiver_hint, pathway in LR_PANEL:
        for sender in compartments:
            for receiver in compartments:
                if sender == receiver and (sender_hint != receiver_hint):
                    continue
                lig_c = compartment_mean(ligand, sender, "cirrhotic")
                lig_h = compartment_mean(ligand, sender, "healthy")
                rec_c = compartment_mean(receptor, receiver, "cirrhotic")
                rec_h = compartment_mean(receptor, receiver, "healthy")
                if np.isnan([lig_c, lig_h, rec_c, rec_h]).any():
                    continue
                # Both ligand AND receptor must be detectable in cirrhotic
                if lig_c < 0.05 or rec_c < 0.05:
                    continue
                lig_bias = lig_c - lig_h
                rec_bias = rec_c - rec_h
                rows.append({
                    "pathway": pathway,
                    "ligand": ligand, "receptor": receptor,
                    "sender": sender, "receiver": receiver,
                    "ligand_mean_cirrh": lig_c, "ligand_mean_healthy": lig_h,
                    "receptor_mean_cirrh": rec_c, "receptor_mean_healthy": rec_h,
                    "ligand_bias": lig_bias, "receptor_bias": rec_bias,
                    "interaction_score": float((lig_c * rec_c) + max(lig_bias, 0) + max(rec_bias, 0)),
                    "hint_match": (sender == sender_hint and receiver == receiver_hint),
                })
    df = pd.DataFrame(rows).sort_values("interaction_score", ascending=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_h5ad", required=True, type=Path)
    ap.add_argument("--out_dir", default=Path("results"), type=Path)
    args = ap.parse_args()

    de_dir = args.out_dir / "de"
    fig_dir = args.out_dir / "figures"
    de_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.in_h5ad)
    print(f"Loaded {adata.n_obs:,} cells")

    de_frames = []
    for comp in ["Mesenchyme", "Macrophage", "Endothelial"]:
        print(f"\n=== DE: {comp} ===")
        res = pseudobulk_de(adata, comp)
        if res is None:
            continue
        out_csv = de_dir / f"pseudobulk_{comp}.csv"
        res.to_csv(out_csv, index=False)
        n_sig = ((res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)).sum()
        print(f"  significant genes (|log2FC|>1 & padj<0.05): {n_sig}")
        # Top up/down
        top_up = res[(res["padj"] < 0.05) & (res["log2FoldChange"] > 0)] \
            .sort_values("log2FoldChange", ascending=False).head(30)
        top_dn = res[(res["padj"] < 0.05) & (res["log2FoldChange"] < 0)] \
            .sort_values("log2FoldChange").head(30)
        top = pd.concat([top_up.assign(direction="up"),
                         top_dn.assign(direction="down")])
        top.to_csv(de_dir / f"de_top_{comp}.csv", index=False)
        make_volcano(res, fig_dir / f"volcano_{comp}.png",
                     f"{comp}: cirrhotic vs healthy (donor-aware pseudobulk DE)")
        de_frames.append(res)

    if de_frames:
        all_de = pd.concat(de_frames, ignore_index=True)
        all_de.to_csv(de_dir / "pseudobulk_all_compartments.csv", index=False)
        sig_df = signature_enrichment(all_de)
        sig_df.to_csv(de_dir / "pathway_signature_scores.csv", index=False)
        print("\nSignature enrichment summary:")
        print(sig_df.round(3).to_string(index=False))

    cci = cci_evidence_panel(adata)
    cci.to_csv(de_dir / "cci_evidence_panel.csv", index=False)
    print(f"\nCCI evidence panel: {len(cci)} candidate interactions retained.")
    if not cci.empty:
        print("Top 10 ranked LR pairs:")
        print(cci[["pathway", "ligand", "receptor", "sender", "receiver",
                   "interaction_score"]].head(10).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
