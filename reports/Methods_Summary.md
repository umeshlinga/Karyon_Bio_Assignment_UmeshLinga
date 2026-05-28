# Methods Summary

This document distills the design and parameter choices made by each script. It is
deliberately compact; the screening-question answers in
`answers/Written_Screening_Questions.md` discuss the trade-offs in more depth.

## 1. Data and curation

- **Source**: GSE136103 RAW tar (~416 MB) from GEO. We extracted only liver
  samples (`*healthy*`, `*cirrhotic*`), skipping PBMC and mouse data, into one
  10x-style folder per sample.
- **Sample sheet**: parsed from the GSM-prefixed folder name
  (`GSM<digits>_<donor>_<fraction>`). Disease/tissue are inferred from the donor
  string. 20 sample-fractions across 10 donors are retained (5 healthy, 5
  cirrhotic).

## 2. QC (script 01)

- Per-sample 10x mtx loaded via `scipy.io.mmread`, gene symbols deduplicated by
  appending the Ensembl ID when collisions occur.
- QC metrics: `n_genes_by_counts`, `total_counts`, `pct_counts_mt` (MT-*),
  `pct_counts_ribo` (RPS*/RPL*), `pct_counts_hb` (HB[ABDEGMQZ]\d?).
- **Thresholds (liver-fibrosis-aware)**:
  - `min_genes >= 200`, `min_counts >= 500`
  - `pct_counts_mt <= 20` (relaxed; stressed HSCs and SAMs legitimately
    elevated)
  - `pct_counts_hb <= 50`, `pct_counts_ribo <= 60`
- **Doublets**: Scrublet via `sc.pp.scrublet` (batch_key=sample_id); flagged
  cells removed before downstream.
- Raw counts stored as `adata.layers["counts"]` before normalization.

## 3. Integration (script 02)

- `sc.pp.normalize_total(target_sum=1e4)` then `sc.pp.log1p`.
- HVG selection: `sc.pp.highly_variable_genes(n_top_genes=3000,
  batch_key="donor")` — donor-aware to avoid one donor dominating the gene set.
- Scale (`max_value=10`) + PCA (40 PCs, arpack).
- Batch correction: Harmony (`harmonypy.run_harmony` on `X_pca`, batch=donor).
  Falls back to BBKNN if Harmony missing; falls back to raw PCA otherwise.
- Neighbors (k=15), UMAP (min_dist=0.3), Leiden (resolution=0.8).
- **Disease is preserved**, not batch-corrected. Donor is the integration
  variable.

## 4. Annotation & compartment validation (script 03)

- Per-cell marker-panel scores via `sc.tl.score_genes` on each panel from
  `utils.LIVER_MARKERS` (HSC quiescent / activated, Kupffer, Monocyte, SAM,
  LSEC, vascular EC, Scar-EC, Cholangiocyte, Hepatocyte, T, B, Plasma, NK, pDC,
  Mast).
- Per-cluster mean scores → cluster labeled by argmax panel, marked "mixed" if
  top–second gap < 0.05, "Unknown" if top score < 0.05.
- Compartments collapsed into Mesenchyme / Macrophage / Endothelial /
  Cholangiocyte / Hepatocyte / T-NK / B-Plasma / DC / Mast.
- **Validation**: per-cell signature scores (SAM, Scar-EC, HSC activation, ECM,
  TGFβ) collapsed to donor means; cirrhotic vs healthy compared with a
  one-sided Mann-Whitney U.

## 5. Donor-aware DE + pathway + CCI (script 04)

- **Pseudobulk**: counts summed per `(donor × compartment)` (min 20 cells per
  group). For each of {Mesenchyme, Macrophage, Endothelial} we run PyDESeq2
  with design `~disease` and contrast `cirrhotic vs healthy`.
- Outputs: `results/de/pseudobulk_<comp>.csv` (full DE table),
  `de_top_<comp>.csv` (top up/down filtered at padj<0.05 & |log2FC|>1), volcano
  PNG.
- **Signature enrichment**: per compartment, for each fibrosis signature, we
  report n genes present, mean log2FC, fraction up-significant, top 5 genes by
  fold change.
- **CCI**: a curated panel of ligand-receptor pairs known to operate in the
  fibrotic niche (PDGF, TGFβ, TWEAK, Osteopontin, chemokines CXCL12/CCL2,
  NOTCH, ECM-integrin, CSF1) is scored by ligand and receptor compartment
  means in cirrhotic vs healthy; only pairs with detectable expression in both
  cirrhotic compartments are retained. Score = `ligand_mean_cirrh *
  receptor_mean_cirrh + max(ligand_bias, 0) + max(receptor_bias, 0)`. This is
  a *hypothesis generator*, not a permutation-based CCI test (see Q7 of the
  written answers for the trade-off discussion).

## 6. Biomarker prioritization (script 05)

- For each significant `(gene, compartment)` row (padj < 0.05 & |log2FC| > 0.5)
  we engineer:
  - `abs_log2FC`, `-log10(padj)`
  - `tau_specificity`: how specific the gene's mean expression is to one
    compartment vs all others
  - `detection_rate_cirrh`: fraction of cirrhotic cells in that compartment
    expressing the gene
  - boolean flags: `pathway_relevance` (in any fibrosis signature),
    `druggability` (curated DGIdb/Open Targets subset), `surface_or_secreted`,
    `literature_fibrosis`
  - `cci_evidence`: max LR interaction score where the gene appears as ligand
    or receptor
- **Rule-based score** (weights derived from the screening-answer rationale):
  `0.20·tau + 0.20·|log2FC| + 0.10·-log10(padj) + 0.10·detection +
  0.10·pathway + 0.15·druggability + 0.10·lit + 0.05·CCI`. Continuous features
  min-max normalized.
- **ML cross-check**: RandomForestRegressor on the same features with
  `literature_fibrosis` as a noisy positive label; predictions min-max
  normalized.
- **Final score = 0.7 · rule + 0.3 · ML**. Top-K kept per unique gene
  (preferring its highest-scoring compartment).
- **Translational category** assigned by a small decision rule based on
  `surface_or_secreted` × `druggability` × `pathway_relevance`.

## 7. Reproducibility & limitations

- Seeded: `random_state=0` for Scrublet, HVG, PCA, UMAP, Leiden, score_genes,
  and the Random Forest.
- All numbers in `reports/Executive_Summary.pdf` are read from
  `results/biomarker_ranked_table.csv`, `compartment_validation.csv`,
  `pathway_signature_scores.csv`, and `qc_per_sample.csv`. No hand-curated
  values.
- For the included run we used `--max_cells_per_sample 1500` (subsample for
  laptop-grade memory). Drop the flag for the full ~95 k-cell run.
- Curated lists (`DRUGGABLE_GENES`, `SURFACE_OR_SECRETED`,
  `LITERATURE_FIBROSIS`) are snapshots; production deployment should swap them
  for live DGIdb / Open Targets / Reactome API queries.
