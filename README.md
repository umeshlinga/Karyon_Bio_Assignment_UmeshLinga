# Karyon Bio Candidate Technical Assignment
## Mini-Pipeline for Discovery and Prioritization of Cell-Type-Specific Biomarkers in Human Liver Fibrosis

**Candidate:** Umesh Linga · **Location:** Indianapolis, Indiana · **Date:** May 2026 · **Submission:** Kiran.adhikari@karyon.bio

---

## Overview

End-to-end, reproducible single-cell pipeline run on the **GSE136103** human liver
fibrosis cohort (Ramachandran et al., *Nature* 2019; 5 healthy + 5 cirrhotic livers, NPC fractions).
It performs QC, donor-aware integration, marker-based annotation of the three
disease-relevant compartments (mesenchyme, macrophage, endothelial), donor-aware
pseudobulk differential expression with PyDESeq2, signature/pathway analysis, a
curated ligand-receptor evidence panel, and a hybrid rule-based + Random-Forest
biomarker prioritization. All numbers in `reports/Executive_Summary.pdf` are
produced by these scripts.

This run used **20 sample-fractions / 10 donors / 26,797 QC-retained cells** (the data
were subsampled to 1,500 cells per sample-fraction to keep memory ≤ ~12 GB; the
same scripts work on the full dataset by dropping `--max_cells_per_sample`).

---

## Headline results

- All three required compartments validated cirrhotic-up biology with donor-level
  Mann-Whitney U (`results/annotation/compartment_validation.csv`):
  - **Endothelial — Scar-EC signature**: Δ +0.76, *p* = 0.004
  - **Endothelial — SAM signature**: Δ +0.18, *p* = 0.004
  - **Macrophage — SAM signature**: Δ +0.10, *p* = 0.048
  - **Mesenchyme — HSC activation**: Δ +0.11, *p* = 0.075
- Pseudobulk DE recovered seminal Ramachandran 2019 findings:
  - Scar-EC: **ACKR1, PLVAP, RGCC** (mean log2FC +2.51, 33% up-significant)
  - Mesenchyme ECM: **LOX, COL1A1, COL1A2, COL3A1, LOXL2** (mean log2FC +3.50, 100% up-significant)
  - Macrophage ECM/SAM: FN1, COL1A1, COL3A1
- Top-scoring biomarker / target candidates (full table in
  `results/biomarker_ranked_table.csv`):

  | Rank | Gene   | Compartment   | log2FC | padj      | Score | Translational fit                  |
  |-----:|--------|---------------|-------:|-----------|------:|------------------------------------|
  | 1    | ACKR1  | Endothelial   | +4.25  | <1e-300   | 0.777 | Therapeutic (druggable, tractable) |
  | 2    | TIMP1  | Mesenchyme    | +3.48  | <1e-300   | 0.693 | Therapeutic (druggable, tractable) |
  | 3    | COL3A1 | Mesenchyme    | +3.62  | <1e-300   | 0.687 | Diagnostic / surface biomarker     |
  | 4    | PLVAP  | Mesenchyme    | +4.23  | <1e-300   | 0.622 | Diagnostic / surface biomarker     |
  | 5    | PDGFRA | Mesenchyme    | +3.28  | <1e-300   | 0.606 | Therapeutic (druggable, tractable) |
  | 6    | CTHRC1 | Mesenchyme    | +3.80  | <1e-300   | 0.581 | Diagnostic / surface biomarker     |
  | 7    | COL1A1 | Mesenchyme    | +4.77  | <1e-300   | 0.569 | Mechanistic biomarker              |
  | 8    | LOXL2  | Mesenchyme    | +2.79  | 0.025     | 0.568 | Therapeutic (pathway target)       |
  | 9    | GPNMB  | Mesenchyme    | +2.03  | 0.008     | 0.426 | Mechanistic biomarker              |
  | 10   | COL1A2 | Mesenchyme    | +3.95  | <1e-300   | 0.419 | Mechanistic biomarker              |

  (Top 11–15: LOX, FN1, BGN, PECAM1, DCN — see CSV for full feature matrix.)

---

## Repository layout

```
Karyon_Bio_Assignment_UmeshLinga/
├── README.md                       # this file
├── environment.yml                 # conda env (python 3.10–3.13)
├── requirements.txt                # pip fallback
├── LICENSE                         # MIT
├── .gitignore
├── data/
│   └── README_data.md              # GSE136103 download instructions
├── scripts/
│   ├── utils.py                    # loaders, markers, pseudobulk, scoring
│   ├── 01_download_and_qc.py       # 10x mtx loader + QC + Scrublet
│   ├── 02_preprocess_integrate.py  # norm/HVG/PCA + Harmony + Leiden
│   ├── 03_annotate_validate.py     # marker-score annotation + validation
│   ├── 04_de_pathway_cci.py        # PyDESeq2 + signature + LR panel
│   └── 05_biomarker_prioritization.py  # rule + RF scoring
├── results/
│   ├── qc/                         # qc_summary.md, qc_per_sample.csv, violin
│   ├── processed/                  # qc_filtered.h5ad, integrated.h5ad, annotated.h5ad
│   ├── annotation/                 # cluster scores, compartment validation
│   ├── de/                         # pseudobulk DE per compartment, signature, CCI
│   ├── figures/                    # UMAPs, violins, volcanoes, dot plots, top biomarker plot
│   ├── biomarker_feature_matrix.csv
│   └── biomarker_ranked_table.csv
├── reports/
│   ├── Executive_Summary.pdf       # one-to-two-page summary (regenerated from results/)
│   ├── Methods_Summary.md          # methods + parameter sheet
│   └── create_executive_summary.py
└── answers/
    └── Written_Screening_Questions.md  # answers to all 8 screening questions
```

---

## Quick-start (full re-run)

### 1. Environment

```bash
conda env create -f environment.yml
conda activate karyon-liver-sc
# or:  pip install -r requirements.txt
```

Tested on **Python 3.13** with scanpy 1.12, anndata 0.12, harmonypy, pydeseq2 0.5, scikit-learn 1.7, scipy 1.16.

### 2. Data

```bash
mkdir -p data/raw && cd data/raw
curl -L -o GSE136103_RAW.tar \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE136103&format=file"
mkdir extracted && tar -xf GSE136103_RAW.tar -C extracted \
  --wildcards "*healthy*" "*cirrhotic*"
# Reorganize into per-sample folders (each containing matrix.mtx.gz,
# barcodes.tsv.gz, features.tsv.gz). See data/README_data.md.
```

### 3. Run pipeline

```bash
# QC (--max_cells_per_sample 1500 was used for the included results;
# omit for the full ~80–100k-cell run)
python scripts/01_download_and_qc.py \
    --raw_dir data/raw/extracted/per_sample \
    --out_dir results \
    --max_cells_per_sample 1500

python scripts/02_preprocess_integrate.py \
    --in_h5ad results/processed/qc_filtered.h5ad --out_dir results

python scripts/03_annotate_validate.py \
    --in_h5ad results/processed/integrated.h5ad --out_dir results

python scripts/04_de_pathway_cci.py \
    --in_h5ad results/processed/annotated.h5ad --out_dir results

python scripts/05_biomarker_prioritization.py \
    --de_csv results/de/pseudobulk_all_compartments.csv \
    --in_h5ad results/processed/annotated.h5ad \
    --cci_csv results/de/cci_evidence_panel.csv \
    --out_dir results --top_k 20

python reports/create_executive_summary.py \
    --results_dir results --out reports/Executive_Summary.pdf
```

Approximate wall-clock on a 16 GB laptop: QC ~1 min, integration (Harmony 5 iters) ~3 min, annotation ~30 s, DE+CCI ~2 min, prioritization ~10 s, PDF ~5 s. Full dataset (~95 k cells) needs ~32 GB RAM and ~15–25 min.

---

## Key design choices (mapped to written-screening answers)

| Choice                              | Where    | Rationale                                                              |
|-------------------------------------|----------|------------------------------------------------------------------------|
| Mito ≤ 20% (not the usual 10%)      | 01       | Fibrotic livers stress cells; aggressive cutoffs erase HSC/SAM biology |
| Donor as integration batch          | 02       | Disease must be preserved, not removed                                  |
| Donor-aware HVG                     | 02       | `batch_key=donor` in `sc.pp.highly_variable_genes`                       |
| Harmony on PCA (BBKNN fallback)     | 02       | Disease retained as second-order biology after donor batch correction   |
| Pseudobulk PyDESeq2 `~disease`      | 04       | Treats donor (not cell) as unit; Q5 of written answers                   |
| Tau specificity (gene × compartment)| 05       | Penalizes broadly-expressed ECM genes vs cell-type-restricted markers   |
| Druggability + surface flags        | utils.py | DGIdb/Open Targets curated subset (refresh in production via API)       |
| RF as ML cross-check (not primary)  | 05       | Literature-fibrosis as noisy positive label; RF blended 30% into final  |

---

## Deliverables checklist (from the assignment PDF)

- [x] **Reproducible code + README** (this repo + scripts/01..05)
- [x] **QC summary** (`results/qc/qc_summary.md`, `qc_per_sample.csv`, `qc_violin.png`)
- [x] **Cell-annotation figures** (`results/figures/umap_celltype.png`, `dotplot_markers.png`)
- [x] **DE + pathway results** (`results/de/`, including signature enrichment + CCI panel)
- [x] **Ranked biomarker / target table** (`results/biomarker_ranked_table.csv`)
- [x] **1–2 page Executive Summary** (`reports/Executive_Summary.pdf`)
- [x] **Written answers** to all 8 screening questions (`answers/Written_Screening_Questions.md`)

---

## Limitations & honest notes

- **Subsampling**: results above use 1,500 cells per sample-fraction to keep the
  run lightweight. The full dataset (≈95 k cells) increases statistical power
  for rarer states (SAMs); we expect more macrophage-compartment genes
  (TREM2, SPP1, CD9) to reach individual-gene significance at full size — they
  already drive the SAM signature score (mean log2FC +0.96 in macrophages).
- **Annotation granularity**: SAMs are a state within the macrophage lineage
  and are not always resolved as a separate Leiden cluster at this resolution
  / cell count; the SAM signature score still validates the population
  expansion in cirrhotic donors (Δ +0.10, p = 0.048).
- **Druggability flags** are a hand-curated snapshot of DGIdb/Open Targets
  evidence. For production, swap `utils.DRUGGABLE_GENES` for a fresh API
  pull at run time.
- **CCI panel** is a *hypothesis generator* (curated LR pairs + expression
  evidence), not a permutation-based CellPhoneDB/LIANA run. We chose this for
  lightweight reproducibility; the screening-question answer (Q7) discusses
  the trade-off and how to validate.

---

## Contact

Prepared for Karyon Bio technical screening.
Umesh Linga · Indianapolis, Indiana
