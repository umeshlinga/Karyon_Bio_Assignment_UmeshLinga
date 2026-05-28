# Data Acquisition and Preparation for GSE136103

## Primary dataset

**GSE136103** — Resolving the fibrotic niche of human liver cirrhosis using
single-cell transcriptomics (Ramachandran et al., *Nature* 2019).

- Page:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136103
- File:  `GSE136103_RAW.tar` (~416 MB, 78 files = 26 samples × 3 files each)

## Step 1 — download

```bash
mkdir -p data/raw && cd data/raw
curl -L -o GSE136103_RAW.tar \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE136103&format=file"
```

Expected size: ~416 MB. Verifies in ~30–60 seconds on a 10 MB/s link.

## Step 2 — extract only the liver samples

The tarball contains liver samples (`*healthy*`, `*cirrhotic*`), 4 PBMC samples
(`*blood*`), and one mouse sample (`*mouse*`). For this assignment we use only
the human liver samples:

```bash
mkdir extracted
tar -xf GSE136103_RAW.tar -C extracted --wildcards "*healthy*" "*cirrhotic*"
```

This produces 63 files in a flat directory, named like:

```
GSM4041150_healthy1_cd45+_barcodes.tsv.gz
GSM4041150_healthy1_cd45+_genes.tsv.gz
GSM4041150_healthy1_cd45+_matrix.mtx.gz
GSM4041151_healthy1_cd45-A_barcodes.tsv.gz
...
```

## Step 3 — reorganize into one folder per sample

`scripts/01_download_and_qc.py` expects per-sample subdirectories, each
containing `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz` (rename
`genes.tsv.gz` -> `features.tsv.gz`):

```bash
cd extracted
mkdir per_sample
for f in *_matrix.mtx.gz; do
    sample="${f%_matrix.mtx.gz}"
    mkdir -p "per_sample/$sample"
    mv "${sample}_matrix.mtx.gz"   "per_sample/$sample/matrix.mtx.gz"
    mv "${sample}_barcodes.tsv.gz" "per_sample/$sample/barcodes.tsv.gz"
    mv "${sample}_genes.tsv.gz"    "per_sample/$sample/features.tsv.gz"
done
```

After this there will be 21 sample directories. `discover_samples` in
`scripts/utils.py` parses the directory name to recover donor, fraction,
disease, and tissue.

## Cohort overview (after extraction)

| Donor       | Fractions                      |
|-------------|--------------------------------|
| healthy1    | cd45+, cd45-A, cd45-B          |
| healthy2    | cd45+, cd45-                   |
| healthy3    | cd45+, cd45-A, cd45-B          |
| healthy4    | cd45+, cd45-                   |
| healthy5    | cd45+                          |
| cirrhotic1  | cd45+, cd45-A, cd45-B          |
| cirrhotic2  | cd45+, cd45-                   |
| cirrhotic3  | cd45+, cd45-A, cd45-B          |
| cirrhotic4  | cd45+, cd45-                   |
| cirrhotic5  | cd45+                          |

That matches the paper's "5 healthy + 5 cirrhotic liver donors" cohort, sorted
into CD45+ (immune) and CD45- (non-immune NPC) fractions.

## Step 4 — run the QC script

```bash
python scripts/01_download_and_qc.py \
    --raw_dir data/raw/extracted/per_sample \
    --out_dir results \
    --max_cells_per_sample 1500   # optional; remove for full dataset
```

## Notes

- The CD45+/CD45- split is the technical batch in this study, captured by the
  `fraction` field in the sample sheet. Donor is used as the integration batch
  (script 02).
- Large processed objects (`*.h5ad`) are git-ignored; regenerate via the
  scripts.
- Optional validation datasets (GSE244832, GSE207310, SCP2154) are not
  bundled; the pipeline is modular and can ingest them via the same loader.
