# Written Screening Questions - Responses
## Karyon Bio Candidate Technical Assignment

**Candidate:** Umesh Linga  
**Location:** Indianapolis, Indiana  
**Focus:** Human liver fibrosis scRNA-seq/snRNA-seq biomarker discovery

All answers are based on hands-on experience with liver and fibrosis single-cell datasets, best practices in the field (scanpy/Seurat ecosystems, pseudobulk methods, multi-sample integration), and specific knowledge of liver fibrosis biology and the GSE136103 dataset.

---

## 01 | Dataset curation and fibrosis-stage harmonization

**Approach:**

1. **Metadata Collection & Standardization**
   - Pull all available clinical/histological metadata from GEO, SRA, supplementary tables, and original publications.
   - Create a unified `samplesheet` with controlled vocabulary:
     - `fibrosis_stage`: Map to common scale — "healthy/F0", "F1", "F2", "F3", "F4/cirrhosis", or binary "advanced_fibrosis" (F3-F4) vs others. Use METAVIR where available; infer from "cirrhosis" labels or Ishak scores.
     - `etiology`: NASH/MASH, viral (HBV/HCV), alcohol, autoimmune, cryptogenic, etc. (critical confounder).
     - `donor_id`, `age`, `sex`, `BMI`, `comorbidities`, `treatment`.
   - Flag samples with incomplete metadata. Use multiple imputation or sensitivity analysis for missing fibrosis scores.

2. **Harmonization Strategy**
   - Use a **tiered labeling**:
     - Tier 1 (high confidence): Direct histological METAVIR or equivalent.
     - Tier 2: Clinical diagnosis of cirrhosis + imaging/labs.
     - Tier 3: Study-defined "fibrotic" vs "healthy".
   - Cross-validate labels by correlating with molecular readouts (e.g., bulk COL1A1 expression or known HSC activation signatures) where possible.
   - For NASH/MASH datasets, align to fibrosis stage rather than steatosis/inflammation alone, as the assignment focuses on fibrosis.

3. **Validation Steps**
   - **Internal consistency**: Check that fibrosis labels correlate with expected biology (e.g., higher ECM gene expression, expansion of myofibroblasts/SAMs in higher stages).
   - **Batch/donor effects**: Visualize metadata vs technical variables (sequencing depth, %mito, library prep).
   - **Inter-dataset alignment**: Use `harmony` or `scVI` with fibrosis_stage as a covariate of interest (not batch). Test preservation of disease signal post-integration (e.g., differential abundance of SAMs or HSC activation score still significant).
   - **Expert review**: If possible, consult hepatopathologist for ambiguous cases.
   - **Sensitivity**: Run key analyses with/without borderline samples or using continuous fibrosis scores if available (e.g., from histology AI).

4. **Practical Tools**
   - Python: `pandas` for curation + `anndata` obs metadata.
   - Automated mapping scripts with dictionaries for common label variants.
   - Version-controlled `metadata_curated.tsv` + provenance (original source + mapping rules).

This ensures downstream biomarker discovery is not confounded by inconsistent staging and maximizes power for F2+ or advanced fibrosis signals.

---

## 02 | QC and preprocessing for liver scRNA-seq/snRNA-seq

**Recommended QC Pipeline (Liver Fibrosis Specific):**

**Filtering thresholds (tuned, not default):**
- **Cells**: 
  - Min genes: 200–500 (lower for some stressed populations).
  - Max genes: 4,000–6,000 (liver cells, especially hepatocytes or large cells, can be high; avoid cutting activated HSCs).
  - Min counts: 500–1,000.
- **Mitochondrial %**: **Relaxed to 15–25%**. In fibrotic liver, stressed, apoptotic, or metabolically altered cells (including activated HSCs and some macrophages) legitimately have higher mito content. Aggressive filtering (e.g., <10%) removes disease-relevant biology. Use adaptive thresholds or model mito as a covariate.
- **Ribosomal % / Hemoglobin**: Filter high ribo if doublet-like; hemoglobin for blood contamination (important in liver dissociations).
- **Doublets**: 
  - `Scrublet` (preferred for Python) or `DoubletFinder` (R). 
  - Run per sample or per batch. Use expected doublet rate ~5–10% for 10x. 
  - **Important for liver**: Hepatocytes and some mesenchymal cells are large/polyploid → can appear doublet-like. Manually inspect suspicious clusters (high marker combination from two lineages) and retain if biologically plausible (e.g., transitional states).
- **Ambient RNA / Soup**: Use `SoupX` (R) or `scAR` / `DecontX`. Critical in fibrotic samples with high debris/ECM or necrotic tissue. Especially important for accurate macrophage and endothelial signals.
- **Genes**: Filter genes expressed in <3–10 cells. Remove obvious contaminants (e.g., high mitochondrial genes if modeling separately).

**Batch & Technical Correction:**
- Identify major technical variables (donor, digestion protocol, 10x chemistry/version, sorting (CD45+/-), sequencing depth).
- **Do not over-correct disease signal.** Use:
  - `scVI` (with `batch_key` and optionally `categorical_covariate_keys` including disease or fibrosis_stage).
  - `Harmony` (with `theta` tuned; validate biology preserved).
  - `scanpy.pp.combat` or `mnnpy` as lighter options.
- **Validation of correction**: 
  - Check that known disease markers (COL1A1 in HSCs, TREM2 in SAMs) remain DE.
  - Proportions of key populations (SAMs, activated mesenchyme) still differ significantly by disease status.
  - Use kBET, LISI, or silhouette scores for batch mixing vs biological separation.

**snRNA-seq Specific:**
- Nuclear vs cytoplasmic differences (e.g., lower detection of some membrane transcripts). Use sn-specific references or account in annotation.
- Often cleaner (less ambient), but may miss some cytoplasmic signals.

**Overall Philosophy:** Prioritize **retaining biological signal** from stressed/diseased cells over "clean" data. Iterate QC with biological validation (e.g., marker gene inspection, known population recovery from the Ramachandran paper).

---

## 03 | Dataset integration without erasing fibrosis biology

**Strategy to Preserve Disease Signal:**

1. **Avoid treating disease as a batch to remove.**
   - Never put `fibrosis_stage` or `disease` purely into the batch correction variable if the goal is to study it.

2. **Recommended Integration Approaches:**
   - **scVI / scArches (best for this use case)**: 
     - Train with `batch_key = "donor"` or technical batch.
     - Include `categorical_covariate_keys = ["disease_status"]` or use conditional VAEs.
     - Or train separate models per disease group and align latent spaces carefully.
   - **Harmony**: Run with technical batches; inspect PCs — ensure top PCs capture disease biology (e.g., HSC activation, SAM signature) before/after.
   - **Seurat RPCA or CCA anchors**: With `k.anchor` tuned; split by disease first if strong imbalance.
   - **Reference-based**: Use a high-quality healthy + cirrhotic reference (e.g., from this GSE or larger atlas) and project query datasets.

3. **Validation that Biology is Preserved:**
   - **Differential abundance / composition**: Test if proportion of SAMs, activated HSCs, or scar-ECs remains significantly different post-integration (using `scCODA`, `milo`, or simple Fisher/permutation tests per donor).
   - **Marker preservation**: Confirm top DE genes from paper (TREM2, CD9 in macrophages; ACKR1, PLVAP in EC; collagen genes in mesenchyme) are still detectable and cell-type specific.
   - **Embedding inspection**: UMAP colored by disease should still show separation or gradients in relevant compartments even after batch correction.
   - **Quantitative metrics**: Compare ARI/NMI for batch vs disease labels pre/post correction. Disease clustering should not collapse.
   - **Downstream robustness**: Key biomarker rankings or pathway enrichments should be stable with/without aggressive correction.

4. **Additional Best Practices:**
   - Integrate **within disease groups** first (healthy livers together, cirrhotic together), then align the two super-clusters if needed.
   - Use **donor as primary batch** (strongest source of variation usually).
   - For very heterogeneous datasets, consider **topic modeling** (LDA or scVI with topics) or **graph-based** methods that are more robust.
   - Always perform **sensitivity analysis**: Report results with minimal vs moderate correction.

This approach ensures we discover true fibrosis-associated states rather than technical artifacts or erased biology.

---

## 04 | Cell-type annotation and validation in fibrotic liver

**Scenario Analysis (the described cluster):**

The cluster expresses classic **activated stromal / myofibroblast markers**: `COL1A1`, `COL3A1`, `ACTA2` (αSMA), `TAGLN`, `PDGFRB`, `LUM`, `DCN`. Strong enrichment in F3/F4 samples → highly suggestive of **activated hepatic stellate cells (aHSCs) or myofibroblasts** contributing to scar formation.

**Validation Workflow:**

1. **Marker Gene Deep Dive**
   - Quiescent HSC markers (should be low or absent): `LRAT`, `RBP1`, `PPARG`, `CYGB`, `REEP6`, `RELN`, `GATA4/6`, `LHX2`.
   - Activated / myofibroblast: High `ACTA2`, `COL1A1/2`, `COL3A1`, `LOX`, `LOXL2`, `TIMP1`, `POSTN`, `SPARC`, `CTHRC1`, `PDGFRB`, `ITGA11`, etc.
   - Distinguish HSC vs Portal Fibroblasts (PFs):
     - HSC-enriched: Often perisinusoidal, express certain genes (e.g., `GPC3`, `DBH` subpopulations in some studies).
     - PFs: More portal tract associated; may express `ELN` (elastin), `WT1` (in some), different ECM profile.
     - In advanced fibrosis, lineages converge; subclustering + trajectory analysis helpful.

2. **Subclustering & Trajectory**
   - Re-cluster the stromal compartment at higher resolution.
   - Pseudotime / RNA velocity (scVelo) or trajectory inference (Monocle, PAGA, slingshot) to see if cells progress from quiescent-like → activated → myofibroblast states.
   - Check for multiple stromal states (e.g., central vs portal vein associated, or different activation trajectories).

3. **Cross-Dataset & Reference Validation**
   - Map to public liver atlases (Human Liver Cell Atlas, GSE136103 original annotations, or larger integrated objects).
   - Use `CellTypist`, `scArches`, or `Azimuth` with liver reference.
   - Compare marker overlap with published aHSC signatures from Ramachandran, Payen, etc.

4. **Functional & Spatial Validation**
   - **Pathway enrichment**: ECM organization, collagen biosynthesis, TGF-β signaling, focal adhesion, wound healing — expected for aHSCs/myofibroblasts.
   - **Spatial transcriptomics** (if available or orthogonal): Check localization to fibrotic septa/scars vs parenchyma. aHSCs and myofibroblasts are scar-associated.
   - **Proportion test**: Significantly higher abundance or activation score in F3/F4 vs F0/F1/healthy (donor-aware stats).
   - **Literature/orthogonal**: Correlate with bulk fibrosis gene signatures or IHC/IF for ACTA2+ / PDGFRB+ cells in scars.

5. **Conclusion for this cluster**
   - Most likely represents **activated hepatic stellate cells / myofibroblast-like cells** (the main ECM-producing population in liver fibrosis).
   - Could be mixed with portal fibroblast-derived myofibroblasts in advanced stages.
   - Label as **"Activated HSC / Myofibroblast"** or **"Scar-associated Mesenchyme (SAMes)"** with subclusters if resolved.
   - Report confidence and any ambiguity transparently.

This rigorous validation prevents misannotation that could derail biomarker discovery (e.g., attributing HSC markers to wrong lineage).

---

## 05 | Donor-aware differential expression and biomarker discovery

**Why simple cell-level DE is dangerous:**

- **Pseudoreplication**: Cells from the same donor are highly correlated (shared genetics, microenvironment, technical factors). Treating 10,000 cells from one donor as 10,000 independent samples massively inflates statistical power and produces many false positives.
- **Donor imbalance**: One or two outlier donors can drive apparent "significant" genes that are not generalizable.
- **Especially problematic** for macrophages and endothelial cells in fibrosis, where donor-to-donor heterogeneity in etiology, inflammation, and fibrosis stage is high.
- Leads to poor reproducibility and biomarkers that fail validation.

**Recommended Statistical Strategy (Donor-aware):**

1. **Primary Recommendation: Pseudobulk Aggregation**
   - Aggregate counts to **donor × cell_type** level (sum or median).
   - Use robust methods: `DESeq2`, `edgeR` (with TMM), or `limma-voom`.
   - Model: `~ fibrosis_stage + etiology + sex + age + donor` or use `duplicateCorrelation` / random effects for donor.
   - In Python: `pyDESeq2` or `decoupler` + statsmodels; or call R via `rpy2`.
   - Filter low-count genes/features properly.

2. **Advanced / Mixed-Effect Models**
   - `MAST` (single-cell aware hurdle model) with donor as random effect.
   - `muscat` or `dream` (variancePartition) frameworks in R/Bioconductor for multi-sample multi-celltype DE.
   - GLMMs via `glmmTMB` or Python equivalents.

3. **Best Practices for Biomarkers**
   - Require **consistency across donors**: e.g., gene significant in ≥70% of donors or use meta-analysis (combine p-values with Fisher's or Stouffer's method).
   - Report **donor-level effect sizes** and heterogeneity (I² statistic).
   - Use **cell-type-specific pseudobulk** for macrophages and endothelial cells separately.
   - For trajectory or state discovery within compartment: Still aggregate or use donor-stratified clustering.
   - Multiple testing: Stricter FDR or focus on effect size + specificity rather than p-value alone.

4. **Implementation Tip**
   - In the pipeline (`04_de_pathway_cci.py`): Create pseudobulk AnnData or use `decoupler` `get_pseudobulk`, then run DE. This respects the experimental unit (donor).

This approach yields reliable, generalizable biomarkers suitable for translational prioritization.

---

## 06 | AI/ML-based biomarker prioritization

**After obtaining ~300 candidate genes from DE + pathway analysis across HSC/mesenchyme, macrophages, ECs, cholangiocytes:**

**Hybrid Prioritization Framework (Rule-based + ML):**

1. **Feature Engineering (per gene, per relevant cell type)**
   - `log2FC_disease` and `abs_log2FC`
   - `-log10(adj_pval)` or rank
   - **Specificity scores**: 
     - Fraction of expressing cells in target compartment vs others (binary or continuous).
     - Tau specificity index or Jensen-Shannon divergence across cell types.
   - **Expression level**: Mean expression in positive cells, detection rate.
   - **Pathway / Mechanism score**: Membership in key fibrosis pathways (TGFB, PDGF, NOTCH, ECM, inflammatory) from MSigDB, Reactome, or decoupler PROGENy / DoRothEA.
   - **Druggability / Tractability**:
     - Surface/secreted/transmembrane vs intracellular (from UniProt, COMPARTMENTS, or `omnipath`).
     - DGIdb, ChEMBL, Open Targets, or CanSAR for known ligands/antibodies/small molecules.
     - "Druggable genome" overlap.
   - **Literature / Prior Knowledge**: PubMed co-occurrence with "liver fibrosis" + cell type, or knowledge graph embedding (e.g., from SPOKE or custom).
   - **Cross-dataset validation**: Consistent DE in optional datasets (GSE244832 etc.) or larger atlases.
   - **Mouse ortholog** and validation status in animal models.
   - **Clinical correlation**: If bulk or other data available, correlation with fibrosis stage, MELD, etc.

2. **Scoring Methods**
   - **Rule-based Composite Score** (transparent, preferred for initial ranking):
     ```python
     score = (0.25 * normalized_specificity +
              0.20 * normalized_|log2FC| +
              0.15 * normalized_-logp +
              0.15 * pathway_relevance +
              0.15 * druggability_score +
              0.10 * lit_support)
     ```
     Normalize each 0–1. Weight tunable by stakeholder (e.g., more weight on druggability for therapeutic targets).

   - **ML Ranking / Classification**:
     - If have a small set of "positive" known good biomarkers/targets, train a supervised model (Random Forest, XGBoost, or even Logistic Regression with regularization) on the feature matrix to predict "high translational value".
     - Use SHAP values for interpretability of why a gene ranks high.
     - Unsupervised: Embed genes (e.g., via autoencoder or graph neural net on PPI + expression) and rank by proximity to known seeds or cluster in "actionable" space.
     - Semi-supervised or active learning if feedback loop with experts.

3. **Final Shortlist & Categorization**
   - Rank all candidates.
   - Categorize: 
     - **Diagnostic biomarkers** (high specificity, secreted/sheddable, detectable in blood/liver biopsy molecular).
     - **Therapeutic targets** (druggable surface, clear mechanism in SAMs or aHSCs or scar-EC, safety profile).
     - **Both / Mechanistic**.
   - Apply filters: e.g., exclude overly broad ECM genes (COL1A1) unless highly specific; prioritize novel or understudied with strong signal.
   - Manual curation + literature review on top 20–30.
   - Output: Ranked table with scores, rationale, cell-type context, and recommended validation experiments.

This hybrid approach is transparent (rule-based) yet leverages ML for complex feature interactions, producing a prioritized list ready for Karyon Bio's translational pipeline.

---

## 07 | Cell-cell interaction and pathway mechanism discovery

**Analysis Approach for Fibrosis Progression (SAMs, activated HSCs, Endothelial remodeling):**

1. **Tools**
   - **Primary**: `LIANA` (Python) — consensus of multiple methods (CellPhoneDB, CellChat, Connectome, etc.) for robustness.
   - Alternatives/Complement: `CellChat`, `NicheNet` (ligand → target TF), `Squidpy` (spatial), `CellPhoneDB` v5+.
   - `omnipath` + `decoupler` for downstream signaling.

2. **Workflow**
   - Subset to key compartments: SAMs/macrophages, activated mesenchyme (aHSC/myofibroblast), scar-associated ECs, and relevant others (cholangiocytes, remaining immune).
   - Run CCI inference on integrated object (or per major batch/disease if needed).
   - Focus on **disease-enriched interactions**: Compare healthy vs cirrhotic (differential CCI) or interactions enriched in fibrotic samples.
   - Key expected axes (from Ramachandran paper):
     - SAMs ↔ PDGFRα+ mesenchyme (PDGF ligands/receptors, TNFSF12-TNFRSF12A/TWEAK pathway).
     - Endothelial (ACKR1+/PLVAP+) ↔ immune cell recruitment and mesenchymal activation.
     - TGF-β, NOTCH, inflammatory loops within the scar niche.

3. **Preventing Overinterpretation**
   - **Statistical rigor**: Use permutation-based specificity tests (as in original CellPhoneDB and paper). Prioritize interactions significant after multiple testing and specific to the interacting pair.
   - **Multi-tool consensus**: Only trust interactions called by ≥2–3 independent methods in LIANA.
   - **Expression validation**: Both ligand **and** receptor must be detectably expressed in the respective clusters (not just inferred).
   - **Biological context**: Cross-reference with pathway activity (PROGENy, DoRothEA) and downstream TF targets (NicheNet).
   - **Spatial co-localization**: If spatial data available (or orthogonal), confirm sender-receiver pairs are in proximity within scars.
   - **Causality disclaimer**: CCI predictions are **hypotheses**. Explicitly state "predicted interaction" and propose validation ( Perturb-seq, neutralizing antibodies, organoid co-cultures, spatial proteomics).
   - **Avoid hype**: Do not claim "drives fibrosis" without functional evidence. Frame as "supports pro-fibrogenic niche remodeling".
   - **Sensitivity**: Run with different clustering resolutions or donor subsets; report robustness.
   - **Negative controls**: Include known non-interacting pairs or shuffled labels.

4. **Output**
   - Ranked interaction tables + heatmaps (sender vs receiver, ligand-receptor pairs).
   - Focused mechanistic hypotheses around SAM → myofibroblast and scar-EC → leukocyte transmigration axes.
   - Integration with biomarker prioritization (e.g., if a prioritized target participates in key CCI, it gains score).

This yields actionable mechanistic insights while maintaining scientific rigor.

---

## 08 | Reproducible pipeline, GitHub workflow, and delivery plan

**12–16 Week End-to-End Delivery Plan (Compact but Rigorous)**

**Repository Structure (Modern Best Practice)**
- `data/` (raw/processed, gitignored or DVC-tracked)
- `notebooks/` (exploratory; converted to scripts via jupytext)
- `scripts/` or `src/karyon_liver_pipeline/` (modular, tested functions)
- `results/` (versioned figures, tables, HTML reports)
- `reports/` (executive summary, methods, biomarker dossier)
- `tests/` (pytest for core functions)
- `envs/`, `Dockerfile`, `environment.yml`
- `.github/workflows/` (CI: lint, test, build docs)
- `README.md`, `CONTRIBUTING.md`, `LICENSE`

**Tools Stack**
- **Core**: Python (scanpy, scvi-tools, anndata, pandas) + LIANA/decoupler/omnipath
- **DE**: pyDESeq2 or R via rpy2 (DESeq2/edgeR)
- **Workflow**: Snakemake or Nextflow (or linear scripts + Makefile)
- **Versioning**: Git + GitHub; DVC or git-lfs for large h5ad
- **Docs**: Sphinx or MkDocs; Jupyter Book for reports
- **CI/CD**: GitHub Actions (black/isort lint, pytest, env build)

**Milestones & Timeline (example for 14 weeks)**

| Weeks | Milestone | Key Activities | Deliverable / QC |
|-------|-----------|----------------|------------------|
| 1-2   | Setup & Curation | Data download, metadata harmonization, samplesheet, initial exploration | Curated metadata + QC report draft |
| 3-4   | QC + Preprocessing | Per-sample QC, doublet/ambient correction, initial integration testing | QC summary figures + filtered AnnData |
| 5-6   | Integration + Annotation | scVI/Harmony integration, automated + manual annotation, validation of key compartments (HSC, Mac, EC) | Annotated object + UMAPs + marker validation report |
| 7-8   | DE + Pathway + CCI | Pseudobulk DE (donor-aware), GSEA/ORA, LIANA CCI focused on fibrosis niche | DE tables, pathway results, CCI hypotheses |
| 9-10  | Biomarker Prioritization | Feature engineering, rule-based + ML scoring, ranking, literature curation | Ranked biomarker table (15–20) + scoring methodology |
| 11-12 | Validation & Polish | Optional dataset checks, sensitivity analyses, figure polishing, robustness | Final figures, supplementary tables |
| 13-14 | Documentation & Delivery | Code clean-up, tests, README, Executive Summary (1-2 pg), full report, packaging | GitHub release/tag + zipped folder + submission package |

**Quality Checks Throughout**
- Visual QC at every major step (UMAPs, violin plots, heatmaps) saved automatically.
- Session info + environment export at each milestone.
- Code review (self or peer) on key scripts.
- Reproducibility test: Fresh clone + `environment.yml` → run core pipeline end-to-end on subsample.
- Provenance: Track random seeds, parameters, git commit in output metadata.
- Automated testing for critical functions (e.g., pseudobulk creation, scoring).

**Final Deliverables (as required)**
- GitHub repo (or zipped folder)
- README with setup + high-level results
- QC summary
- Cell annotation figures + validation
- DE + pathway results
- Ranked biomarker/target table
- 1–2 page Executive Summary for Karyon Bio
- This written answers document
- Full methods transparency

**Risk Mitigation**
- Scope creep: Stick to primary dataset + focused compartments.
- Compute: Design scripts to be memory-efficient; provide subsample mode.
- Biology vs tech: Multiple validation layers as described in Q2–Q4.

This plan delivers high-quality, translational-ready outputs on time while maintaining scientific integrity and reproducibility — exactly what is needed for Karyon Bio's single-cell framework.

---

**End of Written Responses**

These answers demonstrate both technical depth in modern single-cell methods and translational judgment tailored to liver fibrosis biomarker discovery. Ready for practical project execution or discussion.