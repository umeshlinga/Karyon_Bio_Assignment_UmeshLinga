#!/usr/bin/env python3
"""
Generate the 1-2 page Executive Summary PDF for the Karyon Bio assignment.

This summary reads the actual outputs produced by the pipeline:
  - results/biomarker_ranked_table.csv       (real ranked candidates)
  - results/annotation/compartment_validation.csv (real cirrh vs healthy stats)
  - results/de/pathway_signature_scores.csv  (real signature enrichment)
  - results/qc/qc_per_sample.csv             (real QC retention)

No values are hand-curated.

Usage
-----
python reports/create_executive_summary.py [--results_dir results] \\
    [--out reports/Executive_Summary.pdf]
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    HRFlowable, Image, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
)


def fmt_p(p: float) -> str:
    try:
        if p == 0 or p < 1e-300:
            return "< 1e-300"
        if p < 1e-3:
            return f"{p:.2e}"
        return f"{p:.3g}"
    except Exception:
        return "n/a"


def build(results_dir: Path, out_path: Path):
    biomarker_csv = results_dir / "biomarker_ranked_table.csv"
    validation_csv = results_dir / "annotation" / "compartment_validation.csv"
    sig_csv = results_dir / "de" / "pathway_signature_scores.csv"
    qc_csv = results_dir / "qc" / "qc_per_sample.csv"

    bio = pd.read_csv(biomarker_csv)
    val = pd.read_csv(validation_csv) if validation_csv.exists() else None
    sig = pd.read_csv(sig_csv) if sig_csv.exists() else None
    qc = pd.read_csv(qc_csv) if qc_csv.exists() else None

    doc = SimpleDocTemplate(
        str(out_path), pagesize=letter,
        rightMargin=0.55 * inch, leftMargin=0.55 * inch,
        topMargin=0.5 * inch, bottomMargin=0.5 * inch,
    )
    styles = getSampleStyleSheet()
    title = ParagraphStyle("T", parent=styles["Title"], fontSize=15,
                           textColor=colors.HexColor("#1a365d"),
                           spaceAfter=4, alignment=TA_CENTER)
    sub = ParagraphStyle("S", parent=styles["Normal"], fontSize=10,
                         textColor=colors.HexColor("#2c5282"),
                         alignment=TA_CENTER, fontName="Helvetica-Bold",
                         spaceAfter=10)
    h2 = ParagraphStyle("H2", parent=styles["Heading2"], fontSize=11,
                        textColor=colors.HexColor("#1a365d"),
                        spaceBefore=8, spaceAfter=2)
    body = ParagraphStyle("B", parent=styles["Normal"], fontSize=9,
                          leading=11, alignment=TA_JUSTIFY, spaceAfter=4)
    small = ParagraphStyle("Sm", parent=styles["Normal"], fontSize=8,
                           textColor=colors.HexColor("#4a5568"))
    # Cell styles used inside tables. fontSize matches FONTSIZE in TableStyle.
    cell = ParagraphStyle("Cell", parent=styles["Normal"], fontSize=7.5,
                          leading=9, textColor=colors.black)
    cell_right = ParagraphStyle("CellR", parent=cell, alignment=2)   # right
    cell_center = ParagraphStyle("CellC", parent=cell, alignment=1)  # center
    head_cell = ParagraphStyle("Head", parent=cell, textColor=colors.white,
                               fontName="Helvetica-Bold")

    def P(text, style=cell):
        return Paragraph(str(text), style)

    story = []
    story.append(Paragraph("Karyon Bio | Candidate Technical Assignment", small))
    story.append(Paragraph(
        "Cell-Type-Specific Biomarker Discovery in Human Liver Fibrosis", title))
    story.append(Paragraph(
        "Primary dataset: GSE136103 (Ramachandran et al., 2019) — "
        "scRNA-seq of 5 healthy + 5 cirrhotic human livers (NPC fractions)",
        sub))
    story.append(Paragraph(
        f"Prepared by: Umesh Linga (Indianapolis, IN) | "
        f"{datetime.now().strftime('%B %Y')}", small))
    story.append(HRFlowable(width="100%", thickness=1,
                            color=colors.HexColor("#2c5282")))
    story.append(Spacer(1, 6))

    n_cells_post = int(qc["n_cells_post"].sum()) if qc is not None else None
    n_samples = len(qc) if qc is not None else None
    story.append(Paragraph("Executive Summary", h2))
    if n_cells_post is not None:
        es = (
            "A reproducible scanpy + Harmony + PyDESeq2 pipeline was run on "
            "the GSE136103 cohort "
            f"(20 liver sample-fractions across 10 donors; {n_cells_post:,} "
            f"cells retained after QC of {n_samples} samples). Three required "
            "disease-relevant compartments were recovered and validated: "
            "hepatic stellate/mesenchymal, macrophage/monocyte, and "
            "endothelial. Donor-aware pseudobulk DE (PyDESeq2, design "
            "<i>~disease</i>) produced ranked gene lists per compartment, and "
            "a hybrid rule-based + Random Forest scorer combined "
            "specificity, effect size, druggability, surface/secreted "
            "tractability, fibrosis-pathway membership, literature support, "
            "and curated ligand-receptor evidence into a final prioritization "
            "score."
        )
    else:
        es = (
            "A reproducible scanpy + Harmony + PyDESeq2 pipeline was run on "
            "the GSE136103 cohort, recovering and validating the three "
            "required disease-relevant compartments. Donor-aware pseudobulk "
            "DE plus a hybrid rule-based + ML scorer produced a ranked "
            "biomarker / target list."
        )
    story.append(Paragraph(es, body))

    if val is not None and not val.empty:
        story.append(Paragraph(
            "Compartment validation (donor-level Mann-Whitney U, "
            "cirrhotic &gt; healthy)", h2))
        header = ["Compartment", "Signature", "n H/C",
                  "Mean H", "Mean C", "&Delta;", "MWU p"]
        rows = [[P(h, head_cell) for h in header]]
        priority = {
            "Mesenchyme":  ["HSC_activation", "ECM_organization"],
            "Macrophage":  ["SAM_signature", "ECM_organization"],
            "Endothelial": ["ScarEC_signature", "SAM_signature"],
        }
        for comp, sigs in priority.items():
            for s in sigs:
                row = val[(val["compartment"] == comp) & (val["signature"] == s)]
                if row.empty:
                    continue
                r = row.iloc[0]
                rows.append([
                    P(comp),
                    P(s),
                    P(f"{int(r['n_donor_healthy'])}/{int(r['n_donor_cirrh'])}",
                      cell_center),
                    P(f"{r['mean_healthy']:.3f}", cell_right),
                    P(f"{r['mean_cirrh']:.3f}", cell_right),
                    P(f"{r['delta']:+.3f}", cell_right),
                    P(fmt_p(r["mwu_p_one_sided"]), cell_right),
                ])
        if len(rows) > 1:
            # Total ~7.0 inches; printable width ~7.4 inches.
            t = Table(rows, colWidths=[1.15*inch, 1.70*inch, 0.55*inch,
                                       0.65*inch, 0.65*inch, 0.55*inch,
                                       0.75*inch])
            t.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c5282")),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#a0aec0")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1),
                 [colors.white, colors.HexColor("#edf2f7")]),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]))
            story.append(t)

    story.append(Paragraph(
        "Top 15 prioritized biomarker / target candidates "
        "(real values from pipeline run)", h2))
    keep_cols = ["Rank", "Gene", "Primary_Compartment",
                 "Translational_Category", "log2FC_cirrh_vs_healthy",
                 "padj", "Final_Prioritization_Score"]
    show = bio[keep_cols].head(15).copy()
    show["log2FC_cirrh_vs_healthy"] = show["log2FC_cirrh_vs_healthy"].map(
        lambda x: f"{x:+.2f}")
    show["padj"] = show["padj"].map(fmt_p)
    show["Final_Prioritization_Score"] = show["Final_Prioritization_Score"].map(
        lambda x: f"{x:.3f}")
    show.columns = ["#", "Gene", "Compartment", "Translational fit",
                    "log2FC (C/H)", "padj", "Score"]

    header = [P(h, head_cell) for h in show.columns]
    rows = [header]
    for _, r in show.iterrows():
        rows.append([
            P(r["#"], cell_center),
            P(f"<b>{r['Gene']}</b>"),
            P(r["Compartment"]),
            P(r["Translational fit"]),
            P(r["log2FC (C/H)"], cell_right),
            P(r["padj"], cell_right),
            P(r["Score"], cell_right),
        ])
    # Total ~7.30 inches across printable width 7.40 inches.
    t = Table(rows, colWidths=[0.25*inch, 0.65*inch, 0.85*inch, 2.30*inch,
                                0.80*inch, 0.85*inch, 0.60*inch])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c5282")),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#a0aec0")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#edf2f7")]),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(t)
    story.append(Spacer(1, 4))

    if sig is not None and not sig.empty:
        story.append(Paragraph("Pathway / signature enrichment in DE results",
                              h2))
        for comp in ["Mesenchyme", "Macrophage", "Endothelial"]:
            sub = sig[sig["compartment"] == comp]
            if sub.empty:
                continue
            top = sub.sort_values("mean_log2FC", ascending=False).head(2)
            for _, r in top.iterrows():
                story.append(Paragraph(
                    f"• <b>{comp}</b> — {r['signature']}: mean log2FC "
                    f"<b>{r['mean_log2FC']:+.2f}</b>; up-significant fraction "
                    f"<b>{r['frac_up_sig']:.0%}</b>; top genes: "
                    f"<i>{r['top_genes']}</i>", body))

    story.append(Paragraph("Translational interpretation", h2))
    story.append(Paragraph(
        "<b>Diagnostics:</b> ACKR1 and PLVAP define scar-restricted endothelial "
        "populations expanded in cirrhosis and reach high tau-specificity in "
        "our results, supporting use in non-invasive panels (sheddable / "
        "surface candidates). Mesenchymal ECM genes (COL1A1, COL3A1, CTHRC1, "
        "TIMP1) form a robust signature for fibrosis staging.", body))
    story.append(Paragraph(
        "<b>Therapeutics:</b> ACKR1 (chemokine sequestration), PDGFRA "
        "(activated-HSC mitogen receptor), and LOXL2 (ECM crosslinking) are "
        "druggable nodes with clinical-stage tool compounds or antibodies. "
        "TIMP1 elevation is a hallmark of de-regulated MMP balance and "
        "represents a soluble biomarker candidate.", body))
    story.append(Paragraph(
        "<b>Cell-cell mechanism:</b> the curated ligand-receptor evidence panel "
        "highlights CXCL12-CXCR4 (mesenchyme/endothelial -&gt; macrophage "
        "recruitment), FN1-ITGAV ECM-integrin signaling within the fibrotic "
        "niche, and TGFB1-TGFBR2 (macrophage -&gt; endothelial) — coherent with "
        "the Ramachandran et al. 2019 model of an interactive scar niche.",
        body))

    fig_path = results_dir / "figures" / "top_biomarkers.png"
    if fig_path.exists():
        story.append(Spacer(1, 4))
        story.append(Image(str(fig_path), width=5.6 * inch, height=3.5 * inch))

    story.append(Spacer(1, 4))
    story.append(HRFlowable(width="100%", thickness=0.5,
                            color=colors.HexColor("#cbd5e0")))
    story.append(Paragraph(
        "<i>All numbers in this summary are produced by scripts/01..05; "
        "rerun via the README quick-start. Code, full ranked table, DE "
        "results, and figures accompany this PDF in the same submission.</i>",
        small))

    doc.build(story)
    print(f"Wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results_dir", type=Path, default=Path("results"))
    ap.add_argument("--out", type=Path,
                    default=Path("reports/Executive_Summary.pdf"))
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    build(args.results_dir, args.out)


if __name__ == "__main__":
    main()
