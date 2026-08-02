#!/usr/bin/env python3
"""
generate_case_study_pdf.py — one-page case-study PDF (BUILD_PLAN.md §13 deliverables checklist).

Numbers here are transcribed from the real results files (data/processed/evaluation/results.tsv,
data/processed/recalibration/results.tsv, data/processed/pgs_catalog/stageB3_descriptive_comparison.tsv)
as of the run documented in METHODS.md — re-run and update by hand if those results change; this
script doesn't read them live, since a one-page summary should present a fixed, reviewed set of
headline numbers rather than silently drifting with every pipeline re-run.

Usage:
    python3 reports/final_writeup/generate_case_study_pdf.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

OUT_PATH = Path(__file__).resolve().parent / "case_study_summary.pdf"

fig = plt.figure(figsize=(8.5, 11))
fig.patch.set_facecolor("white")

y = 0.97

def text(s, size=10, weight="normal", color="black", dy=0.022, x=0.06, style="normal"):
    global y
    fig.text(x, y, s, fontsize=size, weight=weight, color=color, style=style, va="top", wrap=True)
    y -= dy

text("PRS Cross-Ancestry Portability & Recalibration Benchmark", size=15, weight="bold", dy=0.03)
text("Case Study Summary", size=11, color="#444444", dy=0.035)

text("THE PROBLEM", size=10, weight="bold", color="#1a4d8f", dy=0.024)
text("Polygenic risk scores (PRS) trained on predominantly European-ancestry GWAS are well", size=9.3)
text("documented to underperform when applied to other ancestries. This project measures that", size=9.3)
text("gap with a controlled, synthetic ground truth on real 1000 Genomes multi-ancestry", size=9.3)
text("genotypes (GCTA-simulated phenotype, known causal variants and heritability) — and", size=9.3)
text("checks the same pattern against a real, published PGS Catalog score.", size=9.3, dy=0.032)

text("HEADLINE FINDING", size=10, weight="bold", color="#1a4d8f", dy=0.024)
text("A simple clumping+thresholding PRS (PRSice-2), tuned to its best EUR threshold, beats both", size=9.3)
text("Bayesian genome-wide methods (PRS-CSx, LDpred2-auto) in the EUR held-out set — but loses to", size=9.3)
text("PRS-CSx in every single one of 8 non-EUR (ancestry × scenario) comparisons, and to", size=9.3)
text("LDpred2-auto in 6 of 8. No method wins universally (see Scope & Limitations):", size=9.3, dy=0.026)

# table
col_x = [0.08, 0.24, 0.38, 0.56, 0.74]
headers = ["Scenario", "Ancestry", "PRSice-2 (best)", "PRS-CSx", "LDpred2-auto"]
rows = [
    ["1 (equal)", "AFR", "0.0077", "0.0275", "0.0117"],
    ["1", "AMR", "0.0098", "0.0594", "0.0136"],
    ["1", "EAS", "0.0008", "0.0165", "0.0015"],
    ["1", "SAS", "0.0080", "0.0279", "0.0118"],
    ["2 (varying)", "AFR", "0.0081", "0.0230", "0.0130"],
    ["2", "AMR", "0.0009", "0.0480", "0.0065"],
    ["2", "EAS", "0.0090", "0.0238", "0.0251"],
    ["2", "SAS", "0.0161", "0.0411", "0.0112"],
]
for cx, h in zip(col_x, headers):
    fig.text(cx, y, h, fontsize=8.3, weight="bold", va="top")
y -= 0.017
for r in rows:
    for cx, val in zip(col_x, r):
        fig.text(cx, y, val, fontsize=8.3, va="top", family="monospace")
    y -= 0.0155
y -= 0.01
fig.text(0.06, y, "R² against the true simulated phenotype. PRS-CSx pattern holds in all 8/8 non-EUR comparisons; LDpred2-auto in 6/8.",
          fontsize=8, style="italic", color="#555555", va="top")
y -= 0.032

text("RECALIBRATION: WHAT IT FIXES AND WHAT IT DOESN'T", size=10, weight="bold", color="#1a4d8f", dy=0.024)
text("Empirical per-ancestry recentering/rescaling recovers 100% of calibration (mean/scale", size=9.3)
text("matched to the EUR reference in every case) and 0% of discriminative accuracy (R²", size=9.3)
text("identical before/after to 4 decimal places, in all 95 result rows) — by mathematical", size=9.3)
text("necessity, since linear recentering cannot change a Pearson correlation. This makes a", size=9.3)
text("distinction the field warns is often conflated concrete rather than asserted.", size=9.3, dy=0.032)

text("REAL-SCORE CROSS-CHECK (DESCRIPTIVE ONLY)", size=10, weight="bold", color="#1a4d8f", dy=0.024)
text("A real published height PGS (PGS000297) applied to the same genotypes: 0 of 5", size=9.3)
text("super-populations' score ranking matches published national height statistics — a real,", size=9.3)
text("small-scale illustration of the same portability problem, not a validated accuracy claim", size=9.3)
text("(1000 Genomes has no phenotype data to validate against).", size=9.3, dy=0.032)

text("SCOPE & LIMITATIONS", size=10, weight="bold", color="#1a4d8f", dy=0.024)
text("chr21+chr22 smoke-test scope (not genome-wide); LDpred2-auto's LD reference is chr22-only,", size=8.6, color="#333333")
text("narrower than PRS-CSx's official multi-population panel, a plausible driver of its weaker", size=8.6, color="#333333")
text("cross-ancestry showing above; single simulation replicate; small per-ancestry sample sizes.", size=8.6, color="#333333")
text("Full detail: METHODS.md and reports/final_writeup/TECHNICAL_WRITEUP.md.", size=8.6, color="#333333", dy=0.03)

fig.text(0.06, 0.03, "PRS Cross-Ancestry Portability & Recalibration Benchmark  —  github.com/scholaryzeacademy/prs-cross-ancestry-portability-benchmark",
          fontsize=7.5, color="#888888", va="bottom")

with PdfPages(OUT_PATH) as pdf:
    pdf.savefig(fig)
plt.close(fig)
print(f"wrote {OUT_PATH}")
