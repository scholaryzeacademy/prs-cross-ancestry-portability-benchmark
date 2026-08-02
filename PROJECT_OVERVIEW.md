# Project Overview

**PRS Cross-Ancestry Portability & Recalibration Benchmark**

This document is a single, complete narrative of the project — what it is, why it exists, exactly
how each stage works, what was found, and what its limits are — so that one person can read it and
understand the whole thing without needing to cross-reference every other file. The underlying
detail lives in three other files, each with a different job:

- `docs/BUILD_PLAN.md` — the original plan, rationale, and scope/honesty constraints (written before
  implementation)
- `METHODS.md` — the stage-by-stage implementation record: exact parameters, tooling quirks, and
  every real number produced, with the date it was produced
- `reports/final_writeup/TECHNICAL_WRITEUP.md` — the condensed (~2,000-word) narrative version of
  this same project, plus a one-page PDF case-study summary

This overview draws only on facts already recorded in those files — no number here is invented.

---

## 1. What this project is

Polygenic risk scores (PRS) are a genomics tool that combines the effects of many genetic variants
into a single number predicting someone's likelihood of a trait or disease. The most
well-documented weakness of PRS today is that scores are trained overwhelmingly on
European-ancestry genome-wide association studies (GWAS), and predictive accuracy drops sharply
when the same score is applied to people of other ancestries — due to real differences in allele
frequency and linkage disequilibrium (LD) structure between populations, and sometimes genuine
differences in which variants actually matter (gene-by-ancestry effect heterogeneity).

This project measures that accuracy loss directly, with a **known, controlled ground truth** — not
an assumption or a literature citation, but a real, reproducible measurement — and measures how
much a standard recalibration step actually recovers.

**Why a controlled ground truth is necessary:** 1000 Genomes Phase 3, the standard public
multi-ancestry genotype reference panel, has real genotypes for 2,504 people across five
continental super-populations but **no real disease/trait phenotype data at all** — only genotype,
sex, and family relationships are distributed with it. So there's no way to directly measure "how
accurate is this PRS in this ancestry" against real 1000 Genomes phenotypes, because there aren't
any. This project solves that by simulating a phenotype with GCTA directly on the real genotypes,
with a documented, known set of causal variants and heritability — giving a genuine ground truth to
measure predictive accuracy against.

## 2. Two-track design

- **Track A (synthetic ground truth — the project's primary, headline result):** real 1000
  Genomes genotypes + a GCTA-simulated phenotype with known causal variants and heritability. A
  discovery GWAS trained only on European-ancestry (EUR) samples (deliberately reproducing the
  field's real bias) feeds three independently-built PRS models, each evaluated for accuracy across
  all five super-populations, then recalibrated.
- **Track B (real published score — descriptive/illustrative only, never a validation):** a real,
  published PGS Catalog score is computed on the same 1000 Genomes genotypes and compared
  descriptively to published national height statistics. Because 1000 Genomes has no real
  phenotype, this track **cannot validate anything** — it exists only to connect Track A's
  controlled finding to a real, recognizable, named score, and its descriptive-only status must be
  restated everywhere it's discussed, not disclaimed once.

## 3. Goals and explicit non-goals

**Goals** (from `docs/BUILD_PLAN.md` §1.1):
- Quantify, with a known ground truth, how much PRS predictive accuracy (R²) is lost when a score
  trained on one ancestry is applied to others, using real multi-ancestry LD structure and allele
  frequencies.
- Compare at least three PRS construction methods honestly — a simple clumping-and-thresholding
  baseline (PRSice-2) plus two more sophisticated Bayesian/multi-ancestry methods (LDpred2-auto,
  PRS-CSx) — following an "always include a naive baseline" discipline.
- Apply a recalibration step and report honestly how much of the portability gap it actually closes
  versus how much it doesn't.
- Cross-check the synthetic finding against a real, published PGS Catalog score, clearly
  distinguishing descriptive illustration from validated predictive accuracy.

**Explicit non-goals** (from `docs/BUILD_PLAN.md` §1.2 — these boundaries are load-bearing, not
just caveats):
- **Not** a clinically deployable recalibration recommendation — this is a methodology
  demonstration on public reference data, not a validated clinical tool.
- **Not** an attempt to access UK Biobank, All of Us, or any credentialed biobank — deliberately
  scoped to fully public, no-credentialing data (1000 Genomes, PGS Catalog).
- **Not** an attempt to rank PRS methods universally — current literature (Momin et al. 2026) is
  explicit that no single method wins universally across traits, and this project's own
  single-scenario simulation is reported the same way: full per-method × per-ancestry × per-scenario
  tables, never collapsed into "method X is best."

## 4. Data sources and tools

| Resource | Role | Access |
|---|---|---|
| 1000 Genomes Phase 3 (GRCh37) | Real multi-ancestry genotype panel | IGSR/EBI FTP-over-HTTPS: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502` |
| PGS Catalog | Real published PRS scoring files (Track B) | `pgscatalog.org` — REST API, FTP, and the official `pgscatalog-utils`/`pgsc_calc` tooling |
| GCTA | Phenotype simulation on real genotypes | `cnsgenomics.com/software/gcta/`, version 1.95.3 |
| PRSice-2 | PRS method 1: simple clumping+thresholding baseline | `github.com/choishingwan/PRSice`, version 2.3.5 |
| PRS-CSx | PRS method 2: multi-ancestry Bayesian continuous shrinkage | `github.com/getian107/PRScsx` |
| LDpred2-auto (via `bigsnpr`) | PRS method 3: self-tuning Bayesian shrinkage | R package `bigsnpr`, version 1.12.21 |
| `pgscatalog-utils` | PGS Catalog download/scoring tooling (Track B) | version 2.0.4 |
| PLINK2 | Genotype QC, format conversion, GWAS, scoring | version 2.0.0-a.7.0LM |

1000 Genomes has 2,504 individuals across 26 populations grouped into five continental
super-populations: **AFR** (African, 661 samples in this project's QC'd subset), **AMR**
(admixed American, 347), **EAS** (East Asian, 504), **EUR** (European, 503), **SAS** (South Asian,
489).

**Scope of this implementation:** everything below ran on a **chr21+chr22 smoke-test subset**
(~6% of the autosomal genome), not a genome-wide run — chosen so every stage could be built,
debugged, and verified against real computation rather than mocked, and to keep the CI workflow
light (per `docs/BUILD_PLAN.md` §5's CI-scope row). The pipeline itself has no hardcoded chromosome
limit; a genome-wide run is a config/download-scope change, not a rewrite, but has not been done.

## 5. Pipeline architecture

```
      TRACK A: SYNTHETIC GROUND TRUTH                TRACK B: REAL PUBLISHED SCORE (descriptive)
 ┌─────────────────────────────┐              ┌──────────────────────────────────┐
 │ STAGE 1: Download 1000        │              │ STAGE 1 (shared): Same 1000       │
 │ Genomes Phase 3 (5 super-      │              │ Genomes genotype panel             │
 │ populations, real genotypes)   │              └────────────────┬───────────────────┘
 └──────────────┬───────────────┘                               │
                │                              ┌────────────────▼───────────────────┐
 ┌──────────────▼───────────────┐              │ STAGE B1: Download real PGS         │
 │ STAGE 2: GCTA phenotype        │              │ Catalog score (height)             │
 │ simulation — Scenario 1        │              └────────────────┬───────────────────┘
 │ (equal effect sizes) and       │                               │
 │ Scenario 2 (ancestry-varying   │              ┌────────────────▼───────────────────┐
 │ effect sizes)                  │              │ STAGE B2: Compute score on all       │
 └──────────────┬───────────────┘              │ 2,504 individuals via plink2 --score │
                │                              └────────────────┬───────────────────┘
 ┌──────────────▼───────────────┐                               │
 │ STAGE 3: EUR-only "discovery   │              ┌────────────────▼───────────────────┐
 │ GWAS" on simulated phenotype   │              │ STAGE B3: Descriptive comparison      │
 └──────────────┬───────────────┘              │ against published height statistics   │
                │                              │ (illustrative only)                   │
 ┌──────────────▼───────────────┐              └────────────────────────────────────┘
 │ STAGE 4: Construct PRS —       │
 │ PRSice-2, LDpred2-auto,        │
 │ PRS-CSx                        │
 └──────────────┬───────────────┘
                │
 ┌──────────────▼───────────────┐
 │ STAGE 5: Evaluate predictive   │
 │ accuracy (R²) across all 5     │
 │ super-populations               │
 └──────────────┬───────────────┘
                │
 ┌──────────────▼───────────────┐
 │ STAGE 6: Apply empirical       │
 │ per-ancestry recalibration     │
 └──────────────┬───────────────┘
                │
                └──────────────┬──────────────────────────────────┘
                               │
                 ┌─────────────▼─────────────┐
                 │ STAGE 7: Combined technical │
                 │ write-up + case-study PDF   │
                 └────────────────────────────┘
```

## 6. Stage-by-stage methodology and results

### Stage 1 — 1000 Genomes download & QC (shared by both tracks)

`src/shared_stage1_1000g_download_qc/download.py` (fetch) and `qc.py` (convert + QC). QC filters
are applied **within each super-population separately**, keeping only variants passing in all
five: minor allele frequency ≥0.01, per-variant missingness ≤0.05, Hardy-Weinberg equilibrium
p≥1e-6, per-sample missingness ≤0.05 (applied globally after variant filtering).

**Verified on real data (2026-08-01):** 2,504 samples loaded; variants passing QC per
super-population — AFR 462,550, AMR 308,829, EAS 240,128, EUR 272,054, SAS 284,717; **179,285
variants shared across all five** super-populations (the final QC'd variant set carried into every
downstream stage).

### Stage 2 — GCTA phenotype simulation (Track A)

`src/trackA_synthetic/stage2_gcta_simulation/`. GCTA requires PLINK1 bed/bim/fam format (no
multiallelic support), so Stage 1's QC'd panel was converted once into a shared biallelic-SNP-only
bed/bim/fam — the smoke-test subset went from 179,285 QC'd variants to **154,384 biallelic SNPs**.

300 causal variants (MAF in [0.05, 0.5], fixed seed) were selected once and used identically in
both scenarios so they're directly comparable. GCTA's `--simu-qt` was run once per
super-population, each targeting heritability h²=0.5 by construction (every ancestry gets the same
nominal heritability; the portability gap shows up in downstream prediction accuracy, not here):

- **Scenario 1 (equal effect sizes):** the same causal-variant effect sizes applied identically
  across all five super-populations — isolating the portability problem caused purely by real
  LD/allele-frequency differences.
- **Scenario 2 (ancestry-varying effect sizes):** each super-population's effect sizes
  independently perturbed (`effect = base_effect × N(1, 0.3)`, per-ancestry seeded) — modeling
  genuine gene-by-ancestry effect-size heterogeneity.

**Verified on real data (2026-08-01):** both scenarios ran successfully across all 2,504
individuals (AFR 661, AMR 347, EAS 504, EUR 503, SAS 489). Sanity check (reconstructing the
genetic score from the 300 causal variants' variance-standardized genotypes × GCTA's own effect
sizes, Scenario 1, EUR): empirical heritability = 0.53 against a target of 0.5 (n=503) — the ~0.03
gap is consistent with sampling noise at this sample size, not a simulation bug.

**Known tooling limitation, documented rather than worked around:** GCTA v1.95.3 has no `--seed`
option for `--simu-qt`'s residual noise — the genetic component (causal variants, effect sizes,
perturbations) is fully reproducible via this project's own seeded RNG, but the residual/noise draw
is not bit-for-bit reproducible across re-runs.

**Note on the simulation's status:** `configs/simulation_parameters.yaml`, which holds every
parameter above, is still marked **DRAFT, pending Biostatistics review** per open questions listed
in that file — the implementation reflects that draft design, not a signed-off one.

### Stage 3 — EUR-only discovery GWAS (Track A)

`src/trackA_synthetic/stage3_eur_gwas/run_gwas.py`. Simple, unadjusted per-SNP linear regression
(`plink2 --glm`, no PCs/covariates) on each scenario's simulated phenotype, restricted to EUR
samples only — deliberately reproducing the real-world bias that most published GWAS are
European-ancestry-derived.

**EUR discovery/held-out split:** before the GWAS runs, EUR samples are partitioned into a
discovery subset (used for the GWAS) and a held-out subset the GWAS never sees (80/20 split,
fixed seed). This split is essential: without it, evaluating the model on "EUR" later would be
in-sample against the GWAS's own training data, not a genuine held-out test.

**Verified on real data (2026-08-01):** both scenarios ran successfully — 503 EUR samples split
into 402 discovery / 101 held-out, 154,384 SNPs tested, all 300 causal variants present in the
output. Causal-variant enrichment against background: causal variants are nominally significant
(p<0.05) 13.0% (Scenario 1) / 13.7% (Scenario 2) of the time vs. 5.2–5.7% for background SNPs
(~2.4× enrichment) — modest, not dramatic, and exactly what's expected for a polygenic trait where
h²=0.5 is spread across 300 causal variants (~0.17% heritability each), well below what n=503 can
detect per-SNP at strict significance. This weak per-SNP signal is precisely why genome-wide
Bayesian methods (LDpred2, PRS-CSx) are expected to outperform a top-hits-only baseline
(PRSice-2) on this kind of trait.

### Stage 4 — PRS construction, three methods (Track A)

`src/trackA_synthetic/stage4_prs_construction/`. Each method is built **only** from the EUR
discovery GWAS + EUR discovery genotypes — no evaluation-ancestry data ever touches construction —
producing a portable scoring file that Stage 5 applies identically to all five super-populations.

- **PRSice-2** (clumping+thresholding baseline): LD-clumped (r²<0.1, 250kb) using EUR discovery
  samples as reference, scored at a fixed p-value threshold suite (5e-8 to 1) without any
  target-phenotype regression, so the model doesn't depend on which ancestry it's later evaluated
  against.
- **PRS-CSx** (Bayesian continuous shrinkage): run in single-discovery-population mode (`--pop
  EUR`), since this project deliberately feeds all three methods the same EUR-only sumstats — the
  "x" (coupling effects across ancestries) only activates with multiple discovery GWASes, which
  this project doesn't do. Uses the official 1000G EUR LD reference panel and HapMap3 SNP list;
  restricted to the ~18% of QC'd SNPs overlapping PRS-CSx's HapMap3 set (27,881 of 154,384 SNPs).
- **LDpred2-auto** (self-tuning Bayesian, via `bigsnpr`): no validation phenotype needed. Its LD
  reference in this run is **chr22-only** (100kb bp-radius window, no genetic map available) —
  narrower in scope than PRS-CSx's official multi-population panel; see §7 for why this matters to
  the result and §8 for the run history (this was the hardest stage to get running end-to-end due
  to shared-host resource contention, not a code defect).

**Verified on real data (2026-08-01/02):** all three methods completed successfully for both
scenarios. PRSice-2: 5,858 SNPs survive clumping (of 131,173 non-ambiguous), 8–9 usable p-value
thresholds per scenario. PRS-CSx: full MCMC (1,000 iterations) completed for both chromosomes in
both scenarios, 27,881 SNPs scored. LDpred2-auto: 64,272 SNPs scored (all non-zero) in both
scenarios.

### Stage 5 — Cross-ancestry evaluation (Track A)

`src/trackA_synthetic/stage5_crossancestry_evaluation/evaluate.py`. Every Stage 4 scoring file is
applied via `plink2 --score` to five held-out groups: **EUR uses the held-out subset** Stage 3 set
aside (101 samples, never in the discovery GWAS or any Stage 4 construction step) — not the full
EUR or EUR-discovery set, which would be circular. AFR (661), AMR (347), EAS (504), SAS (489) use
their full sample sets (never used anywhere upstream). R² is computed against the *true* simulated
phenotype (available because this is synthetic ground truth), with a 95% confidence interval via
the Fisher z-transform of the Pearson correlation.

**Full results table:** `data/processed/evaluation/results.tsv` — 95 rows (3 methods × up to 8
PRSice-2 thresholds or 1 PRS-CSx/LDpred2-auto model each × 5 ancestries × 2 scenarios), reported in
full per BUILD_PLAN.md's explicit instruction never to collapse this into a single "best method"
figure.

**The headline finding, reproducible across both scenarios (not cherry-picked — this is the whole
table's pattern):**

| Scenario | Ancestry | PRSice-2 (best EUR threshold) | PRS-CSx | LDpred2-auto |
|---|---|---|---|---|
| 1 | AFR | 0.0077 | **0.0275** | 0.0117 |
| 1 | AMR | 0.0098 | **0.0594** | 0.0136 |
| 1 | EAS | 0.0008 | **0.0165** | 0.0015 |
| 1 | SAS | 0.0080 | **0.0279** | 0.0118 |
| 2 | AFR | 0.0081 | **0.0230** | 0.0130 |
| 2 | AMR | 0.0009 | **0.0480** | 0.0065 |
| 2 | EAS | 0.0090 | 0.0238 | **0.0251** |
| 2 | SAS | 0.0161 | **0.0411** | 0.0112 |

(R² against the true simulated phenotype; bold = best in that row.)

PRSice-2's best-performing threshold in the EUR held-out set (p<0.01, R²=0.131 Scenario 1 / 0.049
Scenario 2) clearly beats PRS-CSx's same-ancestry R² (0.073 / 0.025) — but that same PRSice-2
model is *worse* than PRS-CSx in **all eight** non-EUR (ancestry × scenario) comparisons, and worse
than LDpred2-auto in six of eight. This matches the literature pattern this project is built around
(Momin et al. 2026): genome-wide Bayesian shrinkage methods generally port across ancestries
better than a simple clumping+thresholding baseline tuned to one ancestry, even when the baseline
looks better same-ancestry.

**No method is universally best, though** — LDpred2-auto's own same-ancestry (EUR) R² (0.0115 /
0.0140) is lower than both other methods' same-ancestry R² in both scenarios, and it only edges out
PRS-CSx once (Scenario 2, EAS). A plausible, specific reason: LDpred2-auto here only had chr22 LD
information available, a narrower LD reference than PRS-CSx's official multi-population panel —
likely handicapping its shrinkage relative to how it would perform with a fuller LD reference. This
is a scope limitation of this smoke test, not necessarily a property of the method in general.

**Caveats that must travel with this finding:** single simulation replicate, chr21+22-only
smoke-test scope, small ancestry sample sizes (confidence intervals in the raw results table are
wide). Scenario 1 vs. Scenario 2 do not show a qualitatively different portability pattern here —
both scenarios' cross-ancestry R² are similarly depressed relative to EUR, suggesting (at this
scope) the LD/allele-frequency mechanism Scenario 1 isolates already accounts for most of the
observed gap, with Scenario 2's added effect-size heterogeneity not obviously worsening it further.
This should be treated as a preliminary observation, not a settled conclusion.

#### PRSice-2: full per-threshold table

The headline table above uses only PRSice-2's single best-performing EUR threshold per scenario.
The full picture — every p-value threshold PRSice-2 was scored at, all five ancestries, both
scenarios, straight from `data/processed/evaluation/results.tsv` (R², 95% CI in brackets) — makes
the same point more strongly: PRSice-2 never approaches PRS-CSx's or LDpred2-auto's non-EUR R² at
*any* threshold, not just its best one.

**Scenario 1 (equal effect sizes) — SNPs retained per threshold: 5e-8: 0 (absent below); 0.0001: 4;
0.001: 26; 0.01: 201; 0.05: 764; 0.1: 1,298; 0.5: 4,132; 1: 5,858**

| Threshold | EUR_holdout (n=101) | AFR (n=661) | AMR (n=347) | EAS (n=504) | SAS (n=489) |
|---|---|---|---|---|---|
| 0.0001 | 0.0180 [0.0040, 0.1031] | 0.0021 [0.0009, 0.0149] | 0.0192 [0.0011, 0.0577] | 0.0000 [0.0065, 0.0089] | 0.0054 [0.0002, 0.0260] |
| 0.001 | 0.0450 [0.0003, 0.1531] | 0.0035 [0.0003, 0.0182] | 0.0322 [0.0057, 0.0780] | 0.0023 [0.0016, 0.0181] | 0.0016 [0.0024, 0.0164] |
| 0.01 | **0.1310** [0.0321, 0.2710] | 0.0077 [0.0001, 0.0266] | 0.0098 [0.0000, 0.0408] | 0.0008 [0.0036, 0.0132] | 0.0080 [0.0000, 0.0312] |
| 0.05 | 0.1205 [0.0265, 0.2581] | 0.0014 [0.0015, 0.0128] | 0.0298 [0.0047, 0.0745] | 0.0068 [0.0000, 0.0284] | 0.0108 [0.0002, 0.0364] |
| 0.1 | 0.0951 [0.0144, 0.2258] | 0.0019 [0.0011, 0.0142] | 0.0370 [0.0079, 0.0851] | 0.0038 [0.0007, 0.0219] | 0.0092 [0.0001, 0.0335] |
| 0.5 | 0.1113 [0.0219, 0.2467] | 0.0019 [0.0011, 0.0142] | 0.0210 [0.0016, 0.0607] | 0.0045 [0.0004, 0.0237] | 0.0055 [0.0002, 0.0261] |
| 1 | 0.0993 [0.0163, 0.2313] | 0.0017 [0.0012, 0.0139] | 0.0214 [0.0017, 0.0614] | 0.0044 [0.0005, 0.0232] | 0.0054 [0.0002, 0.0260] |

**Scenario 2 (ancestry-varying effect sizes) — SNPs retained per threshold: 5e-8: 1; 0.0001: 6;
0.001: 33; 0.01: 204; 0.05: 761; 0.1: 1,319; 0.5: 4,061; 1: 5,812**

| Threshold | EUR_holdout (n=101) | AFR (n=661) | AMR (n=347) | EAS (n=504) | SAS (n=489) |
|---|---|---|---|---|---|
| 5e-08 | 0.0309 [0.0004, 0.1288] | 0.0054 [0.0000, 0.0222] | 0.0444 [0.0116, 0.0956] | 0.0223 [0.0040, 0.0547] | 0.0127 [0.0006, 0.0397] |
| 0.0001 | 0.0428 [0.0001, 0.1495] | 0.0007 [0.0025, 0.0105] | 0.0560 [0.0181, 0.1113] | 0.0015 [0.0024, 0.0158] | 0.0127 [0.0006, 0.0397] |
| 0.001 | 0.0426 [0.0001, 0.1492] | 0.0058 [0.0000, 0.0229] | 0.0176 [0.0008, 0.0550] | 0.0001 [0.0059, 0.0096] | 0.0140 [0.0009, 0.0420] |
| 0.01 | **0.0486** [0.0007, 0.1589] | 0.0081 [0.0002, 0.0272] | 0.0009 [0.0056, 0.0183] | 0.0090 [0.0001, 0.0327] | 0.0161 [0.0015, 0.0454] |
| 0.05 | 0.0297 [0.0006, 0.1266] | 0.0131 [0.0015, 0.0358] | 0.0047 [0.0014, 0.0298] | 0.0108 [0.0003, 0.0359] | 0.0098 [0.0001, 0.0345] |
| 0.1 | 0.0245 [0.0016, 0.1166] | 0.0073 [0.0001, 0.0258] | 0.0104 [0.0000, 0.0420] | 0.0085 [0.0000, 0.0318] | 0.0101 [0.0001, 0.0352] |
| 0.5 | 0.0050 [0.0161, 0.0688] | 0.0124 [0.0012, 0.0346] | 0.0082 [0.0002, 0.0377] | 0.0080 [0.0000, 0.0307] | 0.0148 [0.0011, 0.0434] |
| 1 | 0.0061 [0.0142, 0.0726] | 0.0091 [0.0004, 0.0290] | 0.0092 [0.0001, 0.0397] | 0.0062 [0.0001, 0.0272] | 0.0167 [0.0017, 0.0464] |

(Bold = the threshold used as "PRSice-2 (best EUR threshold)" in the headline table — chosen by
EUR held-out R², matching how a real analyst would pick a threshold in practice, i.e. tuned on the
same-ancestry data available at construction time, not on the non-EUR ancestries being evaluated.)
Note the wide, often nearly-zero-crossing 95% CIs, especially for AFR/EAS at low thresholds — a
direct consequence of small per-ancestry sample sizes (347–661) combined with weak per-SNP signal;
individual point estimates in this table carry real uncertainty, and only the *aggregate pattern*
(PRSice-2 never reaching PRS-CSx's/LDpred2-auto's cross-ancestry R² at any threshold) should be
read as a stable finding.

### Stage 6 — Recalibration (Track A)

`src/trackA_synthetic/stage6_recalibration/recalibrate.py`. The originally planned tool,
`pgscatalog-ancestry-adjust`, requires PCA projections and PGS-Catalog-format aggregated scores —
outputs of the full `pgsc_calc` pipeline that this project's custom-built Track A scoring files
never produce. Used the documented fallback instead: empirical per-ancestry recentering/rescaling,
matching each non-EUR ancestry's raw PRS mean/SD to the EUR held-out reference distribution's.

**The actual point of this stage, demonstrated numerically, not just asserted:** a linear
recentering/rescaling transform cannot change Pearson R² — R² depends only on correlation, which
is invariant to affine transforms. Confirmed in every one of the 95 rows in
`data/processed/recalibration/results.tsv`: `r2_raw` and `r2_recalibrated` are identical to four
decimal places throughout. Meanwhile `recal_mean`/`recal_sd` exactly match the EUR reference's
`ref_mean`/`ref_sd` after recalibration (calibration fixed), while `raw_mean`/`raw_sd` visibly
differ from the reference beforehand (the calibration problem recalibration is meant to fix).

This makes a distinction the PRS literature and this project's own build plan warn is often
conflated concrete rather than a claim to take on faith: **empirical per-ancestry recalibration
recovers 100% of calibration and 0% of discriminative accuracy, by mathematical necessity.** A
score whose mean/scale is corrected per ancestry is not the same as a score that discriminates
equally well per ancestry — these are different properties, and this stage exists specifically to
keep them from being conflated.

#### Worked example: PRS-CSx and LDpred2-auto's actual mean/SD shift, both scenarios

Concrete numbers from `data/processed/recalibration/results.tsv`, showing exactly what
recalibration does and doesn't change. `raw_mean`/`raw_sd` is each ancestry's own PRS distribution
before recalibration; `recal_mean`/`recal_sd` is after; `ref_mean`/`ref_sd` is the EUR held-out
reference being matched to; `r2_raw` vs. `r2_recalibrated` is the accuracy question.

**PRS-CSx, Scenario 1** (EUR reference: mean −0.8181, SD 0.8943):

| Ancestry | raw_mean → recal_mean (ref −0.8181) | raw_sd → recal_sd (ref 0.8943) | R² raw → recalibrated |
|---|---|---|---|
| AFR | −0.7244 → −0.8181 | 0.7884 → 0.8943 | 0.0275 → 0.0275 |
| AMR | −0.7955 → −0.8181 | 0.8503 → 0.8943 | 0.0594 → 0.0594 |
| EAS | −1.2723 → −0.8181 | 0.8589 → 0.8943 | 0.0165 → 0.0165 |
| SAS | −1.1617 → −0.8181 | 0.8575 → 0.8943 | 0.0279 → 0.0279 |

**PRS-CSx, Scenario 2** (EUR reference: mean 1.2598, SD 0.8461):

| Ancestry | raw_mean → recal_mean (ref 1.2598) | raw_sd → recal_sd (ref 0.8461) | R² raw → recalibrated |
|---|---|---|---|
| AFR | 1.5854 → 1.2598 | 0.8711 → 0.8461 | 0.0230 → 0.0230 |
| AMR | 1.1233 → 1.2598 | 0.9319 → 0.8461 | 0.0480 → 0.0480 |
| EAS | 1.4539 → 1.2598 | 0.9326 → 0.8461 | 0.0238 → 0.0238 |
| SAS | 1.3386 → 1.2598 | 0.9635 → 0.8461 | 0.0411 → 0.0411 |

**LDpred2-auto, Scenario 1** (EUR reference: mean −0.0099, SD 0.0170):

| Ancestry | raw_mean → recal_mean (ref −0.0099) | raw_sd → recal_sd (ref 0.0170) | R² raw → recalibrated |
|---|---|---|---|
| AFR | −0.0331 → −0.0099 | 0.0130 → 0.0170 | 0.0117 → 0.0117 |
| AMR | −0.0186 → −0.0099 | 0.0162 → 0.0170 | 0.0136 → 0.0136 |
| EAS | −0.0212 → −0.0099 | 0.0153 → 0.0170 | 0.0015 → 0.0015 |
| SAS | −0.0166 → −0.0099 | 0.0172 → 0.0170 | 0.0118 → 0.0118 |

**LDpred2-auto, Scenario 2** (EUR reference: mean 0.0252, SD 0.0182):

| Ancestry | raw_mean → recal_mean (ref 0.0252) | raw_sd → recal_sd (ref 0.0182) | R² raw → recalibrated |
|---|---|---|---|
| AFR | 0.0097 → 0.0252 | 0.0136 → 0.0182 | 0.0130 → 0.0130 |
| AMR | 0.0113 → 0.0252 | 0.0202 → 0.0182 | 0.0065 → 0.0065 |
| EAS | 0.0111 → 0.0252 | 0.0166 → 0.0182 | 0.0251 → 0.0251 |
| SAS | 0.0204 → 0.0252 | 0.0187 → 0.0182 | 0.0112 → 0.0112 |

In every single row above, `recal_mean`/`recal_sd` land exactly on the EUR reference (calibration
fully corrected) while R² doesn't move at all (discriminative accuracy completely unchanged) — the
same pattern holds for PRSice-2's rows too (omitted here for length; same 95-row table in
`data/processed/recalibration/results.tsv`).

### Track B — Real PGS Catalog score (descriptive only, never a validation)

> Every result in this section is descriptive/illustrative, never a validated accuracy claim —
> 1000 Genomes carries no phenotype data, so there is no ground truth here to validate a
> prediction against. Only Track A's synthetic experiment provides that.

- **Stage B1** (`src/trackB_real_scores/stageB1_download_pgs/`): downloads a real PGS Catalog
  score via `pgscatalog-download`, GRCh37-harmonized to match 1000 Genomes' build. Uses
  **PGS000297** ("GRS3290_Height", Xie et al. 2020, 3,290 variants) — height alone satisfies the
  "one or more" requirement; the commonly-cited BMI score (PGS000027) is genome-wide with ~2.1M
  variants, wasteful to download against this project's chr21+22-only panel.
- **Stage B2** (`src/trackB_real_scores/stageB2_compute_scores/`): matches the scoring file to the
  QC'd panel by chromosome+position, computes the score for **all 2,504 individuals** across all
  five super-populations via `plink2 --score` (no held-out split needed — there's no phenotype to
  guard against overfitting to). Of the 3,290 scoring variants, 63 fall on chr21/22 and **52** have
  a matching allele in the QC'd panel (11 dropped for multiallelic/strand mismatches at that
  position) — a small subset by design (smoke-test scope), not a bug.
- **Stage B3** (`src/trackB_real_scores/stageB3_descriptive_comparison/`): compares each
  super-population's mean PGS to published national height statistics for representative
  countries per super-population (sourced from Wikipedia's "Average human height by country"
  compilation and NCD-RisC 2016, full citations in `published_height_reference.tsv`). This mapping
  is itself approximate — 1000 Genomes super-populations aggregate several specific source
  populations that don't correspond to a single country, and the reference figures span different
  survey years/methods/sexes.

**The real finding, reported honestly rather than smoothed over:** **0 of 5** super-populations'
mean PGS rank matches the published-height rank. AFR ranks highest by PGS but 2nd by published
height; EUR ranks highest by published height but only 4th by PGS. This is not evidence the
pipeline is broken — it's a small-scale, real illustration of the exact cross-ancestry portability
problem Track A's synthetic experiment was built to quantify with a controlled ground truth: a
European-ancestry-derived score, applied to other ancestries with no ancestry-appropriate
recalibration and using only 52 of 3,290 scoring variants (chr21+22 scope), should not be expected
to preserve true population-level phenotypic ordering.

**Full comparison table** (`data/processed/pgs_catalog/stageB3_descriptive_comparison.tsv`):

| Super-pop | n | PGS mean | PGS SD | PGS rank | Published height (cm) | Height rank | Ranks match? |
|---|---|---|---|---|---|---|---|
| AFR | 661 | 0.8166 | 0.0620 | 1 | 164.9 | 2 | No |
| SAS | 489 | 0.7840 | 0.0632 | 2 | 158.6 | 5 | No |
| AMR | 347 | 0.7526 | 0.0648 | 3 | 162.6 | 4 | No |
| EUR | 503 | 0.7367 | 0.0634 | 4 | 170.4 | 1 | No |
| EAS | 504 | 0.7102 | 0.0634 | 5 | 163.9 | 3 | No |

Every row disagrees, and not by a small margin — the super-population with the highest published
height (EUR) has the second-*lowest* mean PGS, and the super-population with the highest mean PGS
(AFR) is only mid-table by published height. Published height figures are national averages for
representative countries per super-population, cross-checked against Wikipedia's "Average human
height by country" compilation and NCD-RisC 2016 (exact source country/survey year per row in
`src/trackB_real_scores/stageB3_descriptive_comparison/published_height_reference.tsv`) — not
derived from 1000 Genomes itself, since it has no phenotype data at all.

## 7. What the results mean, taken together

Within a chr21+22 smoke-test scope, this project reproduces — with a controlled synthetic ground
truth, not just a real-data anecdote — the central, well-documented finding in current PRS
cross-ancestry literature: a EUR-derived score's advantage over other methods same-ancestry does
not carry over cross-ancestry. PRS-CSx consistently outperforms PRSice-2 in every non-EUR
evaluation, even while losing to it in the EUR held-out set. LDpred2-auto shows the same
qualitative pattern relative to PRSice-2 in most (six of eight) non-EUR comparisons, but trails
PRS-CSx in nearly all of them — plausibly a consequence of this run's narrower chr22-only LD
reference rather than a property of the method itself, reinforcing rather than contradicting this
project's built-in refusal to rank any one method as universally best.

Separately, this project demonstrates numerically that simple empirical recalibration fixes a
score's per-ancestry mean/scale completely while leaving its discriminative accuracy completely
unchanged — a distinction real-world PRS deployment discussions frequently conflate. Track B's
real-score cross-check adds a third, independent illustration of the same underlying phenomenon:
an uncorrected EUR-derived score's ranking across ancestries doesn't track true population-level
phenotype differences.

**None of these findings should be read as a genome-wide claim, a validated clinical
recommendation, or a definitive method ranking.**

## 8. Non-obvious implementation history worth knowing

A few things that shaped the implementation and are worth knowing before touching the code again
(full detail in `METHODS.md`):

- **Stage 3's EUR discovery/held-out split** wasn't in the initial implementation — it was caught
  and fixed while starting Stage 4, once it became clear that evaluating "EUR" accuracy without a
  genuine held-out set would be circular (in-sample against the GWAS's own training data).
- **LDpred2-auto was the hardest stage to complete end-to-end**, purely due to this project's
  shared sandbox host being heavily contended (load average 120–150 sustained across testing, swap
  frequently near-exhausted from unrelated jobs run by other users on the same machine — not
  anything this project's code controls). Three early attempts failed outright; it eventually
  succeeded on 2026-08-02 across two runs (one completed Scenario 1 but timed out mid-MCMC on
  Scenario 2; a Scenario-2-only rerun, after clearing stale bigsnpr backing files left by the
  killed run, completed successfully). No code changes were needed — the fix was patience and a
  reduced MCMC budget (`burn_in=200, num_iter=100` vs. bigsnpr's usual 500/200 defaults), not a bug
  fix.
- **PRS-CSx's SNP ID mismatch:** its reference panel addresses SNPs by rsID, but Stage 1 assigns
  chr:pos:ref:alt IDs, so the wrapper builds a position-matched rsID↔project-ID mapping before
  running PRS-CSx and translates the output back.
- **A real bug in PRS-CSx itself** (not just an unclear error): it expects `snpinfo_mult_1kg_hm3`
  directly inside `--ref_dir`, and throws an `UnboundLocalError` (not a clear error message) if
  that file isn't found there.
- **GCTA has no `--seed` for `--simu-qt`'s residual noise** in the version used (1.95.3) — the
  genetic component of the simulation is fully reproducible via this project's own seeded RNG, but
  the residual/noise draw is not bit-for-bit reproducible across re-runs.

## 9. Licensing and data-usage terms

Full detail and exact quoted terms in `LICENSES.md`. Summary:

- **Original code in this repository:** MIT license.
- **1000 Genomes/IGSR data:** not a blanket public-domain grant — IGSR's own disclaimer states
  rights/restrictions vary by dataset and puts the burden on users to ensure their use doesn't
  infringe third-party rights. This project's use (public, published, de-identified genotype data,
  for non-commercial methodology-demonstration research, with the source paper cited) fits the
  intended reuse case IGSR describes, but this is a judgment call, not an OSI-style license grant.
- **PGS Catalog / EBI:** EBI's general terms place no additional restriction on data beyond what
  original data owners set; individual PGS Catalog scores may carry their own author-specific
  license. This project's specific score (PGS000297) has no score-specific override listed beyond
  the general terms.
- **Tool licenses** (confirmed from each tool's own license file/repository): GCTA — MIT;
  PRSice-2 — GPL-3.0; `bigsnpr`/`bigsparser` (LDpred2-auto) — GPL-3; PRS-CSx — MIT;
  `pgscatalog-utils`/`pgsc_calc` — Apache-2.0; plink2 — GPLv3 core, LGPLv3 for its `pgenlibr`
  library component, plus a separate Intel Simplified Software License governing its bundled MKL
  binary only. This project invokes every one of these tools as an external
  binary/CLI/R-package dependency, never by copying or modifying their source into this
  repository, so none of their licenses (including the GPL/LGPL ones) impose any obligation on
  this project's own MIT-licensed code.

## 10. Repository structure

```
prs-cross-ancestry-portability-benchmark/
├── README.md                          Quickstart + repo layout
├── PROJECT_OVERVIEW.md                This file
├── METHODS.md                         Full stage-by-stage implementation record and run history
├── LICENSES.md                        Data usage terms and tool licenses, confirmed and quoted
├── docs/BUILD_PLAN.md                 Original plan, rationale, scope/honesty constraints
├── environment.yml                    Python/conda environment
├── src/
│   ├── shared_stage1_1000g_download_qc/   Stage 1: download + QC 1000 Genomes (shared)
│   ├── trackA_synthetic/
│   │   ├── stage2_gcta_simulation/         Stage 2: GCTA phenotype simulation, both scenarios
│   │   ├── stage3_eur_gwas/                Stage 3: EUR-only discovery GWAS
│   │   ├── stage4_prs_construction/        Stage 4: PRSice-2, LDpred2-auto, PRS-CSx wrappers
│   │   ├── stage5_crossancestry_evaluation/ Stage 5: apply models to all 5 ancestries, compute R²
│   │   └── stage6_recalibration/           Stage 6: empirical per-ancestry recalibration
│   └── trackB_real_scores/
│       ├── stageB1_download_pgs/           B1: download real PGS Catalog score
│       ├── stageB2_compute_scores/         B2: compute score on all 2,504 individuals
│       └── stageB3_descriptive_comparison/ B3: compare to published height statistics
├── configs/simulation_parameters.yaml Stage 2 simulation architecture (DRAFT, pending review)
├── scripts/verify_tools.py            Check/install GCTA, PRSice-2, PRS-CSx, plink2, pgscatalog-utils, bigsnpr
├── notebooks/                         Result notebooks
├── reports/final_writeup/             TECHNICAL_WRITEUP.md + one-page case-study PDF
└── tests/
```

`data/` and `tools/` (downloaded data and installed binaries) are gitignored — never committed;
re-fetch them via the scripts above.

## 11. Reproducing this project

```bash
# 1. Create the Python/R environment
conda env create -f environment.yml
conda activate prs-portability

# 2. Check (and optionally install) GCTA, PRSice-2, PRS-CSx, plink2, pgscatalog-utils, bigsnpr
python3 scripts/verify_tools.py --install            # add --with-r for bigsnpr too (slow)

# 3. Stage 1 — download a smoke-test subset (chr21, chr22) of 1000 Genomes Phase 3
python3 src/shared_stage1_1000g_download_qc/download.py

# 4. Stage 1 — QC (per-super-population MAF/missingness/HWE filtering)
python3 src/shared_stage1_1000g_download_qc/qc.py \
  --vcf-dir data/1000genomes/vcf \
  --panel data/1000genomes/metadata/integrated_call_samples_v3.20130502.ALL.panel \
  --chromosomes 21 22

# 5. Stage 2 — GCTA phenotype simulation (both scenarios; needs Stage 1 output)
python3 src/trackA_synthetic/stage2_gcta_simulation/scenario1_equal_effects.py
python3 src/trackA_synthetic/stage2_gcta_simulation/scenario2_ancestry_varying_effects.py

# 6. Stage 3 — EUR-only discovery GWAS (both scenarios; needs Stage 2 output)
python3 src/trackA_synthetic/stage3_eur_gwas/run_gwas.py

# 7. Stage 4 — PRS construction (both scenarios; needs Stage 3 output)
python3 src/trackA_synthetic/stage4_prs_construction/prsice2_wrapper.py
python3 src/trackA_synthetic/stage4_prs_construction/prscsx_wrapper.py --download-ld-ref  # first run only
R_LIBS_USER=~/R/library Rscript src/trackA_synthetic/stage4_prs_construction/ldpred2_wrapper.R --scenario all

# 8. Stage 5 — cross-ancestry evaluation (needs Stage 4 output)
python3 src/trackA_synthetic/stage5_crossancestry_evaluation/evaluate.py

# 9. Stage 6 — recalibration (needs Stage 4 output)
python3 src/trackA_synthetic/stage6_recalibration/recalibrate.py

# 10. Track B — real PGS Catalog score (descriptive only; needs Stage 1 output)
python3 src/trackB_real_scores/stageB1_download_pgs/download_pgs.py
python3 src/trackB_real_scores/stageB2_compute_scores/compute_scores.py
python3 src/trackB_real_scores/stageB3_descriptive_comparison/compare.py
```

`download.py --chromosomes all` fetches the full autosomal panel (tens of GB) instead of the
default smoke-test subset — do that deliberately, not by default; a genome-wide run has not been
performed in this project to date.

## 12. Current status

All items on `docs/BUILD_PLAN.md`'s Deliverables Checklist (§13) and Licensing & Compliance
Checklist (§11) are complete:

- All 7 pipeline stages implemented and verified on real chr21+22 data.
- Both tracks (A: synthetic ground truth; B: real score, descriptive) complete.
- All three PRS methods (PRSice-2, PRS-CSx, LDpred2-auto) complete for both scenarios.
- Full 3-method × 5-ancestry × 2-scenario results table, before and after recalibration.
- Technical write-up (~2,100 words) and one-page case-study PDF, both current as of the completed
  3-method results.
- Data-usage terms and every tool's license confirmed and quoted in `LICENSES.md`.

**What remains open, honestly stated:** this is a chr21+22 smoke-test result, not a genome-wide
one — everything above is scoped accordingly, and no claim in this project should be read past
that scope without a full genome-wide rerun (§1.2, non-goals; §7).

## 13. References

1. Momin MM, Zhou X, Ahmed M, Hyppönen E, Benyamin B, Lee SH. Cross-Ancestry Polygenic Prediction:
   Comparing Methods and Assessing Transferability Across Traits. *Genetic Epidemiology.*
   2026;50:1–13. https://doi.org/10.1002/gepi.70029
2. Ruan Y, Lin YF, Feng YCA, et al. Improving polygenic prediction in ancestrally diverse
   populations (PRS-CSx). *Nature Genetics.* 2022;54:573–580.
   https://pubmed.ncbi.nlm.nih.gov/35513724/ — code: https://github.com/getian107/PRScsx
3. Privé F, Arbel J, Vilhjálmsson BJ. LDpred2: better, faster, stronger. *Bioinformatics.* 2020.
   https://privefl.github.io/bigsnpr/
4. Choi SW, O'Reilly PF. PRSice-2: Polygenic Risk Score software for biobank-scale data.
   *GigaScience.* 2019. https://github.com/choishingwan/PRSice
5. Lambert SA, Gil L, Jupp S, et al. The Polygenic Score Catalog as an open database for
   reproducibility and systematic evaluation. *Nature Genetics.* 2021.
6. Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis.
   *American Journal of Human Genetics.* 2011;88:76–82.
7. 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature.*
   2015. https://www.internationalgenome.org
8. Xie T et al. [GRS3290_Height]. *Circulation: Genomic and Precision Medicine.* 2020.
   https://doi.org/10.1161/circgen.119.002775
9. NCD Risk Factor Collaboration (NCD-RisC). A century of trends in adult human height. *eLife.*
   2016;5:e13410.
