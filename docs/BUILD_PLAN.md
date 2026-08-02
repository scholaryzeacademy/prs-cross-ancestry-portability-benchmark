# BUILD_PLAN.md
## Polygenic Risk Score Cross-Ancestry Portability & Recalibration Benchmark
### Real Multi-Ancestry Genotypes + Controlled Synthetic Ground Truth + Real Published Scores

---

> **Project type:** Second-wave open-source portfolio piece (Genomics & Variant Interpretation vertical)
> **Status:** Planning
> **License target:** MIT/Apache-2.0 for original code; third-party data/tool licenses reviewed individually (see §11)

---

## 0. One-Paragraph Summary

This project has the same two-track structure as the Clinical/RWE build plan, for the same underlying reason: **1000 Genomes Phase 3 — the standard, fully public, multi-ancestry genotype reference panel — has real genotypes but no real disease/trait phenotypes**, so a genuine predictive-accuracy validation needs a synthetic phenotype with a known, controlled genetic architecture layered on top of real genotypes. **Track A** uses GCTA to simulate quantitative phenotypes with documented causal variants and heritability directly on real 1000 Genomes genotypes across five super-populations (AFR, AMR, EAS, EUR, SAS), letting us precisely measure how much predictive accuracy each PRS method loses when a score trained on one ancestry (mimicking the field's real, well-documented European-ancestry GWAS bias) is applied to others — and how much a recalibration step recovers. **Track B** applies real, published PGS Catalog scores to the same real genotypes as a descriptive, illustrative cross-check, explicitly flagged as non-validating given the lack of real phenotypes in 1000 Genomes. This project directly operationalizes the single most active, most publicly discussed generalization-gap problem in genomics today: PRS built predominantly on European-ancestry GWAS systematically underperforming everywhere else [Momin et al., *Genetic Epidemiology*, 2026].

---

## 1. Goals & Non-Goals

### 1.1 Goals
- Quantify, with a known ground truth, how much predictive accuracy (R²) is lost when a PRS trained on one ancestry group is applied to others, using real multi-ancestry LD structure and allele frequencies from 1000 Genomes.
- Compare at least three PRS construction methods — a simple clumping-and-thresholding baseline (PRSice-2), and two more sophisticated Bayesian/multi-ancestry methods (LDpred2, PRS-CSx) — honestly, following the same "include a naive baseline" discipline used in every other project in this portfolio.
- Apply and evaluate a recalibration step (ancestry-adjustment) and report, honestly, how much of the portability gap it actually closes versus how much it doesn't.
- Cross-check the synthetic-ground-truth findings against real, published PGS Catalog scores applied to the same genotype panel, clearly distinguishing descriptive illustration from validated predictive accuracy.

### 1.2 Non-Goals
- We are **not** claiming to produce clinically deployable PRS recalibration recommendations — this is a methodology demonstration on public reference data, not a validated clinical tool.
- We are **not** attempting to access UK Biobank, All of Us, or any other biobank requiring a data-access application — this project is deliberately scoped to fully public, no-credentialing data (1000 Genomes, PGS Catalog), which is both a genuine constraint and, honestly stated as such, itself part of the project's transparent scope.
- We are **not** trying to definitively rank PRS methods in general — current literature is explicit that no single method wins universally across traits [Momin et al., 2026], and our own single-scenario simulation should not be overgeneralized into a universal claim either.

---

## 2. Definition of Done

1. Real 1000 Genomes Phase 3 genotypes across all five super-populations are ingested and quality-controlled.
2. At least two synthetic phenotype scenarios are generated via GCTA with documented causal variants and heritability: one with **identical effect sizes across ancestries** (isolating the pure LD/allele-frequency portability problem) and one with **deliberately different effect sizes across ancestries** (modeling genuine effect-size heterogeneity), mirroring the two-scenario sensitivity-analysis discipline used in the Clinical/RWE build plan.
3. Three PRS methods (PRSice-2 clumping+thresholding, LDpred2, PRS-CSx) are each trained on a EUR-only "discovery GWAS" (deliberately mimicking real-world European-ancestry GWAS bias) and evaluated for predictive accuracy in held-out EUR, AFR, AMR, EAS, and SAS samples.
4. A recalibration/ancestry-adjustment step is applied and its effect on predictive accuracy and calibration is reported honestly, including if it doesn't fully close the gap.
5. At least one real, published PGS Catalog score (e.g., for height or BMI) is applied to the same 1000 Genomes genotypes via `pgsc_calc`/`pgscatalog-utils`, with results presented explicitly as descriptive/illustrative rather than a validated accuracy claim.
6. Public GitHub repo, technical write-up, and case-study PDF, per the standard deliverables checklist.

---

## 3. Architecture Overview

```
      TRACK A: SYNTHETIC GROUND TRUTH                TRACK B: REAL PUBLISHED SCORES (descriptive)
 ┌─────────────────────────────┐              ┌──────────────────────────────────┐
 │ STAGE 1: Download 1000        │              │ STAGE 1 (shared): Same 1000       │
 │ Genomes Phase 3 (5 super-      │              │ Genomes genotype panel             │
 │ populations, real genotypes)   │              └────────────────┬───────────────────┘
 └──────────────┬───────────────┘                               │
                │                              ┌────────────────▼───────────────────┐
 ┌──────────────▼───────────────┐              │ STAGE B1: Download real PGS         │
 │ STAGE 2: GCTA phenotype        │              │ Catalog scores (e.g., height, BMI) │
 │ simulation — Scenario 1        │              │ via pgscatalog-utils                │
 │ (equal effect sizes) and       │              └────────────────┬───────────────────┘
 │ Scenario 2 (ancestry-varying   │                               │
 │ effect sizes)                  │              ┌────────────────▼───────────────────┐
 └──────────────┬───────────────┘              │ STAGE B2: Compute scores on 1000     │
                │                              │ Genomes samples via pgsc_calc         │
 ┌──────────────▼───────────────┐              └────────────────┬───────────────────┘
 │ STAGE 3: EUR-only "discovery   │                               │
 │ GWAS" on simulated phenotype   │              ┌────────────────▼───────────────────┐
 └──────────────┬───────────────┘              │ STAGE B3: Descriptive comparison      │
                │                              │ against published population-level    │
 ┌──────────────▼───────────────┐              │ trait statistics (illustrative only)  │
 │ STAGE 4: Construct PRS —       │              └────────────────────────────────────┘
 │ PRSice-2, LDpred2, PRS-CSx     │
 └──────────────┬───────────────┘
                │
 ┌──────────────▼───────────────┐
 │ STAGE 5: Evaluate predictive   │
 │ accuracy (R²) across all 5     │
 │ super-populations               │
 └──────────────┬───────────────┘
                │
 ┌──────────────▼───────────────┐
 │ STAGE 6: Apply recalibration/  │
 │ ancestry-adjustment            │
 │ (pgscatalog-ancestry-adjust)   │
 └──────────────┬───────────────┘
                │
                └──────────────┬──────────────────────────────────┘
                               │
                 ┌─────────────▼─────────────┐
                 │ STAGE 7: Combined Honest    │
                 │ Report: synthetic ground-   │
                 │ truth results + descriptive │
                 │ real-score cross-check      │
                 └────────────────────────────┘
```

---

## 4. Data Sources & Tools (Verified)

| Resource | Role | Access | Notes |
|---|---|---|---|
| **1000 Genomes Phase 3 (or 30x high-coverage resequencing)** | Real multi-ancestry genotype panel | `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502` (Phase 3, GRCh37); also via IGSR Data Portal, Globus, or Aspera | 2,504 individuals (or 3,202 in the high-coverage 30x resequenced release) across 26 populations grouped into 5 continental super-populations: AFR, AMR, EAS, EUR, SAS. **Confirmed: no disease/complex-trait phenotype data is distributed with this panel** — only genotypes, sex, and family relationships — which is exactly why Track A requires simulated phenotypes. |
| **PGS Catalog** | Real published PRS scoring files (Track B) | `pgscatalog.org` — REST API, FTP (`ftp://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS######/ScoringFiles/`), and the official `pgscatalog-utils`/`pgsc_calc` tooling | Confirmed pip-installable: `pip install pgscatalog-utils`, providing `pgscatalog-download`, `pgscatalog-match`, `pgscatalog-aggregate`, and — directly relevant to this project — **`pgscatalog-ancestry-adjust`**, an official tool specifically for "adjust[ing] calculated PGS in the context of genetic ancestry." `pgsc_calc` is the accompanying Nextflow pipeline for reproducible score calculation on user genotype data. |
| **GCTA** | Phenotype simulation on real genotypes | `cnsgenomics.com/software/gcta/` — free, open source | `--simu-qt` (quantitative trait) and `--simu-cc` (case-control) flags simulate phenotypes from real genotype data given a specified causal-variant list and heritability — a well-established, widely used approach in the population genetics literature specifically for this kind of method-validation exercise. |
| **PRS-CSx** | PRS construction method 1 (multi-ancestry Bayesian) | `github.com/getian107/PRScsx` | Couples genetic effects across ancestries via a shared continuous shrinkage prior; identified as a top performer for highly polygenic traits in the most current cross-ancestry method comparison [Momin et al., 2026]. |
| **LDpred2** | PRS construction method 2 (Bayesian, single- or multi-ancestry via `bigsnpr`) | R package `bigsnpr` (CRAN) | Standard, well-maintained Bayesian PRS method; identified alongside GBLUP as a top performer for highly polygenic traits like height/BMI [Momin et al., 2026]. |
| **PRSice-2** | PRS construction method 3 (simple clumping + thresholding baseline) | `github.com/choishingwan/PRSice` | The deliberately simple baseline method in this project — following the "always include a naive baseline" discipline used across the whole portfolio; notably, current literature finds PRSice actually performs *best* for less polygenic traits like cholesterol, which is itself a useful, non-obvious finding to reproduce and report [Momin et al., 2026]. |

**Action item before build starts:** confirm the exact current license terms for any specific PGS Catalog scores selected for Track B — the Catalog's overall access is governed by EBI's terms of use, but individual scores may carry their own specific license (CC or otherwise), and this must be checked per-score, not assumed uniform across the Catalog.

---

## 5. Tech Stack

| Layer | Choice | Rationale |
|---|---|---|
| Language | Python 3.10+ and R 4.x (for `bigsnpr`/LDpred2) | PRS methodology tooling spans both ecosystems; budget for both environments |
| Environment management | Conda + `environment.yml` for Python; a documented `renv`/CRAN lockfile for the R side | Keeps the dual-language stack reproducible |
| Genotype QC/manipulation | `PLINK`/`PLINK2` | The universal standard for genotype data manipulation, filtering, and format conversion across every tool in this stack |
| Phenotype simulation | `GCTA` | As above — confirmed standard tool for this exact use case |
| PRS Catalog tooling | `pgscatalog-utils` (`pip install pgscatalog-utils`), `pgsc_calc` (Nextflow) | Official, confirmed-working tooling; use rather than hand-rolling scoring-file parsing |
| Cross-ancestry PRS | `PRS-CSx` (Python, MCMC-based) | GPU not required — this is a CPU-bound MCMC method; budget realistic wall-clock time (can be slow on large variant sets) |
| Bayesian PRS | `LDpred2` via `bigsnpr` (R) | — |
| Simple baseline PRS | `PRSice-2` (C++/R) | — |
| Statistical evaluation | `scikit-learn`/`statsmodels` (Python) or base R for R² computation, confidence intervals | Standard |
| Reporting | Jinja2-templated HTML + notebook, consistent with the rest of the portfolio | — |
| CI | GitHub Actions smoke test on a small chromosome-arm subset with a handful of simulated causal variants | Full-genome PRS-CSx/LDpred2 runs are too heavy for CI — run those manually/on-demand |

---

## 6. Pipeline Stage Detail

### Track A

**Stage 1 — Download and QC 1000 Genomes Phase 3**
- Download VCFs via FTP/Globus, convert to PLINK binary format, apply standard QC filters (minor allele frequency, missingness, Hardy-Weinberg equilibrium within each super-population)
- **Owner:** Bioinformatics

**Stage 2 — Simulate Phenotypes (Two Scenarios)**
- Select a documented set of causal SNPs (e.g., a few hundred, spread across multiple chromosomes) and a target heritability (e.g., h² = 0.5, a realistic value for a moderately polygenic trait)
- **Scenario 1 (equal effect sizes across ancestries):** run GCTA `--simu-qt` with the same causal-variant effect sizes applied uniformly across all individuals regardless of super-population — this isolates the portability problem caused purely by real differences in allele frequency and LD structure across ancestries, holding biology constant
- **Scenario 2 (ancestry-varying effect sizes):** perturb the causal-variant effect sizes differently per super-population before simulation, modeling genuine gene-by-ancestry effect heterogeneity — document the exact perturbation scheme
- **Owner:** Biostatistics (owns the ground-truth architecture design, exactly mirroring the RWE project's Stage A)

**Stage 3 — EUR-Only Discovery GWAS**
- Run a GWAS (simple linear regression per-SNP, or `PLINK --linear`) on the simulated phenotype using **only the EUR super-population samples** — deliberately reproducing the real-world bias that most published GWAS are European-ancestry-derived
- **Owner:** Bioinformatics

**Stage 4 — Construct PRS via Three Methods**
- Feed the EUR-only GWAS summary statistics into PRSice-2 (clumping+thresholding), LDpred2, and PRS-CSx, producing three independently-constructed PRS models
- **Owner:** AI Engineer + Bioinformatics

**Stage 5 — Evaluate Cross-Ancestry Predictive Accuracy**
- Apply each of the three PRS models to held-out samples from **all five** super-populations (including a held-out EUR subset as the "easy," same-ancestry baseline case)
- Compute R² (or an appropriate pseudo-R²/AUC for a case-control scenario if that variant is also run) per method per ancestry group, with confidence intervals
- Repeat for both Scenario 1 and Scenario 2, and explicitly compare — does the portability gap look different when it's purely an LD/frequency effect (Scenario 1) versus when real effect-size heterogeneity is also present (Scenario 2)?
- **Owner:** Biostatistics (final sign-off)

**Stage 6 — Apply Recalibration**
- Apply `pgscatalog-ancestry-adjust` (or implement a standard empirical ancestry-matched recentering/rescaling approach if the tool's assumptions don't cleanly fit our synthetic setup) to the cross-ancestry PRS predictions
- Report, honestly, how much of the Stage 5 portability gap this recovers — and be explicit if it recovers calibration (the score's mean/scale being appropriate per ancestry) without fully recovering discriminative accuracy (R²), since these are different things and conflating them would be a real analytical error worth avoiding
- **Owner:** Biostatistics

### Track B

**Stage B1 — Download Real PGS Catalog Scores**
- Use `pgscatalog-download` to retrieve one or more real, well-powered PGS Catalog scores for a common, well-studied trait (height and/or BMI, consistent with the specific traits emphasized in current cross-ancestry PRS literature)
- **Owner:** Bioinformatics

**Stage B2 — Compute Scores on 1000 Genomes Samples**
- Run `pgsc_calc` (or the underlying `pgscatalog-match`/`pgscatalog-aggregate` steps directly) to compute the real published score for every 1000 Genomes individual
- **Owner:** Bioinformatics

**Stage B3 — Descriptive Cross-Check**
- Compare the distribution of computed scores across super-populations to published, external, population-level trait statistics from the literature (not from 1000 Genomes itself, since it lacks phenotypes) — explicitly framed as an illustrative plausibility check, not a validated accuracy measurement
- **State this limitation clearly and prominently** in the write-up — this track exists to connect Track A's synthetic findings to a real, recognizable, named PGS Catalog score, not to independently validate predictive accuracy
- **Owner:** Genetics domain expert + Biostatistics

---

## 7. Repository Structure

```
prs-cross-ancestry-portability-benchmark/
├── README.md
├── METHODS.md
├── LICENSES.md
├── environment.yml
├── renv.lock                        # R-side reproducibility
├── .github/workflows/ci.yml
├── src/
│   ├── shared_stage1_1000g_download_qc/
│   ├── trackA_synthetic/
│   │   ├── stage2_gcta_simulation/
│   │   │   ├── scenario1_equal_effects.py
│   │   │   └── scenario2_ancestry_varying_effects.py
│   │   ├── stage3_eur_gwas/
│   │   ├── stage4_prs_construction/
│   │   │   ├── prsice2_wrapper.py
│   │   │   ├── ldpred2_wrapper.R
│   │   │   └── prscsx_wrapper.py
│   │   ├── stage5_crossancestry_evaluation/
│   │   └── stage6_recalibration/
│   └── trackB_real_scores/
│       ├── stageB1_download_pgs/
│       ├── stageB2_compute_scores/
│       └── stageB3_descriptive_comparison/
├── notebooks/
│   ├── 01_1000g_qc_summary.ipynb
│   ├── 02_scenario1_results.ipynb
│   ├── 03_scenario2_results.ipynb
│   ├── 04_recalibration_results.ipynb
│   └── 05_trackB_descriptive_crosscheck.ipynb
├── configs/
│   └── simulation_parameters.yaml    # causal variants, effect sizes, heritability — fully documented
├── tests/
└── reports/
    └── final_writeup/
```

---

## 8. Milestones & Timeline

| Phase | Duration | Deliverable |
|---|---|---|
| **Phase 0 — Setup & Diligence** | Week 1 | 1000 Genomes download and QC pipeline working; GCTA, PRSice-2, LDpred2/bigsnpr, PRS-CSx, and pgscatalog-utils all installed and verified; simulation architecture (Stage 2) designed and reviewed by Biostatistics |
| **Phase 1 — Synthetic Phenotype Generation** | Week 2 | Stage 2 complete for both scenarios |
| **Phase 2 — Discovery GWAS & PRS Construction** | Weeks 3–4 | Stages 3–4 complete across all three methods |
| **Phase 3 — Cross-Ancestry Evaluation** | Weeks 4–5 | Stage 5 complete; this is the phase producing the project's headline result — protect time here |
| **Phase 4 — Recalibration** | Week 6 | Stage 6 complete |
| **Phase 5 — Track B (Real Scores)** | Week 6–7 (can run partly in parallel with Phase 3–4) | Stages B1–B3 complete |
| **Phase 6 — Reporting & Publication** | Weeks 7–8 | Combined write-up; repo finalized; publication |

**Total: ~8 weeks**, comparable to the first-wave genomics project — this one has real statistical/simulation-design weight in Stage 2 and real compute weight in Stage 4 (PRS-CSx's MCMC sampling can be slow), so it isn't a fast build despite using well-established tools.

---

## 9. Validation & QA Strategy

1. **Sanity-check the GCTA simulation itself first** — confirm that a simple GWAS + PRS run within the EUR-only held-out set (same-ancestry train/test) recovers close to the expected R² given the simulated heritability, before trusting any cross-ancestry comparison; if the same-ancestry baseline doesn't behave sensibly, the whole downstream comparison is unreliable.
2. **Compare Scenario 1 vs. Scenario 2 explicitly** — if the portability gap looks similar in both, that's an interesting finding suggesting LD/frequency differences alone (not effect-size heterogeneity) drive most of the portability loss in this simulation; if they look very different, that's equally worth reporting and discussing.
3. **Report per-method, per-ancestry, per-scenario results as a full table, not a single summary statistic** — collapsing this into "method X is best" would repeat exactly the overclaiming mistake flagged in every other project in this portfolio; the honest finding from current literature is that no method wins universally, and our own simulation should be reported with the same care.
4. **Be explicit about Track B's descriptive-only status** in every place it's referenced in the write-up — this is a repeated-emphasis requirement, not a single disclaimer, following the same discipline used for the Clinical/RWE project's synthetic-data framing.
5. **Document every GCTA simulation parameter** (causal variant list, effect sizes, heritability, per-ancestry perturbation scheme for Scenario 2) in `configs/simulation_parameters.yaml` so the entire exercise is independently reproducible.

---

## 10. Compute & Infrastructure Requirements

- **No GPU required** — PRS-CSx, LDpred2, and PRSice-2 are all CPU-bound (PRS-CSx uses MCMC sampling, which can be slow but doesn't need a GPU).
- **1000 Genomes Phase 3 data is a few GB in VCF form** — modest by genomics standards, manageable on standard cloud storage or even a capable local machine.
- **PRS-CSx's MCMC sampling is the most likely compute bottleneck** — budget realistic wall-clock time during Phase 0 by running it on a small chromosome subset first, before committing to genome-wide timing estimates for the full pipeline.
- **DevOps/Cloud team member's core contribution:** managing the dual Python/R environment reproducibly (Docker or a well-documented Conda + `renv` combination), and provisioning adequate CPU (not GPU) compute for the PRS-CSx MCMC stage.

---

## 11. Licensing & Compliance Checklist

- [x] Confirm 1000 Genomes/IGSR data usage terms — not a blanket public-domain grant; see `LICENSES.md` for the exact disclaimer language and this project's specific use-case judgment
- [x] Confirm PGS Catalog's overall EBI terms of use, and check the specific license of any individual scoring file used in Track B, since some carry their own CC or non-commercial terms distinct from the Catalog's general terms — see `LICENSES.md`
- [x] Confirm GCTA, PRSice-2, LDpred2/`bigsnpr`, PRS-CSx, and `pgscatalog-utils`/`pgsc_calc` license terms (check each repository's LICENSE file directly) — see `LICENSES.md`'s Tool Licenses table
- [x] Cite Momin et al. 2026, Ruan et al. 2022 (PRS-CSx), and the PGS Catalog's own publications clearly as the methodological foundation this project builds on — see `reports/final_writeup/TECHNICAL_WRITEUP.md` References
- [x] Choose and apply an explicit open-source license for original code — MIT, see `LICENSE`

---

## 12. Team & Role Allocation

| Role | Primary Responsibility |
|---|---|
| **Bioinformatics** | Stage 1 (1000 Genomes QC), Stage 3 (GWAS), Track B's data retrieval/computation |
| **Biostatistics** | Stage 2 (simulation architecture — the project's technical core), Stage 5 (evaluation), Stage 6 (recalibration) — final sign-off on every quantitative claim |
| **AI Engineer** | Stage 4 (PRS method implementation/wrapping, especially PRS-CSx's MCMC configuration) |
| **Genetics domain expert** | Reviews the causal-variant/heritability simulation design for biological plausibility, leads Track B's descriptive interpretation, ensures the health-equity framing of the write-up is handled with appropriate care and accuracy |
| **Software Engineer** | Repository architecture, dual-language (Python/R) tooling integration, Stage 7 reporting |
| **DevOps/Cloud** | Environment reproducibility, compute provisioning for PRS-CSx |
| **Project owner** | Owns Definition of Done and the decision of how much emphasis the health-equity framing gets in external communications, given its sensitivity |

---

## 13. Deliverables Checklist

- [x] Public GitHub repository with complete documentation
- [x] `METHODS.md` documenting the full simulation architecture (causal variants, effect sizes, heritability, both scenarios) and every PRS method's configuration
- [x] `LICENSES.md` per §11 (license-confirmation sub-items above still open)
- [x] Working CI smoke test (Stage 1 only, per §5's CI-scope row; Stages 2+ run manually/on-demand as intended)
- [x] Full cross-ancestry predictive accuracy results table: 3 methods × 5 ancestries × 2 scenarios, before and after recalibration — all three methods (PRSice-2, PRS-CSx, LDpred2) completed 2026-08-02, see `METHODS.md` Stage 4/5
- [x] Recalibration effectiveness analysis, reported honestly including any residual gap
- [x] Track B descriptive cross-check results, clearly labeled as illustrative
- [x] ~1,500–2,500 word technical write-up — `reports/final_writeup/TECHNICAL_WRITEUP.md` (1,888 words)
- [x] One-page case-study-style PDF summary — `reports/final_writeup/case_study_summary.pdf`

---

## References

1. **Momin MM, Zhou X, Ahmed M, Hyppönen E, Benyamin B, Lee SH.** Cross-Ancestry Polygenic Prediction: Comparing Methods and Assessing Transferability Across Traits. *Genetic Epidemiology.* 2026;50:1–13. https://doi.org/10.1002/gepi.70029

2. **Ruan Y, Lin YF, Feng YCA, et al.** Improving polygenic prediction in ancestrally diverse populations (PRS-CSx). *Nature Genetics.* 2022;54:573–580. https://pubmed.ncbi.nlm.nih.gov/35513724/ — code: https://github.com/getian107/PRScsx

3. **Privé F, Arbel J, Vilhjálmsson BJ.** LDpred2: better, faster, stronger. *Bioinformatics.* 2020. Package: https://privefl.github.io/bigsnpr/

4. **Choi SW, O'Reilly PF.** PRSice-2: Polygenic Risk Score software for biobank-scale data. *GigaScience.* 2019. https://github.com/choishingwan/PRSice

5. **Lambert SA, Gil L, Jupp S, et al.** The Polygenic Score Catalog as an open database for reproducibility and systematic evaluation. *Nature Genetics.* 2021. https://pmc.ncbi.nlm.nih.gov/articles/PMC11165303/ — tools: https://github.com/PGScatalog/pgscatalog_utils

6. **1000 Genomes Project Consortium.** A global reference for human genetic variation. *Nature.* 2015. Data: https://www.internationalgenome.org

7. **Yang J, Lee SH, Goddard ME, Visscher PM.** GCTA: a tool for genome-wide complex trait analysis. *American Journal of Human Genetics.* 2011;88:76–82. https://cnsgenomics.com/software/gcta/

8. **Polygenic risk score translation across diverse populations.** *Frontiers in Cardiovascular Medicine.* 2026. https://doi.org/10.3389/fcvm.2026.1870807

9. **Polygenic risk score portability for common diseases across genetically diverse populations.** *Human Genomics.* 2024. https://doi.org/10.1186/s40246-024-00664-y

---

*This BUILD_PLAN.md is a living document. The health-equity framing of this project is genuinely sensitive — Stage 2's simulation design and the write-up's interpretation of results deserve unhurried review, and the team should be comfortable stating clearly what the synthetic simulation does and does not demonstrate about real clinical PRS deployment.*
