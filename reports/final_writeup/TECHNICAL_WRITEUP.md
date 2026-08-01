# PRS Cross-Ancestry Portability & Recalibration Benchmark: Technical Write-Up

**Status:** Track A Stages 1–6 (LDpred2 pending — see Limitations); Track B Stages B1–B3 complete.
Full parameters, tooling quirks, and per-stage run records: `METHODS.md`. Full project plan and
honesty/scope constraints: `docs/BUILD_PLAN.md`.

## 1. Motivation

Polygenic risk scores (PRS) built predominantly on European-ancestry GWAS are among the most
publicly discussed generalization-gap problems in genomics: a score trained on one ancestry
systematically underperforms in others, largely due to real differences in allele frequency and
linkage disequilibrium (LD) structure across populations, and sometimes due to genuine
gene-by-ancestry effect-size heterogeneity. Quantifying *how much* accuracy is lost, and how much
a recalibration step recovers, requires a validation with a known, controlled ground truth — real
phenotype data can't disentangle these effects, because we don't know the "true" genetic
architecture behind any real trait.

This project addresses that with a two-track design. **Track A** uses GCTA to simulate a
quantitative phenotype with a documented set of causal variants and heritability directly on real
1000 Genomes Phase 3 genotypes across five super-populations (AFR, AMR, EAS, EUR, SAS), giving a
genuine, controlled measurement of predictive-accuracy loss and recalibration recovery. **Track
B** applies a real, published PGS Catalog score to the same genotypes as a descriptive,
illustrative cross-check — explicitly non-validating, since 1000 Genomes carries no phenotype
data to validate against.

## 2. Methods

### 2.1 Data & scope

Real 1000 Genomes Phase 3 genotypes, QC'd (MAF≥0.01, genotype missingness≤0.05, HWE p≥1e-6 —
each filter applied *within* every super-population separately, keeping only variants passing in
all five) and restricted to biallelic SNPs for GCTA compatibility. This implementation is a
**chr21+chr22 smoke-test scope** (154,384 QC'd variants, 2,504 individuals: AFR 661, AMR 347, EAS
504, EUR 503, SAS 489) — the same scope this project's CI workflow uses, deliberately smaller than
a genome-wide run so every stage could be built, debugged, and verified against real computation
rather than mocked. The pipeline itself has no hardcoded chromosome limits; a genome-wide run is a
config change, not a rewrite.

### 2.2 Track A: synthetic ground truth

1. **Phenotype simulation (GCTA `--simu-qt`):** 300 causal variants (MAF 0.05–0.5, fixed seed),
   target heritability h²=0.5. **Scenario 1** applies identical effect sizes across all five
   super-populations, isolating the portability problem caused purely by real LD/allele-frequency
   differences. **Scenario 2** perturbs each super-population's effect sizes independently
   (`effect = base × N(1, 0.3)`, per-ancestry seeded), modeling genuine effect-size heterogeneity.
   Sanity-checked against a reconstructed genetic score: empirical heritability in EUR was 0.53
   against a target of 0.5 (n=503) — within sampling noise.
2. **Discovery GWAS:** simple, unadjusted per-SNP linear regression (`plink2 --glm`), restricted
   to an 80/20 split of EUR samples (402 discovery / 101 held out) — the held-out EUR subset is
   never used for GWAS or PRS construction, only for Stage 5's same-ancestry evaluation baseline.
3. **PRS construction**, three methods, each built only from the EUR discovery GWAS + EUR
   discovery genotypes (no evaluation-ancestry data ever touches construction):
   - **PRSice-2** (clumping+thresholding baseline): LD-clumped (r²<0.1, 250kb) on EUR discovery
     samples, scored at a fixed p-value threshold suite (5e-8 to 1) without target-phenotype
     regression, so the model doesn't depend on which ancestry it's later evaluated against.
   - **PRS-CSx** (Bayesian continuous-shrinkage): run in single-discovery-population mode (EUR
     only, matching this project's design — see §3), using the official 1000G EUR LD reference
     panel, restricted to the ~18% of our SNPs overlapping PRS-CSx's HapMap3 SNP set.
   - **LDpred2-auto** (Bayesian, `bigsnpr`): implemented, self-tuning (no validation phenotype
     needed), but **not completed end-to-end in this environment** — see Limitations.
4. **Cross-ancestry evaluation:** every scoring model applied via `plink2 --score` to held-out EUR
   plus the full AFR/AMR/EAS/SAS sample sets, R² computed against the true simulated phenotype
   (available because this is synthetic ground truth) with a 95% CI (Fisher z-transform).
5. **Recalibration:** empirical per-ancestry recentering/rescaling of each ancestry's PRS
   distribution to the EUR reference's mean/SD (BUILD_PLAN's documented fallback —
   `pgscatalog-ancestry-adjust` requires `pgsc_calc` pipeline outputs this project's custom
   scoring files don't produce).

### 2.3 Track B: real PGS Catalog cross-check

A real published height score (PGS000297, "GRS3290_Height", Xie et al. 2020, 3,290 variants,
GRCh37-harmonized) applied to all 2,504 individuals via `plink2 --score` — 52 variants overlap
this project's chr21+22 scope. Compared descriptively against published national height
statistics for representative countries per super-population (full citations in
`src/trackB_real_scores/stageB3_descriptive_comparison/published_height_reference.tsv`).

## 3. Results

### 3.1 Track A: does the EUR-derived model port?

Full results: `data/processed/evaluation/results.tsv` (85 rows: 2 methods × up to 8 PRSice-2
thresholds or 1 PRS-CSx model × 5 ancestries × 2 scenarios — no result is collapsed into a single
"best method" figure, per this project's reporting discipline).

**The real, reproducible finding across both scenarios:** PRSice-2's best EUR-tuned threshold
(p<0.01) clearly beats PRS-CSx in the EUR held-out set (R²=0.131 vs. 0.073, Scenario 1; 0.049 vs.
0.025, Scenario 2) — but that same PRSice-2 model loses to PRS-CSx in **all eight** non-EUR
(ancestry × scenario) comparisons:

| Scenario | Ancestry | PRSice-2 (best EUR threshold) | PRS-CSx |
|---|---|---|---|
| 1 (equal effects) | AFR | 0.0077 | **0.0275** |
| 1 | AMR | 0.0098 | **0.0594** |
| 1 | EAS | 0.0008 | **0.0165** |
| 1 | SAS | 0.0080 | **0.0279** |
| 2 (ancestry-varying) | AFR | 0.0081 | **0.0230** |
| 2 | AMR | 0.0009 | **0.0480** |
| 2 | EAS | 0.0090 | **0.0238** |
| 2 | SAS | 0.0161 | **0.0411** |

This matches the literature pattern this project is built around: a genome-wide Bayesian
shrinkage method ports across ancestries better than a simple clumping+thresholding baseline
tuned to one ancestry, even when the baseline looks better same-ancestry. Scenario 1 vs. Scenario
2 do not show a qualitatively different portability pattern here — both scenarios' cross-ancestry
R² are similarly depressed relative to EUR, suggesting (at this scope) the LD/allele-frequency
mechanism Scenario 1 isolates already accounts for most of the observed gap, with Scenario 2's
added effect-size heterogeneity not obviously worsening it further. This should be treated as a
preliminary observation, not a settled conclusion (see Limitations).

### 3.2 Recalibration: what it does and doesn't fix

Every one of the 85 (scenario × method × threshold × ancestry) rows in
`data/processed/recalibration/results.tsv` shows `r2_raw` and `r2_recalibrated` identical to four
decimal places, while `recal_mean`/`recal_sd` move to exactly match the EUR reference's after
recalibration (versus visibly differing before). This isn't a null result — it's the expected,
provable behavior of linear recentering/rescaling: R² depends only on correlation, which is
invariant to affine transforms of either variable. Empirical per-ancestry recalibration recovers
**100% of calibration and 0% of discriminative accuracy**, by mathematical necessity. This makes
concrete a distinction the PRS literature (and this project's own build plan) warns is often
conflated: a score whose mean/scale is corrected per ancestry is not the same as a score that
discriminates equally well per ancestry.

### 3.3 Track B: descriptive cross-check (not a validation)

**Every number in this section is descriptive/illustrative, never a validated accuracy claim** —
1000 Genomes has no phenotype data, so there's no ground truth here to validate a prediction
against; only Track A's synthetic experiment provides that. The real, honestly-reported finding:
**0 of 5** super-populations' mean PGS ranks match the published-height ranks (full table:
`data/processed/pgs_catalog/stageB3_descriptive_comparison.tsv`). AFR ranks highest by PGS but
2nd by published height; EUR ranks highest by published height but only 4th by PGS. Rather than
being evidence of a broken pipeline, this is a small-scale, real illustration of the exact PRS
cross-ancestry portability problem Track A's synthetic experiment exists to quantify with a
controlled ground truth: a European-ancestry-derived score, applied to other ancestries with no
ancestry-appropriate recalibration and using only 52 of 3,290 scoring variants (chr21+22 scope),
should not be expected to preserve true population-level phenotypic ordering.

## 4. Limitations

- **LDpred2 is implemented but not run end-to-end.** Three attempts (full chr21+22 panel;
  chr22-only with a reduced LD window and a 25-minute cap, which completed LD computation in ~14
  min but not MCMC sampling; chr22-only with a further-reduced MCMC budget and a 45-minute cap,
  externally killed with swap fully exhausted) all failed on this specific shared sandbox host
  (load average 130–150 sustained on 32 cores during testing). The script
  (`src/trackA_synthetic/stage4_prs_construction/ldpred2_wrapper.R`) is code-reviewed and should
  run as-is on a less contended host; the three-method comparison BUILD_PLAN.md calls for is
  currently two of three.
- **Smoke-test scope.** chr21+22 only (~6% of the autosomal genome), 300 causal variants, single
  simulation replicate. Effect sizes, portability-gap magnitudes, and the Scenario 1 vs. 2
  comparison should not be over-generalized to genome-wide behavior without a full run.
- **Small ancestry sample sizes** (347–661 per group) mean the confidence intervals in the raw
  results table are wide; the *direction* of the PRSice-2-vs-PRS-CSx portability pattern is
  consistent across all eight comparisons, but individual R² point estimates carry real
  uncertainty.
- **PRS-CSx ran in single-discovery-population mode**, not its namesake multi-ancestry coupling —
  this project deliberately feeds all three Stage 4 methods only the EUR-only discovery GWAS
  (mimicking real-world GWAS bias, per BUILD_PLAN.md's design), so PRS-CSx here is functionally
  equivalent to its predecessor, PRS-CS.
- **Track B's height score only covers 52 of 3,290 variants** at this project's genomic scope,
  and the published-height reference statistics are national averages for representative
  countries, not exact matches to 1000 Genomes' specific constituent source populations — both
  narrow how much weight the Track B finding can bear beyond "illustrative."

## 5. Conclusion

Within a chr21+22 smoke-test scope, this project reproduces — with a controlled synthetic ground
truth, not just real-data anecdote — the central, well-documented finding in current PRS
cross-ancestry literature: a EUR-derived score's advantage over other methods same-ancestry does
not carry over cross-ancestry, and a Bayesian genome-wide method (PRS-CSx) here consistently
outperforms a simple clumping+thresholding baseline (PRSice-2) in every non-EUR evaluation, even
while losing to it in the EUR held-out set. Separately, this project demonstrates numerically
(not just by assertion) that simple empirical recalibration fixes a score's per-ancestry
mean/scale completely while leaving its discriminative accuracy completely unchanged — a
distinction real-world PRS deployment discussions frequently conflate. The real-score Track B
cross-check adds a third, independent illustration of the same underlying phenomenon: an
uncorrected EUR-derived score's ranking across ancestries doesn't track true population-level
phenotype differences. None of these findings should be read as a genome-wide claim, a validated
clinical recommendation, or a definitive method ranking — see BUILD_PLAN.md §1.2 for this
project's explicit non-goals, and METHODS.md for the full parameter and run record behind every
number reported here.

## References

1. Momin MM, Zhou X, Ahmed M, Hyppönen E, Benyamin B, Lee SH. Cross-Ancestry Polygenic
   Prediction: Comparing Methods and Assessing Transferability Across Traits. *Genetic
   Epidemiology.* 2026;50:1–13. https://doi.org/10.1002/gepi.70029
2. Ruan Y, Lin YF, Feng YCA, et al. Improving polygenic prediction in ancestrally diverse
   populations (PRS-CSx). *Nature Genetics.* 2022;54:573–580.
   https://github.com/getian107/PRScsx
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
