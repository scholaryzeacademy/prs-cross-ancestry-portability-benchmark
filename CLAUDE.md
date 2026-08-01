# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project

PRS Cross-Ancestry Portability & Recalibration Benchmark — a portfolio project quantifying how much
polygenic risk score (PRS) predictive accuracy is lost when a score trained on European-ancestry GWAS
is applied to other ancestries, and how much a recalibration step recovers. Full plan: `docs/BUILD_PLAN.md`
(read it before making architectural decisions — this file only summarizes it).

**Current status:** Planning/early build. Only `docs/BUILD_PLAN.md` and `scripts/download_data.py` exist
so far; the `src/`, `notebooks/`, `configs/`, `tests/`, `reports/` layout in §7 of the build plan is the
target structure, not yet built.

## Two-track structure

- **Track A (synthetic ground truth, the primary result):** real 1000 Genomes Phase 3 genotypes (AFR,
  AMR, EAS, EUR, SAS) + GCTA-simulated phenotypes with known causal variants/heritability, in two
  scenarios (equal effect sizes across ancestries vs. ancestry-varying effect sizes). A EUR-only
  discovery GWAS feeds three PRS methods (PRSice-2, LDpred2, PRS-CSx), evaluated for R² across all five
  super-populations, then recalibrated.
- **Track B (real published scores, descriptive only):** real PGS Catalog scores (e.g. height, BMI)
  computed on the same 1000 Genomes genotypes via `pgsc_calc`/`pgscatalog-utils`. 1000 Genomes has no
  real phenotypes, so this track is explicitly non-validating — always label it descriptive/illustrative,
  never as an accuracy result, wherever it's discussed in code comments, notebooks, or write-ups.

## Non-negotiable framing (from BUILD_PLAN.md §1.2, §9)

- Never claim clinically deployable recommendations — this is a methodology demonstration on public data.
- Never rank PRS methods as universally "best." Report per-method × per-ancestry × per-scenario results
  as full tables; current literature (Momin et al. 2026) is explicit that no method wins universally.
- Distinguish recalibration's effect on *calibration* (mean/scale per ancestry) from its effect on
  *discriminative accuracy* (R²) — conflating the two is a flagged analytical error.
- Re-state Track B's descriptive-only status every place it's referenced, not just once.
- The health-equity framing is sensitive; treat Stage 2 (simulation design) and result interpretation
  with care rather than defaulting to punchy claims.

## Data & tooling

- Data sources: 1000 Genomes Phase 3 (IGSR/EBI), PGS Catalog (`pgscatalog.org`). No credentialed/biobank
  access (UK Biobank, All of Us) is in scope — deliberately public-data-only.
- Tools: PLINK/PLINK2 (QC), GCTA (`--simu-qt`/`--simu-cc` phenotype simulation), PRSice-2
  (clumping+thresholding baseline), LDpred2 via R `bigsnpr`, PRS-CSx (MCMC, CPU-bound — no GPU needed but
  budget real wall-clock time), `pgscatalog-utils`/`pgsc_calc` (PGS Catalog tooling, incl.
  `pgscatalog-ancestry-adjust` for recalibration).
- Stack spans Python 3.10+ and R 4.x; expect both `environment.yml` and an `renv`/CRAN lockfile as the
  project matures.
- `scripts/download_data.py` is a stdlib-only fetcher (no third-party deps; `requirements.txt` documents
  this) for genotypes, PGS scoring files, and tool binaries/source. Key usage:
  ```
  python3 scripts/download_data.py --list                                   # show resolved catalog
  python3 scripts/download_data.py                                          # light smoke-test set (chr21/22, no LD panels)
  python3 scripts/download_data.py --categories all --chromosomes all --prscsx-ld eur
  ```
  Downloads are idempotent (size/Content-Length checked, skip-if-present, `--force` to redo) and land in
  `data/`, which is gitignored — never commit downloaded data or tool binaries.

## Conventions

- Every PRS method construction/evaluation step should include a naive baseline for comparison (PRSice-2
  fills that role among the three methods) — this discipline runs through the whole portfolio.
- Document every GCTA simulation parameter (causal variants, effect sizes, heritability, per-ancestry
  perturbation scheme) in `configs/simulation_parameters.yaml` once it exists — the simulation must be
  independently reproducible from that file alone.
- CI should stay light: a small chromosome-arm subset with a handful of simulated causal variants. Full
  genome-wide PRS-CSx/LDpred2 runs are too heavy for CI and are run manually/on-demand.
