# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project

PRS Cross-Ancestry Portability & Recalibration Benchmark — a portfolio project quantifying how much
polygenic risk score (PRS) predictive accuracy is lost when a score trained on European-ancestry GWAS
is applied to other ancestries, and how much a recalibration step recovers. Full plan: `docs/BUILD_PLAN.md`
(read it before making architectural decisions — this file only summarizes it).

**Current status:** All 7 stages implemented and verified on real data (chr21+chr22 smoke-test
subset). Stage 1 (`src/shared_stage1_1000g_download_qc/`) downloads/QCs 1000 Genomes. Stage 2
(`src/trackA_synthetic/stage2_gcta_simulation/`) runs the GCTA phenotype simulation for both scenarios.
Stage 3 (`src/trackA_synthetic/stage3_eur_gwas/`) runs the EUR-only discovery GWAS for both scenarios.
Stage 4 (`src/trackA_synthetic/stage4_prs_construction/`) has all three methods verified working —
PRSice-2, PRS-CSx, and (as of 2026-08-02) LDpred2, which needed two chr22-scoped reruns to get past
this sandbox's resource contention (not a code issue — see METHODS.md for the full run history
before re-attempting anything similar). Stage 5 (`.../stage5_crossancestry_evaluation/`) and Stage 6
(`.../stage6_recalibration/`) are implemented and verified with the full 3-method × 5-ancestry ×
2-scenario results table. Stage 7 (`reports/final_writeup/`) is the combined technical write-up and
case-study PDF. Track B (`src/trackB_real_scores/`, all three stages: download real PGS Catalog
score, compute on all 2,504 samples, compare to published height stats) is implemented and
verified — remember its results are descriptive/illustrative only, never a validated accuracy claim
(see the Non-negotiable framing section below). See `METHODS.md` for the run record and non-obvious
fixes/quirks found while implementing each stage — read it before re-deriving something that
was already debugged there.

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
- `scripts/verify_tools.py --install` checks/installs GCTA, PRSice-2, PRS-CSx, plink2, pgscatalog-utils,
  and bigsnpr without sudo (static binaries into `tools/`, `pip install --user`, R via a user-local
  library — the system R library isn't writable here). Add `--with-r` for bigsnpr (slow: compiles from
  source, and needs the `~/.R/Makevars` gfortran-linker fix documented in this project's memory).
- `src/shared_stage1_1000g_download_qc/download.py` fetches 1000 Genomes VCFs + panel (stdlib-only,
  idempotent, resumable — this sandbox's bandwidth to some hosts is only ~20-40 KB/s, so resumable
  downloads matter). Downloaded data/tools land in `data/` and `tools/`, both gitignored — never commit
  them.

## Conventions

- Every PRS method construction/evaluation step should include a naive baseline for comparison (PRSice-2
  fills that role among the three methods) — this discipline runs through the whole portfolio.
- Every GCTA simulation parameter (causal variants, effect sizes, heritability, per-ancestry
  perturbation scheme) lives in `configs/simulation_parameters.yaml` — it's still marked DRAFT pending
  Biostatistics review (see the open questions listed in that file), so treat its values as provisional,
  not settled, until that review happens.
- CI should stay light: a small chromosome-arm subset with a handful of simulated causal variants. Full
  genome-wide PRS-CSx/LDpred2 runs are too heavy for CI and are run manually/on-demand.
