# PRS Cross-Ancestry Portability & Recalibration Benchmark

Quantifies how much predictive accuracy a polygenic risk score (PRS) trained on
European-ancestry GWAS loses when applied to other ancestries — using real multi-ancestry
1000 Genomes genotypes with a controlled, synthetic (GCTA-simulated) phenotype as ground
truth — and how much an ancestry-recalibration step recovers.

Full plan, rationale, and honesty/scope constraints: **[`docs/BUILD_PLAN.md`](docs/BUILD_PLAN.md)**.
Repo-wide working agreements for contributors (human or AI): **[`CLAUDE.md`](CLAUDE.md)**.

> **Status: Stages 1-6 implemented.** Two tracks are planned — **Track A** (synthetic ground
> truth via GCTA simulation + three PRS methods) and **Track B** (real PGS Catalog scores,
> descriptive only). Stages 1-3 (1000G download/QC, GCTA simulation, EUR discovery GWAS) and
> 5-6 (cross-ancestry evaluation, recalibration) are implemented and verified on real data.
> Stage 4 (PRS construction) has PRSice-2 and PRS-CSx complete; LDpred2 is implemented but not
> yet run end-to-end (blocked by this sandbox's resource contention, not a code issue — see
> METHODS.md). Track B is not yet built.

## Quickstart

```bash
# 1. Create the Python/R environment
conda env create -f environment.yml
conda activate prs-portability

# 2. Check (and optionally install) GCTA, PRSice-2, PRS-CSx, plink2, pgscatalog-utils, bigsnpr
python3 scripts/verify_tools.py --install            # add --with-r to also attempt bigsnpr (slow)

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
```

`download.py --chromosomes all` fetches the full autosomal panel (tens of GB) instead of the
default smoke-test subset — do that deliberately, not by default.

## Repository layout

```
src/shared_stage1_1000g_download_qc/   Stage 1: download + QC 1000 Genomes (shared by both tracks)
src/trackA_synthetic/stage2_gcta_simulation/  Stage 2: GCTA phenotype simulation (both scenarios)
src/trackA_synthetic/stage3_eur_gwas/  Stage 3: EUR-only discovery GWAS (both scenarios)
src/trackA_synthetic/stage4_prs_construction/  Stage 4: PRSice-2, LDpred2, PRS-CSx wrappers
src/trackA_synthetic/stage5_crossancestry_evaluation/  Stage 5: apply models to all 5 ancestries, compute R²
src/trackA_synthetic/stage6_recalibration/  Stage 6: empirical per-ancestry recalibration
src/trackA_synthetic/                  Track A: GCTA simulation -> EUR GWAS -> PRS construction -> evaluation -> recalibration
src/trackB_real_scores/                Track B: real PGS Catalog scores applied to the same genotypes (descriptive only)
configs/simulation_parameters.yaml     Stage 2 simulation architecture (causal variants, effect sizes, heritability) — DRAFT, pending Biostatistics review
scripts/verify_tools.py                Check/install GCTA, PRSice-2, PRS-CSx, plink2, pgscatalog-utils, bigsnpr
notebooks/, reports/                   Result notebooks and the final write-up
```

See `docs/BUILD_PLAN.md` §7 for the full target layout and §12 for stage ownership.

## Data & tooling

All data sources are fully public, no-credential-application-required (1000 Genomes Phase 3,
PGS Catalog) — see `docs/BUILD_PLAN.md` §1.2 and §4. Downloaded data/tools are gitignored
(`data/`, `tools/`); re-fetch them with the scripts above rather than committing them.

## Reporting discipline

This project's honesty constraints (no universal "best method" claims, Track B is descriptive
only, calibration vs. discriminative-accuracy are reported separately) are load-bearing — see
`docs/BUILD_PLAN.md` §9 and `CLAUDE.md` before adding new analysis or write-up content.
