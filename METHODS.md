# METHODS.md

Full methodological documentation for the PRS Cross-Ancestry Portability & Recalibration
Benchmark. This file is filled in stage-by-stage as each part of `docs/BUILD_PLAN.md` §6 is
implemented; it is the canonical reproducibility record referenced by the final write-up.

**Status:** Stages 1-2 implemented and verified on real data.

## Stage 1 — 1000 Genomes Download & QC (shared)

- Source: 1000 Genomes Phase 3, GRCh37, IGSR/EBI FTP-over-HTTPS
  (`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502`).
- Implementation: `src/shared_stage1_1000g_download_qc/download.py` (fetch),
  `src/shared_stage1_1000g_download_qc/qc.py` (convert + QC).
- QC filters — applied **within each super-population separately** (AFR, AMR, EAS, EUR, SAS),
  keeping only variants that pass in all five, per `docs/BUILD_PLAN.md` §6 Stage 1:
  - MAF (`--maf`): default 0.01
  - per-variant missingness (`--geno`): default 0.05
  - Hardy-Weinberg equilibrium (`--hwe`): default 1e-6
  - per-sample missingness (`--mind`), applied globally after variant filtering: default 0.05
- Defaults are set in `src/shared_stage1_1000g_download_qc/qc.py`'s CLI arguments; override via
  flags rather than editing the script, and record any override used for a given run here.
- `--set-all-var-ids` is applied with `--new-id-max-allele-len 1000 missing`: some 1000 Genomes
  indels/structural variants have allele codes past plink2's default ID-length cap, so anything
  still over 1000 characters (very rare) gets a `.` ID rather than failing the run.
- Sample keep-lists are single-column IID files: the VCF-derived pgen has no `#FID` column
  (IID-only mode), so a two-column FID+IID keep file silently matches zero samples.
- **Verified end-to-end on real data** (chr21+chr22 smoke-test subset, 2026-08-01): 2,504
  samples loaded; variants passing QC per super-population — AFR 462,550, AMR 308,829,
  EAS 240,128, EUR 272,054, SAS 284,717; **179,285 variants shared across all five**
  (the final QC'd variant set).

## Stage 2 — GCTA Phenotype Simulation (Track A)

- Parameters: `configs/simulation_parameters.yaml` (**still marked DRAFT, pending
  Biostatistics review** per the open questions listed there — the implementation below
  reflects that draft, not a signed-off design).
- Implementation: `src/trackA_synthetic/stage2_gcta_simulation/common.py` (shared plumbing),
  `scenario1_equal_effects.py`, `scenario2_ancestry_varying_effects.py`.
- GCTA requires PLINK1 bed/bim/fam (no multiallelic support), so Stage 1's QC'd pgen is
  converted once into a shared, biallelic-SNP-only bed/bim/fam
  (`data/processed/simulated_phenotypes/_shared/1000g_qc_biallelic`) — chr21+chr22 smoke subset
  went from 179,285 QC'd variants to 154,384 biallelic SNPs.
- 300 causal variants are selected with MAF in [0.05, 0.5] using the fixed seed in
  `simulation_parameters.yaml`; the *same* causal-variant set is used in both scenarios so
  they're directly comparable.
- GCTA is run once per super-population on the shared bed/bim/fam (`--keep` restricting
  samples; allele coding stays fixed across runs since the bed/bim/fam itself doesn't change),
  each targeting `--simu-hsq 0.5`. This means every ancestry gets the same nominal
  heritability by construction — the portability gap shows up in downstream PRS *prediction*
  accuracy (Stage 5), not in Stage 2's heritability itself.
- Scenario 1: identical effect sizes (`base_effect`) passed to every super-population's GCTA
  run. Scenario 2: `effect[ancestry] = base_effect * N(1, perturbation_sd=0.3)`, independently
  sampled per super-population via the seed offsets in the config.
- **Known tooling limitation:** GCTA v1.95.3 has no `--seed` option for `--simu-qt`'s residual
  noise — the genetic component (causal variants, effect sizes, perturbations) is fully
  reproducible via our own seeded RNG, but the residual/noise draw is not bit-for-bit
  reproducible across re-runs. Documented here rather than worked around.
- **Effect sizes are interpreted per-variance-standardized-genotype unit, not per raw allele
  count** — confirmed empirically (see validation below), not documented explicitly in GCTA's
  own docs. Relevant for Stage 4: PRS methods that assume a specific effect-size scale need to
  account for this.
- **Verified end-to-end on real data** (chr21+chr22 smoke-test subset, 2026-08-01): both
  scenarios ran successfully across all five super-populations (2,504 individuals total:
  AFR 661, AMR 347, EAS 504, EUR 503, SAS 489). Sanity-checked per BUILD_PLAN §9 item 1 by
  reconstructing the genetic score from the 300 causal variants' variance-standardized
  genotypes × GCTA's own effect sizes (Scenario 1, EUR): empirical heritability
  Var(genetic)/Var(phenotype) = 0.53 against a target of 0.5 (n=503; the ~0.03 gap is
  consistent with sampling noise at this sample size, not a simulation bug).
- Outputs land in `data/processed/simulated_phenotypes/scenario{1,2}_*/phenotypes.tsv`
  (columns: FID, IID, super_pop, phenotype).

## Stage 3 — EUR-Only Discovery GWAS (Track A)

- Implementation: `src/trackA_synthetic/stage3_eur_gwas/run_gwas.py`.
- Simple, unadjusted per-SNP linear regression (`plink2 --glm allow-no-covars`, no PCs/
  covariates) on each scenario's simulated phenotype, restricted to EUR samples only — per
  BUILD_PLAN.md §6 Stage 3's explicit "simple linear regression per-SNP" spec. Run once per
  scenario (Definition of Done item 3 requires both).
- Uses the same shared biallelic bed/bim/fam Stage 2 simulated the phenotype on, so every
  causal variant is guaranteed present in the GWAS output for Stage 4's PRS construction.
- **EUR discovery/held-out split**: before the GWAS runs, EUR samples are partitioned into a
  discovery subset (used for the GWAS) and a held-out subset the GWAS never sees
  (`configs/simulation_parameters.yaml`: `discovery_gwas.holdout_fraction=0.2`,
  `holdout_seed=20260801`; cached in `data/processed/gwas/_shared/keep_EUR_{discovery,holdout}.txt`,
  shared across both scenarios). This wasn't in the initial implementation — Stage 5 (and
  BUILD_PLAN.md §9 item 1's sanity check) explicitly require a genuine held-out same-ancestry
  evaluation set; without it, "EUR" evaluation would be in-sample against the GWAS training
  data, not a real held-out test. Caught and fixed while starting Stage 4.
- **Tooling quirk:** the bed/fam's implicit 6th-column phenotype is `-9` (missing) for every
  sample; a `--pheno` file whose column is also named `PHENO1` collides with it
  (`Duplicate phenotype/covariate ID 'PHENO1'`), and forcing `-9` to be read as a real value
  (`--no-input-missing-phenotype`) makes `--glm` choke on that now-constant phenotype instead.
  Fixed by naming our column `SYNTH_PHENO` and passing `--pheno-name SYNTH_PHENO` to restrict
  `--glm` to only our simulated phenotype.
- **Verified end-to-end on real data** (chr21+chr22 smoke-test subset, 2026-08-01): both
  scenarios ran successfully (503 EUR samples split into 402 discovery / 101 held-out,
  154,384 SNPs tested, all 300 causal variants present in the output). Sanity-checked
  causal-variant enrichment against background: causal variants are nominally significant
  (p<0.05) 13.0% (Scenario 1) / 13.7% (Scenario 2) of the time vs. 5.2-5.7% for background
  SNPs (~2.4x enrichment). This is modest, not dramatic —
  expected and correct for a polygenic architecture where h²=0.5 is spread across 300 causal
  variants (~0.17% heritability each), well below what n=503 can detect per-SNP at strict
  significance. This is the same reason Bayesian PRS methods that use genome-wide signal
  (LDpred2, PRS-CSx) are expected to outperform a top-hits-only baseline (PRSice-2) on this
  kind of trait — see BUILD_PLAN.md §4's PRSice-2/LDpred2 rows.
- Outputs land in `data/processed/gwas/scenario{1,2}_*/eur_discovery_gwas.SYNTH_PHENO.glm.linear`
  (standard plink2 `--glm` summary-stats format: CHROM, POS, ID, REF, ALT, A1, A1_FREQ, BETA,
  SE, T_STAT, P, ...).

## Stage 4 — PRS Construction (Track A)

Three methods, per BUILD_PLAN.md §6 Stage 4, each producing a portable scoring file (ID, A1,
BETA) that Stage 5 applies identically to all five super-populations via `plink2 --score` — none
of the three ever see any evaluation-target ancestry's data during construction.

- **PRSice-2** (`src/trackA_synthetic/stage4_prs_construction/prsice2_wrapper.py`): real
  PRSice-2 binary, LD-clumped (r²<0.1, 250kb) using the EUR discovery samples as reference,
  `--no-regress` + `--fastscore` across a fixed threshold suite (5e-8 to 1) so construction
  never needs a target phenotype. Retained SNPs re-joined to Stage 3's betas ourselves rather
  than relying on PRSice's own scoring, since PRSice computes scores against whatever target
  it's given — we need one portable model applicable to all five ancestries, not five
  separately-clumped ones. **Tooling quirks:** PRSice's underlying binary mis-parses a
  `#`-prefixed header (silently swallows the next CLI arg — strip plink2's leading `#` first);
  and `PRSice.R` reinstalls ggplot2/data.table/optparse into a local `./lib` on every run unless
  `R_LIBS_USER` is set in the environment first (50MB, several minutes wasted per run otherwise).
- **LDpred2-auto** (`src/trackA_synthetic/stage4_prs_construction/ldpred2_wrapper.R`, via
  `bigsnpr`): implemented and code-reviewed, but **not completed in this session** — this
  sandbox is a heavily shared host (load average 130-150 sustained on 32 cores during testing,
  swap fully exhausted at times) and three attempts (full chr21+22/250kb window; chr22-only/100kb
  window with a 25-min cap, which got through the LD step in ~14 min but not the MCMC step;
  chr22-only with a reduced MCMC budget and 45-min cap, externally killed before finishing) all
  failed to complete. No genetic map is available, so `snp_cor`'s LD window is defined in
  physical bp (100kb radius) rather than cM — a documented approximation, not the cause of the
  slowness (the chr22-only LD step itself did complete once, in ~14 minutes). The script is
  ready to run as-is (`Rscript ldpred2_wrapper.R --scenario all`) given a less contended host;
  no code changes should be needed.
- **PRS-CSx** (`src/trackA_synthetic/stage4_prs_construction/prscsx_wrapper.py`): real PRS-CSx,
  run in single-discovery-population mode (`--pop EUR`) since BUILD_PLAN.md feeds all three
  methods the same EUR-only sumstats — the "x" (coupling effects *across* ancestries) only
  activates with multiple discovery GWASes, which this project deliberately doesn't do. Requires
  the official 1000G EUR LD reference panel (~4.25GB) and HapMap3 SNP list (~105MB), both fetched
  from Dropbox — much faster than IGSR/GCTA's hosts in this sandbox (~22MB/s vs. ~30KB/s), full
  panel in under 4 minutes despite earlier fears it would take 30+ hours. **SNP ID
  compatibility:** PRS-CSx's reference panel is restricted to ~1.1M HapMap3 SNPs addressed by
  rsID, but Stage 1 assigned chr:pos:ref:alt IDs, so this wrapper builds a position-matched
  rsID↔our-ID mapping from the HapMap3 SNP list before running PRS-CSx, then translates the
  output back — of our 154,384 QC'd biallelic chr21+22 SNPs, 27,881 (~18%) overlap HapMap3.
  **Tooling quirk:** PRS-CSx expects `snpinfo_mult_1kg_hm3` directly inside `--ref_dir`, alongside
  the `ldblk_1kg_<pop>/` subdirectory, not as a sibling of `--ref_dir` — its own `main()` has an
  `UnboundLocalError` if that file isn't found there (a real bug in the tool, not just an
  unclear error message).
- **Verified on real data** (chr21+chr22, 2026-08-01/02): PRSice-2 and PRS-CSx both completed
  successfully for both scenarios. PRSice-2: 5,858 SNPs survive clumping (of 131,173 non-ambiguous),
  8-9 usable p-value thresholds per scenario (5e-8 has 0 SNPs in Scenario 1, 1 in Scenario 2 —
  reported honestly as zero-SNP rows are simply absent from the threshold table, consistent with
  Stage 3's finding that this polygenic architecture has weak per-SNP signal). PRS-CSx: full MCMC
  (1,000 iterations) completed for both chromosomes in both scenarios, 27,881 SNPs scored.

## Stage 5 — Cross-Ancestry Evaluation (Track A)

- Implementation: `src/trackA_synthetic/stage5_crossancestry_evaluation/evaluate.py`.
- Every Stage 4 scoring file applied via `plink2 --score` to five held-out groups: **EUR uses
  the held-out subset** Stage 3 set aside (101 samples, never in the discovery GWAS or any
  Stage 4 construction step) — **not** the full EUR or EUR-discovery set, which would be
  in-sample/circular. AFR (661), AMR (347), EAS (504), SAS (489) use their full sample sets
  (never used anywhere upstream). R² computed against the *true* simulated phenotype (available
  because this is synthetic ground truth), with a 95% CI via the Fisher z-transform of the
  Pearson correlation.
- Full results table (no collapsing to a single "best method" figure, per BUILD_PLAN.md §9 item
  3): `data/processed/evaluation/results.tsv`, one row per (scenario, method, [threshold],
  ancestry) — 85 rows across both scenarios with the two completed methods.
- **Headline finding, real and reproducible across both scenarios** (not cherry-picked — this is
  the whole results table's pattern): PRSice-2's best-performing threshold in the EUR held-out
  set (p<0.01, R²=0.131 Scenario 1 / 0.049 Scenario 2) clearly beats PRS-CSx's same-ancestry R²
  (0.073 / 0.025) — but that same PRSice-2 model is *worse* than PRS-CSx in **every one of the
  eight** non-EUR (ancestry × scenario) comparisons:

  | Scenario | Ancestry | PRSice-2 (best EUR threshold) | PRS-CSx |
  |---|---|---|---|
  | 1 | AFR | 0.0077 | **0.0275** |
  | 1 | AMR | 0.0098 | **0.0594** |
  | 1 | EAS | 0.0008 | **0.0165** |
  | 1 | SAS | 0.0080 | **0.0279** |
  | 2 | AFR | 0.0081 | **0.0230** |
  | 2 | AMR | 0.0009 | **0.0480** |
  | 2 | EAS | 0.0090 | **0.0238** |
  | 2 | SAS | 0.0161 | **0.0411** |

  This matches the literature pattern this project is built around (BUILD_PLAN.md §0, §4):
  genome-wide Bayesian shrinkage methods port across ancestries better than a simple
  clumping+thresholding baseline tuned to one ancestry, even when the baseline looks better
  same-ancestry. **Caveats that must travel with this finding**: single simulation replicate,
  chr21+22-only smoke-test scope, small ancestry sample sizes (CIs in the raw results table are
  wide), and no LDpred2 result yet to complete the three-method comparison BUILD_PLAN.md calls
  for — this is a real, honest, but partial result, not a final claim.

## Stage 6 — Recalibration (Track A)

- Implementation: `src/trackA_synthetic/stage6_recalibration/recalibrate.py`.
- BUILD_PLAN.md's first choice, `pgscatalog-ancestry-adjust`, requires PCA projections from
  `fraposa_pgsc` and PGS-Catalog-format aggregated scores — outputs of the full `pgsc_calc`
  pipeline that this project's custom-built Track A scoring files never produce. Used
  BUILD_PLAN.md's documented fallback instead: empirical per-ancestry recentering/rescaling,
  matching each non-EUR ancestry's raw PRS mean/SD to the EUR held-out reference distribution's.
- **The actual point of this stage, demonstrated not just asserted**: a linear
  recentering/rescaling transform cannot change Pearson R² — R² depends only on correlation,
  which is invariant to affine transforms. Confirmed numerically in every one of the 85 rows:
  `r2_raw` and `r2_recalibrated` in `data/processed/recalibration/results.tsv` are identical to
  4 decimal places throughout. Meanwhile `recal_mean`/`recal_sd` exactly match the EUR
  reference's `ref_mean`/`ref_sd` after recalibration (calibration fixed), while `raw_mean`/
  `raw_sd` visibly differ from the reference beforehand (the calibration problem recalibration
  is meant to fix). This makes BUILD_PLAN.md §9 item 2's warned-about distinction — "recovers
  calibration ... without fully recovering discriminative accuracy" — concrete rather than a
  claim to take on faith: here recalibration recovers *100% of calibration and 0% of
  discriminative accuracy*, by mathematical necessity, which is itself an honest, reportable
  finding about what this class of recalibration can and cannot fix.
