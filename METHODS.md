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
- **Tooling quirk:** the bed/fam's implicit 6th-column phenotype is `-9` (missing) for every
  sample; a `--pheno` file whose column is also named `PHENO1` collides with it
  (`Duplicate phenotype/covariate ID 'PHENO1'`), and forcing `-9` to be read as a real value
  (`--no-input-missing-phenotype`) makes `--glm` choke on that now-constant phenotype instead.
  Fixed by naming our column `SYNTH_PHENO` and passing `--pheno-name SYNTH_PHENO` to restrict
  `--glm` to only our simulated phenotype.
- **Verified end-to-end on real data** (chr21+chr22 smoke-test subset, 2026-08-01): both
  scenarios ran successfully (503 EUR samples, 154,384 SNPs tested, all 300 causal variants
  present in the output). Sanity-checked causal-variant enrichment against background: causal
  variants are nominally significant (p<0.05) 13.3% (Scenario 1) / 15.0% (Scenario 2) of the
  time vs. 5.6-5.7% for background SNPs (~2.5x enrichment). This is modest, not dramatic —
  expected and correct for a polygenic architecture where h²=0.5 is spread across 300 causal
  variants (~0.17% heritability each), well below what n=503 can detect per-SNP at strict
  significance. This is the same reason Bayesian PRS methods that use genome-wide signal
  (LDpred2, PRS-CSx) are expected to outperform a top-hits-only baseline (PRSice-2) on this
  kind of trait — see BUILD_PLAN.md §4's PRSice-2/LDpred2 rows.
- Outputs land in `data/processed/gwas/scenario{1,2}_*/eur_discovery_gwas.SYNTH_PHENO.glm.linear`
  (standard plink2 `--glm` summary-stats format: CHROM, POS, ID, REF, ALT, A1, A1_FREQ, BETA,
  SE, T_STAT, P, ...).

## Stages 4-6 (Track A) / B1-B3 (Track B)

Not yet implemented. See `docs/BUILD_PLAN.md` §6 for the full stage detail this section will
document once each stage lands.
