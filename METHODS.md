# METHODS.md

Full methodological documentation for the PRS Cross-Ancestry Portability & Recalibration
Benchmark. This file is filled in stage-by-stage as each part of `docs/BUILD_PLAN.md` §6 is
implemented; it is the canonical reproducibility record referenced by the final write-up.

**Status:** Phase 0 (setup). Only Stage 1 is implemented so far.

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

Draft architecture (pending Biostatistics review): `configs/simulation_parameters.yaml`.
Not yet implemented.

## Stages 3-6 (Track A) / B1-B3 (Track B)

Not yet implemented. See `docs/BUILD_PLAN.md` §6 for the full stage detail this section will
document once each stage lands.
