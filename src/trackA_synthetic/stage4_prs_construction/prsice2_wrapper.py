#!/usr/bin/env python3
"""
prsice2_wrapper.py — Stage 4, PRSice-2 clumping+thresholding baseline (BUILD_PLAN.md §6 Stage 4).

Runs real PRSice-2 (not a hand-rolled substitute) using Stage 3's EUR discovery-only GWAS as
the base and the EUR discovery samples as the LD-clumping reference — never the held-out EUR
subset or the other four ancestries, since this baseline model must be constructed independently
of any evaluation target for Stage 5's portability comparison to be meaningful.

Design: rather than relying on PRSice-2's own --pheno/regression-based "best threshold" search
(which would require a target phenotype and couple construction to a specific evaluation
target), this runs PRSice-2 with --no-regress across a fixed, standard suite of p-value
thresholds (--fastscore --bar-levels) and --print-snp to get the post-clumping SNP set. Each
threshold's retained SNPs are then joined back against the GWAS betas to produce a portable
scoring file (SNP, A1, BETA) that Stage 5 applies identically to all five super-populations via
`plink2 --score` — the actual "does a EUR-derived model port to other ancestries" test.

Usage:
    python3 prsice2_wrapper.py                  # both scenarios
    python3 prsice2_wrapper.py --scenario 1
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
SHARED_DIR = REPO_ROOT / "data" / "processed" / "simulated_phenotypes" / "_shared"
BED_PREFIX = SHARED_DIR / "1000g_qc_biallelic"
EUR_DISCOVERY_KEEP_FILE = REPO_ROOT / "data" / "processed" / "gwas" / "_shared" / "keep_EUR_discovery.txt"
GWAS_ROOT = REPO_ROOT / "data" / "processed" / "gwas"
PRS_MODELS_ROOT = REPO_ROOT / "data" / "processed" / "prs_models"

PRSICE_R = REPO_ROOT / "tools" / "prsice2" / "PRSice.R"
PRSICE_BIN = REPO_ROOT / "tools" / "prsice2" / "PRSice_linux"
R_USER_LIB = Path.home() / "R" / "library"  # avoids PRSice.R reinstalling ggplot2/data.table/optparse

# Standard PRS-CT threshold suite. 5e-8/1e-6 are included for completeness even though this
# simulation's polygenic architecture (h2=0.5 spread across 300 causal variants) means no SNP
# is expected to reach those thresholds — see METHODS.md Stage 3 for why, and Stage 5 reports
# honestly on which thresholds end up with zero SNPs rather than silently dropping them.
BAR_LEVELS = ["5e-8", "1e-6", "1e-4", "1e-3", "0.01", "0.05", "0.1", "0.5", "1"]
CLUMP_R2 = "0.1"
CLUMP_KB = "250"

SCENARIOS = {
    "1": "scenario1_equal_effects",
    "2": "scenario2_ancestry_varying_effects",
}


def run(cmd: list, cwd: Path = None) -> None:
    logging.info("+ %s", " ".join(str(c) for c in cmd))
    env = os.environ.copy()
    env["R_LIBS_USER"] = str(R_USER_LIB)
    subprocess.run([str(c) for c in cmd], check=True, cwd=cwd, env=env)


def clean_sumstats(gwas_glm_linear: Path, dest: Path) -> None:
    """PRSice-2's underlying binary mis-parses a '#'-prefixed header (silently swallows the
    next CLI argument), so strip plink2's leading '#' from the CHROM column name."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(gwas_glm_linear) as fin, open(dest, "w") as fout:
        header = fin.readline()
        fout.write(header.lstrip("#"))
        fout.writelines(fin)


def load_sumstats(sumstats_path: Path) -> dict:
    with open(sumstats_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return {row["ID"]: (row["A1"], row["BETA"]) for row in reader}


def parse_threshold_counts(prsice_summary: Path) -> list[tuple[str, int]]:
    """Read <out>.prsice: columns Pheno, Set, Threshold, Num_SNP."""
    thresholds = []
    with open(prsice_summary) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            thresholds.append((row["Threshold"], int(row["Num_SNP"])))
    return thresholds


def parse_clumped_snps(snp_file: Path) -> dict:
    """Read <out>.snp (post-clumping SNP set): columns CHR, SNP, BP, P, <geneset columns>."""
    snp_p = {}
    with open(snp_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            snp_p[row["SNP"]] = float(row["P"])
    return snp_p


def write_scoring_file(snp_ids: list, sumstats: dict, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w") as f:
        f.write("ID\tA1\tBETA\n")
        for snp_id in snp_ids:
            a1, beta = sumstats[snp_id]
            f.write(f"{snp_id}\t{a1}\t{beta}\n")


def run_prsice2_for_scenario(scenario_name: str) -> list:
    gwas_glm_linear = GWAS_ROOT / scenario_name / "eur_discovery_gwas.SYNTH_PHENO.glm.linear"
    if not gwas_glm_linear.exists():
        raise SystemExit(f"{gwas_glm_linear} not found — run src/trackA_synthetic/stage3_eur_gwas/run_gwas.py first.")

    out_dir = PRS_MODELS_ROOT / scenario_name / "prsice2"
    out_dir.mkdir(parents=True, exist_ok=True)

    sumstats_clean = out_dir / "sumstats_clean.tsv"
    clean_sumstats(gwas_glm_linear, sumstats_clean)

    out_prefix = out_dir / "prsice2"
    run([
        "Rscript", PRSICE_R,
        "--prsice", PRSICE_BIN,
        "--base", sumstats_clean,
        "--chr", "CHROM", "--bp", "POS", "--snp", "ID", "--A1", "A1", "--stat", "BETA", "--pvalue", "P", "--beta",
        "--target", BED_PREFIX,
        "--keep", EUR_DISCOVERY_KEEP_FILE,
        "--binary-target", "F",
        "--no-regress",
        "--fastscore", "--bar-levels", ",".join(BAR_LEVELS),
        "--print-snp",
        "--clump-r2", CLUMP_R2, "--clump-kb", CLUMP_KB,
        "--out", out_prefix,
    ], cwd=out_dir)

    sumstats = load_sumstats(sumstats_clean)
    snp_p = parse_clumped_snps(out_prefix.with_suffix(".snp"))
    threshold_counts = parse_threshold_counts(out_prefix.with_suffix(".prsice"))

    scoring_files = []
    for threshold_str, n_snp in threshold_counts:
        if n_snp == 0:
            logging.info("%s: threshold %s has 0 SNPs after clumping — no scoring file written", scenario_name, threshold_str)
            continue
        threshold = float(threshold_str)
        selected = [snp for snp, p in snp_p.items() if p <= threshold]
        dest = out_dir / f"threshold_{threshold_str}.scoring.tsv"
        write_scoring_file(selected, sumstats, dest)
        scoring_files.append((threshold_str, len(selected), dest))
        logging.info("%s: threshold %s -> %d SNPs -> %s", scenario_name, threshold_str, len(selected), dest)

    return scoring_files


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenario", choices=["1", "2", "all"], default="all")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    scenario_ids = ["1", "2"] if args.scenario == "all" else [args.scenario]
    for sid in scenario_ids:
        scoring_files = run_prsice2_for_scenario(SCENARIOS[sid])
        logging.info("PRSice-2 complete for scenario %s: %d usable thresholds", sid, len(scoring_files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
