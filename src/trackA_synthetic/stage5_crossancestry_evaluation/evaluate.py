#!/usr/bin/env python3
"""
evaluate.py — Stage 5: cross-ancestry predictive accuracy evaluation (BUILD_PLAN.md §6 Stage 5).

Applies every Stage 4 PRS model (PRSice-2 at each p-value threshold, LDpred2-auto, PRS-CSx) to
held-out samples from all five super-populations — EUR uses the held-out subset Stage 3 set
aside (never seen by the discovery GWAS or any Stage 4 construction step), AFR/AMR/EAS/SAS use
their full sample sets (never seen at all, the genuine cross-ancestry portability test).

For each (scenario, method, [threshold], ancestry) combination, computes R^2 between the PRS
and the *true* simulated phenotype (available because this is synthetic ground truth — real PRS
evaluations only have a proxy phenotype, not this) with a 95% CI via the Fisher z-transform of
the Pearson correlation. Every combination is reported in one full table — BUILD_PLAN.md §9 item
3 explicitly warns against collapsing this into a single "best method" statistic.

Usage:
    python3 evaluate.py                  # both scenarios, all methods, all ancestries
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
SHARED_DIR = REPO_ROOT / "data" / "processed" / "simulated_phenotypes" / "_shared"
BED_PREFIX = SHARED_DIR / "1000g_qc_biallelic"
GWAS_SHARED_DIR = REPO_ROOT / "data" / "processed" / "gwas" / "_shared"
PRS_MODELS_ROOT = REPO_ROOT / "data" / "processed" / "prs_models"
SIM_PHENO_ROOT = REPO_ROOT / "data" / "processed" / "simulated_phenotypes"
OUT_DIR = REPO_ROOT / "data" / "processed" / "evaluation"

SCENARIOS = {
    "1": "scenario1_equal_effects",
    "2": "scenario2_ancestry_varying_effects",
}

# Ancestry -> keep file. EUR uses the held-out subset (never in the discovery GWAS or Stage 4
# construction); AFR/AMR/EAS/SAS use their full sample sets (never used anywhere upstream).
EVAL_GROUPS = {
    "EUR_holdout": GWAS_SHARED_DIR / "keep_EUR_holdout.txt",
    "AFR": SHARED_DIR / "keep_AFR.txt",
    "AMR": SHARED_DIR / "keep_AMR.txt",
    "EAS": SHARED_DIR / "keep_EAS.txt",
    "SAS": SHARED_DIR / "keep_SAS.txt",
}


def require_plink2() -> str:
    exe = shutil.which("plink2")
    if exe is None:
        raise SystemExit("plink2 not found on PATH. Run scripts/verify_tools.py --install first.")
    return exe


def run(cmd: list) -> None:
    logging.info("+ %s", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True, capture_output=True)


def load_true_phenotype(scenario_name: str) -> dict:
    path = SIM_PHENO_ROOT / scenario_name / "phenotypes.tsv"
    pheno = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pheno[row["IID"]] = float(row["phenotype"])
    return pheno


def compute_prs(plink2: str, scoring_file: Path, keep_file: Path, out_prefix: Path) -> dict:
    run([
        plink2, "--bfile", BED_PREFIX, "--keep", keep_file,
        "--score", scoring_file, "1", "2", "3", "header", "cols=+scoresums",
        "--out", out_prefix,
    ])
    scores = {}
    with open(f"{out_prefix}.sscore") as f:
        header = f.readline().lstrip("#").split()
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            fields = line.split()
            scores[fields[idx["IID"]]] = float(fields[idx["SCORE1_SUM"]])
    return scores


def r2_with_ci(prs: list, pheno: list) -> tuple:
    n = len(prs)
    mean_x, mean_y = sum(prs) / n, sum(pheno) / n
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(prs, pheno))
    var_x = sum((x - mean_x) ** 2 for x in prs)
    var_y = sum((y - mean_y) ** 2 for y in pheno)
    if var_x == 0 or var_y == 0:
        return 0.0, 0.0, 0.0, n
    r = cov / math.sqrt(var_x * var_y)
    r = max(min(r, 0.999999), -0.999999)
    if n <= 3:
        return r ** 2, float("nan"), float("nan"), n
    z = math.atanh(r)
    se = 1 / math.sqrt(n - 3)
    z_lo, z_hi = z - 1.96 * se, z + 1.96 * se
    r_lo, r_hi = math.tanh(z_lo), math.tanh(z_hi)
    # r2's CI bounds follow whichever of r_lo/r_hi is closer to zero vs. farther, handling sign
    r2_candidates = sorted([r_lo ** 2, r_hi ** 2])
    return r ** 2, r2_candidates[0], r2_candidates[1], n


def discover_models(scenario_name: str) -> list:
    """Return (method, threshold_label, scoring_file) for every Stage 4 output found."""
    models = []
    prsice_dir = PRS_MODELS_ROOT / scenario_name / "prsice2"
    for f in sorted(prsice_dir.glob("threshold_*.scoring.tsv")):
        threshold = f.stem.replace("threshold_", "").replace(".scoring", "")
        models.append(("prsice2", threshold, f))

    ldpred2_file = PRS_MODELS_ROOT / scenario_name / "ldpred2" / "ldpred2_auto.scoring.tsv"
    if ldpred2_file.exists():
        models.append(("ldpred2", "-", ldpred2_file))

    prscsx_file = PRS_MODELS_ROOT / scenario_name / "prscsx" / "prscsx.scoring.tsv"
    if prscsx_file.exists():
        models.append(("prscsx", "-", prscsx_file))

    return models


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenario", choices=["1", "2", "all"], default="all")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    plink2 = require_plink2()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    work_dir = OUT_DIR / "_work"
    work_dir.mkdir(exist_ok=True)

    scenario_ids = ["1", "2"] if args.scenario == "all" else [args.scenario]
    rows = []
    for sid in scenario_ids:
        scenario_name = SCENARIOS[sid]
        true_pheno = load_true_phenotype(scenario_name)
        models = discover_models(scenario_name)
        logging.info("%s: %d models found (%s)", scenario_name, len(models),
                     ", ".join(sorted({m[0] for m in models})))

        for method, threshold, scoring_file in models:
            for ancestry, keep_file in EVAL_GROUPS.items():
                out_prefix = work_dir / f"{scenario_name}_{method}_{threshold}_{ancestry}"
                scores = compute_prs(plink2, scoring_file, keep_file, out_prefix)
                common_ids = [i for i in scores if i in true_pheno]
                prs_vals = [scores[i] for i in common_ids]
                pheno_vals = [true_pheno[i] for i in common_ids]
                r2, r2_lo, r2_hi, n = r2_with_ci(prs_vals, pheno_vals)
                rows.append({
                    "scenario": scenario_name, "method": method, "threshold": threshold,
                    "ancestry": ancestry, "n": n, "r2": round(r2, 4),
                    "r2_ci_low": round(r2_lo, 4) if r2_lo == r2_lo else "NA",
                    "r2_ci_high": round(r2_hi, 4) if r2_hi == r2_hi else "NA",
                })
                logging.info("%s | %s %s | %s: n=%d R2=%.4f [%.4f, %.4f]",
                             scenario_name, method, threshold, ancestry, n, r2, r2_lo, r2_hi)

    dest = OUT_DIR / "results.tsv"
    with open(dest, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["scenario", "method", "threshold", "ancestry", "n", "r2", "r2_ci_low", "r2_ci_high"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Stage 5 evaluation complete: %d rows -> %s", len(rows), dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
