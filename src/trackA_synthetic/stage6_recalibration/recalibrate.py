#!/usr/bin/env python3
"""
recalibrate.py — Stage 6: apply recalibration (BUILD_PLAN.md §6 Stage 6).

BUILD_PLAN.md's first choice, `pgscatalog-ancestry-adjust`, requires PCA projections from
`fraposa_pgsc` and PGS-Catalog-format aggregated scores — outputs of the full `pgsc_calc`
Nextflow pipeline that this project's Track A (synthetic, custom-built scoring files) never
produces. BUILD_PLAN.md explicitly anticipates this ("implement a standard empirical
ancestry-matched recentering/rescaling approach if the tool's assumptions don't cleanly fit our
synthetic setup"), so that's what this does: for each (scenario, method, [threshold]) model,
recenter and rescale each non-EUR ancestry's PRS distribution to match the EUR held-out
reference distribution's mean/SD (the standard empirical recalibration approach).

Important, and the actual point of this stage per BUILD_PLAN.md §9 item 2: a linear
recentering/rescaling transform cannot change Pearson R^2 — R^2 depends only on correlation,
which is invariant to affine transforms of either variable. So R^2 before/after recalibration
will come out numerically identical here, *by mathematical necessity*, not because recalibration
"didn't work". What recalibration actually fixes is calibration (the score's mean/scale being
appropriate per ancestry) — captured here by each ancestry's raw vs. recalibrated PRS mean/SD
relative to the EUR reference. Reporting both side by side makes BUILD_PLAN's warned-about
distinction (discriminative accuracy vs. calibration are different things) concrete rather than
asserted.

Usage:
    python3 recalibrate.py                  # both scenarios
"""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import statistics
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
SHARED_DIR = REPO_ROOT / "data" / "processed" / "simulated_phenotypes" / "_shared"
BED_PREFIX = SHARED_DIR / "1000g_qc_biallelic"
GWAS_SHARED_DIR = REPO_ROOT / "data" / "processed" / "gwas" / "_shared"
PRS_MODELS_ROOT = REPO_ROOT / "data" / "processed" / "prs_models"
SIM_PHENO_ROOT = REPO_ROOT / "data" / "processed" / "simulated_phenotypes"
EVAL_DIR = REPO_ROOT / "data" / "processed" / "evaluation"
OUT_DIR = REPO_ROOT / "data" / "processed" / "recalibration"

SCENARIOS = {"1": "scenario1_equal_effects", "2": "scenario2_ancestry_varying_effects"}

EVAL_GROUPS = {
    "EUR_holdout": GWAS_SHARED_DIR / "keep_EUR_holdout.txt",
    "AFR": SHARED_DIR / "keep_AFR.txt",
    "AMR": SHARED_DIR / "keep_AMR.txt",
    "EAS": SHARED_DIR / "keep_EAS.txt",
    "SAS": SHARED_DIR / "keep_SAS.txt",
}
REFERENCE_ANCESTRY = "EUR_holdout"


def require_plink2() -> str:
    exe = shutil.which("plink2")
    if exe is None:
        raise SystemExit("plink2 not found on PATH. Run scripts/verify_tools.py --install first.")
    return exe


def run(cmd: list) -> None:
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


def r2(x: list, y: list) -> float:
    n = len(x)
    mean_x, mean_y = sum(x) / n, sum(y) / n
    cov = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    var_x = sum((a - mean_x) ** 2 for a in x)
    var_y = sum((b - mean_y) ** 2 for b in y)
    if var_x == 0 or var_y == 0:
        return 0.0
    r = cov / (var_x * var_y) ** 0.5
    return r ** 2


def discover_models(scenario_name: str) -> list:
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

        for method, threshold, scoring_file in models:
            # Reference distribution (EUR held-out) computed once per model.
            ref_prefix = work_dir / f"{scenario_name}_{method}_{threshold}_{REFERENCE_ANCESTRY}"
            ref_scores = compute_prs(plink2, scoring_file, EVAL_GROUPS[REFERENCE_ANCESTRY], ref_prefix)
            ref_vals = list(ref_scores.values())
            ref_mean, ref_sd = statistics.mean(ref_vals), statistics.stdev(ref_vals)

            for ancestry, keep_file in EVAL_GROUPS.items():
                out_prefix = work_dir / f"{scenario_name}_{method}_{threshold}_{ancestry}"
                scores = compute_prs(plink2, scoring_file, keep_file, out_prefix)
                common_ids = [i for i in scores if i in true_pheno]
                raw_vals = [scores[i] for i in common_ids]
                pheno_vals = [true_pheno[i] for i in common_ids]
                n = len(common_ids)

                raw_mean, raw_sd = statistics.mean(raw_vals), statistics.stdev(raw_vals)
                # Empirical recalibration: z-score within the ancestry's own distribution, then
                # rescale to the EUR reference's mean/SD (pgscatalog-ancestry-adjust's
                # "empirical" method, applied without needing its PCA/aggregated-score inputs).
                if raw_sd == 0:
                    recal_vals = [ref_mean for _ in raw_vals]
                else:
                    recal_vals = [(v - raw_mean) / raw_sd * ref_sd + ref_mean for v in raw_vals]
                recal_mean, recal_sd = statistics.mean(recal_vals), statistics.stdev(recal_vals)

                r2_raw = r2(raw_vals, pheno_vals)
                r2_recal = r2(recal_vals, pheno_vals)

                rows.append({
                    "scenario": scenario_name, "method": method, "threshold": threshold, "ancestry": ancestry, "n": n,
                    "raw_mean": round(raw_mean, 4), "raw_sd": round(raw_sd, 4),
                    "recal_mean": round(recal_mean, 4), "recal_sd": round(recal_sd, 4),
                    "ref_mean": round(ref_mean, 4), "ref_sd": round(ref_sd, 4),
                    "r2_raw": round(r2_raw, 4), "r2_recalibrated": round(r2_recal, 4),
                })
                logging.info(
                    "%s | %s %s | %s: mean %.3f->%.3f (ref %.3f), sd %.3f->%.3f (ref %.3f), R2 %.4f->%.4f",
                    scenario_name, method, threshold, ancestry, raw_mean, recal_mean, ref_mean,
                    raw_sd, recal_sd, ref_sd, r2_raw, r2_recal,
                )

    dest = OUT_DIR / "results.tsv"
    fieldnames = ["scenario", "method", "threshold", "ancestry", "n", "raw_mean", "raw_sd",
                  "recal_mean", "recal_sd", "ref_mean", "ref_sd", "r2_raw", "r2_recalibrated"]
    with open(dest, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Stage 6 recalibration complete: %d rows -> %s", len(rows), dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
