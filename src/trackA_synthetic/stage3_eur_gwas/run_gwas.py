#!/usr/bin/env python3
"""
run_gwas.py — Stage 3: EUR-only discovery GWAS (BUILD_PLAN.md §6 Stage 3).

Runs a simple, unadjusted per-SNP linear regression GWAS on Stage 2's simulated phenotype,
restricted to EUR samples only — deliberately reproducing the real-world bias that most
published GWAS are European-ancestry-derived. Run separately per scenario, since Scenario 1
and Scenario 2 have different simulated phenotypes (BUILD_PLAN.md Definition of Done item 3:
"trained on a EUR-only discovery GWAS ... Repeat for both Scenario 1 and Scenario 2").

Genotypes are the same shared biallelic bed/bim/fam Stage 2 used to simulate the phenotype
(data/processed/simulated_phenotypes/_shared/), so the causal variants are guaranteed to be
present in the GWAS output for Stage 4's PRS construction to pick up. No PCs/covariates are
used — BUILD_PLAN.md explicitly calls for "a simple linear regression per-SNP" as the
discovery GWAS, and the simulated phenotype has no population-structure confounding to adjust
for (it's pure per-ancestry genetic effect + residual noise, not stratified across ancestries
within the EUR-only regression).

Before running the GWAS, EUR samples are split into a discovery subset (used for the GWAS) and
a held-out subset (config: discovery_gwas.holdout_fraction, default 0.2) that the GWAS never
sees. Stage 5 needs this held-out EUR set for a genuine same-ancestry evaluation — BUILD_PLAN.md
§6 Stage 5 explicitly calls for "a held-out EUR subset as the 'easy,' same-ancestry baseline
case", and §9 item 1's sanity check requires a same-ancestry train/test split. Without holding
samples out, "EUR" evaluation in Stage 5 would be in-sample against the GWAS training data.

Usage:
    python3 run_gwas.py                  # both scenarios
    python3 run_gwas.py --scenario 1     # Scenario 1 only
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
CONFIG_PATH = REPO_ROOT / "configs" / "simulation_parameters.yaml"
SHARED_DIR = REPO_ROOT / "data" / "processed" / "simulated_phenotypes" / "_shared"
BED_PREFIX = SHARED_DIR / "1000g_qc_biallelic"
EUR_KEEP_FILE = SHARED_DIR / "keep_EUR.txt"
GWAS_ROOT = REPO_ROOT / "data" / "processed" / "gwas"
GWAS_SHARED_DIR = GWAS_ROOT / "_shared"
EUR_DISCOVERY_KEEP_FILE = GWAS_SHARED_DIR / "keep_EUR_discovery.txt"
EUR_HOLDOUT_KEEP_FILE = GWAS_SHARED_DIR / "keep_EUR_holdout.txt"

SCENARIOS = {
    "1": "scenario1_equal_effects",
    "2": "scenario2_ancestry_varying_effects",
}


def split_eur_discovery_holdout(config: dict) -> None:
    """Partition EUR samples into a GWAS-discovery subset and a held-out subset the GWAS never
    sees, so Stage 5 has a genuine same-ancestry evaluation set. Idempotent; cached so the same
    split is reused by both scenarios."""
    if EUR_DISCOVERY_KEEP_FILE.exists() and EUR_HOLDOUT_KEEP_FILE.exists():
        return
    GWAS_SHARED_DIR.mkdir(parents=True, exist_ok=True)

    eur_ids = [line.split("\t")[1] for line in EUR_KEEP_FILE.read_text().splitlines() if line.strip()]
    gwas_cfg = config["discovery_gwas"]
    rng = np.random.default_rng(gwas_cfg["holdout_seed"])
    n_holdout = round(len(eur_ids) * gwas_cfg["holdout_fraction"])
    holdout = set(rng.choice(eur_ids, size=n_holdout, replace=False).tolist())
    discovery = [i for i in eur_ids if i not in holdout]

    EUR_DISCOVERY_KEEP_FILE.write_text("".join(f"0\t{i}\n" for i in discovery))
    EUR_HOLDOUT_KEEP_FILE.write_text("".join(f"0\t{i}\n" for i in sorted(holdout)))
    logging.info(
        "split %d EUR samples into %d discovery (GWAS) / %d held-out (Stage 5 eval)",
        len(eur_ids), len(discovery), len(holdout),
    )


def require_plink2() -> str:
    exe = shutil.which("plink2")
    if exe is None:
        raise SystemExit("plink2 not found on PATH. Run scripts/verify_tools.py --install first.")
    return exe


def run(cmd: list) -> None:
    logging.info("+ %s", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True)


PHENO_COL = "SYNTH_PHENO"  # not "PHENO1": plink2's bed/fam already has an implicit PHENO1
# column (the fam file's 6th column, "-9"/missing for all samples here), and a --pheno file
# reusing that name collides with it ("Duplicate phenotype/covariate ID 'PHENO1'").


def write_eur_discovery_pheno_file(phenotypes_tsv: Path, discovery_ids: set, dest: Path) -> int:
    """Extract EUR-discovery-only rows from Stage 2's combined phenotypes.tsv into a plink2
    --pheno file (held-out EUR samples are excluded so the GWAS never sees them)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(phenotypes_tsv) as fin, open(dest, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        fout.write(f"#FID\tIID\t{PHENO_COL}\n")
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if fields[idx["super_pop"]] != "EUR" or fields[idx["IID"]] not in discovery_ids:
                continue
            fout.write(f"{fields[idx['FID']]}\t{fields[idx['IID']]}\t{fields[idx['phenotype']]}\n")
            n += 1
    return n


def run_gwas_for_scenario(plink2: str, scenario_name: str) -> Path:
    phenotypes_tsv = REPO_ROOT / "data" / "processed" / "simulated_phenotypes" / scenario_name / "phenotypes.tsv"
    if not phenotypes_tsv.exists():
        raise SystemExit(
            f"{phenotypes_tsv} not found — run "
            f"src/trackA_synthetic/stage2_gcta_simulation/{scenario_name.split('_', 1)[0]}_*.py first."
        )

    discovery_ids = {
        line.split("\t")[1] for line in EUR_DISCOVERY_KEEP_FILE.read_text().splitlines() if line.strip()
    }

    out_dir = GWAS_ROOT / scenario_name
    pheno_file = out_dir / "eur_discovery_pheno.tsv"
    n_eur = write_eur_discovery_pheno_file(phenotypes_tsv, discovery_ids, pheno_file)
    logging.info("%s: %d EUR discovery samples with simulated phenotype", scenario_name, n_eur)

    out_prefix = out_dir / "eur_discovery_gwas"
    run([
        plink2, "--bfile", BED_PREFIX,
        "--keep", EUR_DISCOVERY_KEEP_FILE,
        "--pheno", pheno_file,
        # The bed/fam's implicit 6th-column phenotype is "-9" (missing) for every sample and
        # would otherwise also get loaded and run through --glm as a second, constant
        # phenotype; --pheno-name restricts --glm to just our simulated phenotype column.
        "--pheno-name", PHENO_COL,
        "--glm", "allow-no-covars",
        "--out", out_prefix,
    ])
    result = out_prefix.with_name(out_prefix.name + f".{PHENO_COL}.glm.linear")
    if not result.exists():
        raise SystemExit(f"expected GWAS output not found: {result}")
    return result


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenario", choices=["1", "2", "all"], default="all")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    plink2 = require_plink2()
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    split_eur_discovery_holdout(config)

    scenario_ids = ["1", "2"] if args.scenario == "all" else [args.scenario]
    for sid in scenario_ids:
        result = run_gwas_for_scenario(plink2, SCENARIOS[sid])
        logging.info("Stage 3 GWAS complete for scenario %s: %s", sid, result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
