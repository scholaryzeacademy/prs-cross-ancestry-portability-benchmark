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

REPO_ROOT = Path(__file__).resolve().parents[3]
SHARED_DIR = REPO_ROOT / "data" / "processed" / "simulated_phenotypes" / "_shared"
BED_PREFIX = SHARED_DIR / "1000g_qc_biallelic"
EUR_KEEP_FILE = SHARED_DIR / "keep_EUR.txt"
GWAS_ROOT = REPO_ROOT / "data" / "processed" / "gwas"

SCENARIOS = {
    "1": "scenario1_equal_effects",
    "2": "scenario2_ancestry_varying_effects",
}


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


def write_eur_pheno_file(phenotypes_tsv: Path, dest: Path) -> int:
    """Extract EUR-only rows from Stage 2's combined phenotypes.tsv into a plink2 --pheno file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(phenotypes_tsv) as fin, open(dest, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        fout.write(f"#FID\tIID\t{PHENO_COL}\n")
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if fields[idx["super_pop"]] != "EUR":
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

    out_dir = GWAS_ROOT / scenario_name
    pheno_file = out_dir / "eur_pheno.tsv"
    n_eur = write_eur_pheno_file(phenotypes_tsv, pheno_file)
    logging.info("%s: %d EUR samples with simulated phenotype", scenario_name, n_eur)

    out_prefix = out_dir / "eur_discovery_gwas"
    run([
        plink2, "--bfile", BED_PREFIX,
        "--keep", EUR_KEEP_FILE,
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
    scenario_ids = ["1", "2"] if args.scenario == "all" else [args.scenario]
    for sid in scenario_ids:
        result = run_gwas_for_scenario(plink2, SCENARIOS[sid])
        logging.info("Stage 3 GWAS complete for scenario %s: %s", sid, result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
