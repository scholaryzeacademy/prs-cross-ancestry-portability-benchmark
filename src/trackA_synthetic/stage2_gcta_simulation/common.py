"""
common.py — shared plumbing for Stage 2 GCTA phenotype simulation (BUILD_PLAN.md §6 Stage 2).

Both scenario1_equal_effects.py and scenario2_ancestry_varying_effects.py call into this
module so the causal-variant set, bed/bim/fam conversion, and per-ancestry GCTA invocation are
identical between scenarios — the only thing that differs is which effect sizes go into the
causal-loci file handed to GCTA.

Design notes
------------
* GCTA requires PLINK1 bed/bim/fam (no multiallelic support), so Stage 1's QC'd pgen is
  converted once into a shared bed/bim/fam restricted to biallelic SNPs.
* Stage 1's pgen has no #FID column (plink2 IID-only mode), so plink2 assigns FID=0 when
  writing bed/fam. GCTA's --keep therefore needs "0<tab>IID" two-column files, not the
  single-column IID files Stage 1 wrote for pgen's --keep.
* GCTA computes genetic value from the causal-loci effect sizes and each run's own sample's
  genotypes/allele frequencies, then picks residual noise variance so heritability matches
  --simu-hsq *within that run's sample*. Running GCTA once per super-population (via --keep on
  one shared bed/bim/fam, so allele coding/effect-allele stays fixed across runs) gives every
  ancestry the same nominal heritability while letting real per-ancestry allele-frequency/LD
  differences (and, in Scenario 2, different effect sizes) drive the portability gap.
* GCTA v1.95.3 has no --seed option for --simu-qt's residual noise — the genetic component
  (causal variants, effect sizes, per-ancestry perturbations) is fully reproducible via our own
  seeded RNG, but the residual draw is not bit-for-bit reproducible across GCTA re-runs. This is
  a genuine tooling limitation, documented rather than worked around.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

import numpy as np
import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
CONFIG_PATH = REPO_ROOT / "configs" / "simulation_parameters.yaml"
QC_PREFIX = REPO_ROOT / "data" / "processed" / "1000g_qc" / "1000g_qc"
SIM_ROOT = REPO_ROOT / "data" / "processed" / "simulated_phenotypes"
SHARED_DIR = SIM_ROOT / "_shared"

SUPER_POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def load_config(path: Path = CONFIG_PATH) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def run(cmd: list) -> None:
    logging.info("+ %s", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True)


def require_plink2() -> str:
    exe = shutil.which("plink2")
    if exe is None:
        raise SystemExit("plink2 not found on PATH. Run scripts/verify_tools.py --install first.")
    return exe


def require_gcta() -> str:
    exe = shutil.which("gcta64") or shutil.which("gcta")
    if exe:
        return exe
    tools_dir = REPO_ROOT / "tools"
    if tools_dir.exists():
        for pattern in ("gcta64", "gcta"):
            for p in sorted(tools_dir.rglob(pattern)):
                if p.is_file():
                    return str(p)
    raise SystemExit("gcta not found on PATH or under ./tools. Run scripts/verify_tools.py --install first.")


def prepare_bed(plink2: str) -> Path:
    """Convert Stage 1's QC'd pgen into a shared, biallelic-SNP-only bed/bim/fam. Idempotent."""
    bed_prefix = SHARED_DIR / "1000g_qc_biallelic"
    if bed_prefix.with_suffix(".bed").exists():
        return bed_prefix
    SHARED_DIR.mkdir(parents=True, exist_ok=True)
    if not QC_PREFIX.with_suffix(".pgen").exists():
        raise SystemExit(
            f"{QC_PREFIX}.pgen not found — run src/shared_stage1_1000g_download_qc/qc.py first."
        )
    run([
        plink2, "--pfile", QC_PREFIX,
        "--max-alleles", "2", "--snps-only", "just-acgt",
        "--make-bed", "--out", bed_prefix,
    ])
    return bed_prefix


def prepare_gcta_keep_files() -> dict:
    """GCTA needs FID+IID keep files ("0<tab>IID", since plink2 assigned FID=0 on bed
    conversion); Stage 1 wrote single-column IID files for pgen's --keep. Idempotent."""
    SHARED_DIR.mkdir(parents=True, exist_ok=True)
    keep_files = {}
    for sp in SUPER_POPULATIONS:
        src = REPO_ROOT / "data" / "processed" / "1000g_qc" / f"samples_{sp}.txt"
        if not src.exists():
            raise SystemExit(f"{src} not found — run src/shared_stage1_1000g_download_qc/qc.py first.")
        dest = SHARED_DIR / f"keep_{sp}.txt"
        with open(src) as fin, open(dest, "w") as fout:
            for line in fin:
                iid = line.strip()
                if iid:
                    fout.write(f"0\t{iid}\n")
        keep_files[sp] = dest
    return keep_files


def select_causal_variants(plink2: str, bed_prefix: Path, config: dict) -> list[str]:
    """Pick n_causal variants with MAF in [maf_min, maf_max] from the shared bed/bim/fam,
    using the fixed random_seed in configs/simulation_parameters.yaml. Cached so both
    scenarios simulate on the exact same causal-variant set."""
    cache = SHARED_DIR / "causal_variants.txt"
    if cache.exists():
        return cache.read_text().split()

    cv_cfg = config["causal_variants"]
    freq_prefix = SHARED_DIR / "freq"
    run([plink2, "--bfile", bed_prefix, "--freq", "--out", freq_prefix])

    eligible = []
    with open(freq_prefix.with_suffix(".afreq")) as f:
        header = f.readline().lstrip("#").split()
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            fields = line.split()
            snp_id = fields[idx["ID"]]
            alt_freq = float(fields[idx["ALT_FREQS"]])
            maf = min(alt_freq, 1 - alt_freq)
            if cv_cfg["maf_min"] <= maf <= cv_cfg["maf_max"]:
                eligible.append(snp_id)

    rng = np.random.default_rng(cv_cfg["random_seed"])
    n_causal = cv_cfg["n_causal"]
    if len(eligible) < n_causal:
        raise SystemExit(f"only {len(eligible)} variants eligible (MAF filter), need {n_causal}")
    chosen = sorted(rng.choice(eligible, size=n_causal, replace=False).tolist())

    cache.write_text("\n".join(chosen) + "\n")
    logging.info("selected %d causal variants (%d eligible by MAF filter)", n_causal, len(eligible))
    return chosen


def generate_base_effects(causal_variants: list[str], seed: int) -> dict:
    rng = np.random.default_rng(seed)
    effects = rng.normal(loc=0.0, scale=1.0, size=len(causal_variants))
    return dict(zip(causal_variants, effects))


def write_causal_loci_file(effects: dict, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w") as f:
        for snp_id, effect in effects.items():
            f.write(f"{snp_id}\t{effect}\n")


def run_gcta_simu(gcta: str, bed_prefix: Path, keep_file: Path, causal_loci_file: Path, hsq: float, out_prefix: Path) -> Path:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    run([
        gcta, "--bfile", bed_prefix, "--keep", keep_file,
        "--simu-qt", "--simu-causal-loci", causal_loci_file,
        "--simu-hsq", hsq, "--simu-rep", "1",
        "--out", out_prefix,
    ])
    return out_prefix.with_suffix(".phen")


def concatenate_phenotypes(phen_by_pop: dict, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w") as fout:
        fout.write("FID\tIID\tsuper_pop\tphenotype\n")
        for sp, phen_path in phen_by_pop.items():
            with open(phen_path) as fin:
                for line in fin:
                    fid, iid, pheno = line.split()
                    fout.write(f"{fid}\t{iid}\t{sp}\t{pheno}\n")
    logging.info("wrote combined phenotype file: %s", dest)
