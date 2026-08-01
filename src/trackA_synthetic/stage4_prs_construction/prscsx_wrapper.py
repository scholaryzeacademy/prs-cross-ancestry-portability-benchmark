#!/usr/bin/env python3
"""
prscsx_wrapper.py — Stage 4, PRS-CSx (BUILD_PLAN.md §6 Stage 4).

Runs real PRS-CSx (Ge et al. 2019 / Ruan et al. 2022) on Stage 3's EUR discovery-only GWAS,
using the official 1000G EUR LD reference panel. BUILD_PLAN.md feeds all three Stage 4 methods
the same EUR-only sumstats, so this runs PRS-CSx in its single-discovery-population mode
(--pop EUR) — which is exactly PRS-CS, the non-multi-ancestry predecessor method; PRS-CSx's
"x" (coupling genetic effects *across* ancestries) only activates when multiple discovery
GWASes are supplied, which this project's design deliberately doesn't do (BUILD_PLAN.md §0: the
whole point is testing how a EUR-only score ports elsewhere, not building a multi-ancestry one).

Requires the official per-population LD reference panel (multi-GB per population; see
download_ld_reference() below) — this is real published tooling, not reproduced here.

SNP ID compatibility: PRS-CSx's reference panel is restricted to ~1.1M HapMap3 SNPs, addressed
by rsID. Our own variant IDs (from Stage 1's `--set-all-var-ids @:#:$r:$a`) are
chr:pos:ref:alt, not rsIDs, so this wrapper joins our variants against the PRS-CSx-provided
snpinfo_mult_1kg_hm3 file (matched by chromosome+position) to get a rsID <-> our-ID mapping,
runs PRS-CSx on the rsID-labeled data, then translates the output back to our own IDs so Stage
5 can apply the resulting scoring file via `plink2 --score` against our own bed/bim/fam exactly
like the PRSice-2/LDpred2 outputs.

Usage:
    python3 prscsx_wrapper.py                  # both scenarios
    python3 prscsx_wrapper.py --scenario 1
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
BIM_PATH = SHARED_DIR / "1000g_qc_biallelic.bim"
GWAS_ROOT = REPO_ROOT / "data" / "processed" / "gwas"
EUR_DISCOVERY_KEEP_FILE = GWAS_ROOT / "_shared" / "keep_EUR_discovery.txt"
PRS_MODELS_ROOT = REPO_ROOT / "data" / "processed" / "prs_models"

PRSCSX_DIR = REPO_ROOT / "tools" / "prscsx" / "PRScsx-master"
# PRS-CSx expects snpinfo_mult_1kg_hm3 directly inside --ref_dir, alongside the per-population
# ldblk_1kg_<pop>/ subdirectories (see PRScsx.py's main(): both are read as ref_dir + '/...').
LD_REF_DIR = REPO_ROOT / "tools" / "prscsx" / "ld_reference" / "1kg"
SNPINFO_PATH = LD_REF_DIR / "snpinfo_mult_1kg_hm3"

CHROMS = ["21", "22"]  # matches this project's smoke-test genotype panel (see BUILD_PLAN.md §5 CI row)
POP = "EUR"

PRSCSX_LD_URLS = {
    "eur": "https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=1",
    "afr": "https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=0",
    "amr": "https://www.dropbox.com/s/uv5ydr4uv528lca/ldblk_1kg_amr.tar.gz?dl=0",
    "eas": "https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=0",
    "sas": "https://www.dropbox.com/s/hsm0qwgyixswdcv/ldblk_1kg_sas.tar.gz?dl=0",
}
SNPINFO_URL = "https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3?dl=1"


def download_ld_reference(pop: str = "eur") -> None:
    """One-time fetch of PRS-CSx's official LD reference panel + HapMap3 SNP info file.
    Several GB per population; not run automatically by run_prscsx_for_scenario()."""
    LD_REF_DIR.mkdir(parents=True, exist_ok=True)
    if not SNPINFO_PATH.exists():
        subprocess.run(["curl", "-L", "-C", "-", "--retry", "999", "--retry-delay", "10",
                         "--retry-all-errors", "-o", str(SNPINFO_PATH), SNPINFO_URL], check=True)
    tarball = LD_REF_DIR / f"ldblk_1kg_{pop}.tar.gz"
    extracted = LD_REF_DIR / f"ldblk_1kg_{pop}"
    if extracted.exists():
        return
    url = PRSCSX_LD_URLS[pop].replace("dl=0", "dl=1")
    subprocess.run(["curl", "-L", "-C", "-", "--retry", "999", "--retry-delay", "10",
                     "--retry-all-errors", "-o", str(tarball), url], check=True)
    subprocess.run(["tar", "xzf", str(tarball)], cwd=LD_REF_DIR, check=True)
    tarball.unlink()


def build_rsid_map() -> dict:
    """{(chr, pos): rsid} restricted to CHROMS, from PRS-CSx's HapMap3 SNP info file."""
    rsid_map = {}
    with open(SNPINFO_PATH) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["CHR"] in CHROMS:
                rsid_map[(row["CHR"], row["BP"])] = row["SNP"]
    return rsid_map


def write_target_bim(rsid_map: dict, dest_prefix: Path) -> dict:
    """Write a rsID-labeled copy of our bim (PRS-CSx only ever reads the .bim text, never the
    matching .bed/.fam) restricted to variants present in the HapMap3 reference. Returns
    {rsid: our_variant_id} for translating PRS-CSx's output back to our own ID convention."""
    dest_prefix.parent.mkdir(parents=True, exist_ok=True)
    rsid_to_ours = {}
    with open(BIM_PATH) as fin, open(dest_prefix.with_suffix(".bim"), "w") as fout:
        for line in fin:
            chrom, our_id, gd, pos, a1, a2 = line.split()
            rsid = rsid_map.get((chrom, pos))
            if rsid is None:
                continue
            fout.write(f"{chrom}\t{rsid}\t{gd}\t{pos}\t{a1}\t{a2}\n")
            rsid_to_ours[rsid] = our_id
    return rsid_to_ours


def write_sumstats(gwas_glm_linear: Path, rsid_map: dict, dest: Path) -> int:
    """PRS-CSx sumstats format: SNP A1 A2 BETA SE (see parse_genet.py's positional parsing)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(gwas_glm_linear) as fin, open(dest, "w") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fout.write("SNP\tA1\tA2\tBETA\tSE\n")
        for row in reader:
            rsid = rsid_map.get((row["#CHROM"], row["POS"]))
            if rsid is None or row["BETA"] == "NA":
                continue
            fout.write(f"{rsid}\t{row['A1']}\t{row['OMITTED']}\t{row['BETA']}\t{row['SE']}\n")
            n += 1
    return n


def count_samples(keep_file: Path) -> int:
    return sum(1 for line in open(keep_file) if line.strip())


def run(cmd: list, cwd: Path = None) -> None:
    logging.info("+ %s", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True, cwd=cwd)


def run_prscsx_for_scenario(scenario_name: str, rsid_map: dict, rsid_to_ours: dict, target_bim_prefix: Path) -> Path:
    gwas_glm_linear = GWAS_ROOT / scenario_name / "eur_discovery_gwas.SYNTH_PHENO.glm.linear"
    if not gwas_glm_linear.exists():
        raise SystemExit(f"{gwas_glm_linear} not found — run src/trackA_synthetic/stage3_eur_gwas/run_gwas.py first.")

    out_dir = PRS_MODELS_ROOT / scenario_name / "prscsx"
    out_dir.mkdir(parents=True, exist_ok=True)

    sumstats_path = out_dir / "sumstats_rsid.tsv"
    n_matched = write_sumstats(gwas_glm_linear, rsid_map, sumstats_path)
    n_gwas = count_samples(EUR_DISCOVERY_KEEP_FILE)
    logging.info("%s: %d HapMap3-matched SNPs, n_gwas=%d", scenario_name, n_matched, n_gwas)

    out_name = "prscsx"
    run([
        sys.executable, PRSCSX_DIR / "PRScsx.py",
        f"--ref_dir={LD_REF_DIR}",
        f"--bim_prefix={target_bim_prefix}",
        f"--sst_file={sumstats_path}",
        f"--n_gwas={n_gwas}",
        f"--pop={POP}",
        f"--chrom={','.join(CHROMS)}",
        f"--out_dir={out_dir}",
        f"--out_name={out_name}",
        "--seed=20260801",
    ], cwd=PRSCSX_DIR)

    # Concatenate per-chromosome posterior effect files and translate rsIDs back to our IDs.
    dest = out_dir / "prscsx.scoring.tsv"
    with open(dest, "w") as fout:
        fout.write("ID\tA1\tBETA\n")
        for chrom in CHROMS:
            eff_file = out_dir / f"{out_name}_{POP}_pst_eff_a1_b0.5_phiauto_chr{chrom}.txt"
            if not eff_file.exists():
                logging.warning("%s: no output for chr%s (%s missing)", scenario_name, chrom, eff_file)
                continue
            with open(eff_file) as fin:
                for line in fin:
                    _chrom, rsid, _bp, a1, _a2, beta = line.split()
                    our_id = rsid_to_ours.get(rsid)
                    if our_id is None:
                        continue
                    fout.write(f"{our_id}\t{a1}\t{beta}\n")
    return dest


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenario", choices=["1", "2", "all"], default="all")
    p.add_argument("--download-ld-ref", action="store_true", help="Fetch the EUR LD reference panel + HapMap3 SNP info first")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.download_ld_ref:
        download_ld_reference("eur")

    if not SNPINFO_PATH.exists() or not (LD_REF_DIR / "ldblk_1kg_eur").exists():
        raise SystemExit(
            "PRS-CSx LD reference not found. Run with --download-ld-ref first "
            f"(expects {SNPINFO_PATH} and {LD_REF_DIR / 'ldblk_1kg_eur'})."
        )

    rsid_map = build_rsid_map()
    target_bim_prefix = PRS_MODELS_ROOT / "_shared" / "prscsx_target"
    rsid_to_ours = write_target_bim(rsid_map, target_bim_prefix)
    logging.info("%d HapMap3-matched target SNPs (chr %s)", len(rsid_to_ours), ",".join(CHROMS))

    scenarios = {"1": "scenario1_equal_effects", "2": "scenario2_ancestry_varying_effects"}
    scenario_ids = ["1", "2"] if args.scenario == "all" else [args.scenario]
    for sid in scenario_ids:
        dest = run_prscsx_for_scenario(scenarios[sid], rsid_map, rsid_to_ours, target_bim_prefix)
        logging.info("PRS-CSx complete for scenario %s: %s", sid, dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
