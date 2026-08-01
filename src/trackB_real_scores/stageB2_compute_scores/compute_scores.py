#!/usr/bin/env python3
"""
compute_scores.py — Stage B2: compute real PGS Catalog scores on 1000 Genomes (BUILD_PLAN.md
§6 Stage B2).

Matches Stage B1's GRCh37-harmonized PGS Catalog scoring file to this project's QC'd 1000
Genomes genotype panel by chromosome+position (not rsID — Stage 1 assigned chr:pos:ref:alt IDs),
builds a `plink2 --score`-compatible scoring file, then computes the score for **all** 2,504
individuals across all five super-populations — unlike Track A's Stage 5, there's no held-out
split here, because Track B is descriptive only (no real phenotype to guard against overfitting
to in Track B's own use — see Stage B3 and BUILD_PLAN.md §6 Stage B3's repeated-disclaimer
requirement).

Usage:
    python3 compute_scores.py                       # PGS000297 (height), default from Stage B1
    python3 compute_scores.py --pgs-id PGS000297
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
SHARED_DIR = REPO_ROOT / "data" / "processed" / "simulated_phenotypes" / "_shared"
BED_PREFIX = SHARED_DIR / "1000g_qc_biallelic"
BIM_PATH = BED_PREFIX.with_suffix(".bim")
PGS_DIR = REPO_ROOT / "data" / "processed" / "pgs_catalog"
OUT_DIR = REPO_ROOT / "data" / "processed" / "pgs_catalog"

PANEL_PATH = REPO_ROOT / "data" / "1000genomes" / "metadata" / "integrated_call_samples_v3.20130502.ALL.panel"


def require_plink2() -> str:
    exe = shutil.which("plink2")
    if exe is None:
        raise SystemExit("plink2 not found on PATH. Run scripts/verify_tools.py --install first.")
    return exe


def load_position_index() -> dict:
    """{(chr, pos): (our_id, a1, a2)} from the shared QC'd bim."""
    index = {}
    with open(BIM_PATH) as f:
        for line in f:
            chrom, our_id, gd, pos, a1, a2 = line.split()
            index[(chrom, pos)] = (our_id, a1, a2)
    return index


def build_scoring_file(pgs_scoring_path: Path, position_index: dict, dest: Path) -> tuple:
    dest.parent.mkdir(parents=True, exist_ok=True)
    n_total, n_matched = 0, 0
    with gzip.open(pgs_scoring_path, "rt") as fin, open(dest, "w") as fout:
        header = None
        fout.write("ID\tEFFECT_ALLELE\tEFFECT_WEIGHT\n")
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                idx = {name: i for i, name in enumerate(header)}
                continue
            n_total += 1
            chrom, pos = fields[idx["hm_chr"]], fields[idx["hm_pos"]]
            effect_allele = fields[idx["effect_allele"]]
            weight = fields[idx["effect_weight"]]
            match = position_index.get((chrom, pos))
            if match is None:
                continue
            our_id, a1, a2 = match
            if effect_allele not in (a1, a2):
                continue  # allele mismatch at this position (e.g. multiallelic/strand issue) — skip rather than guess
            fout.write(f"{our_id}\t{effect_allele}\t{weight}\n")
            n_matched += 1
    return n_total, n_matched


def load_super_pop_labels() -> dict:
    labels = {}
    with open(PANEL_PATH, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            labels[row["sample"]] = row["super_pop"]
    return labels


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pgs-id", default="PGS000297")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    plink2 = require_plink2()
    pgs_scoring_path = PGS_DIR / f"{args.pgs_id}_hmPOS_GRCh37.txt.gz"
    if not pgs_scoring_path.exists():
        raise SystemExit(f"{pgs_scoring_path} not found — run stageB1_download_pgs/download_pgs.py first.")

    position_index = load_position_index()
    scoring_file = OUT_DIR / f"{args.pgs_id}.scoring.tsv"
    n_total, n_matched = build_scoring_file(pgs_scoring_path, position_index, scoring_file)
    logging.info("%s: %d/%d variants matched to the QC'd genotype panel (chr21+22 smoke-test scope)",
                 args.pgs_id, n_matched, n_total)
    if n_matched == 0:
        raise SystemExit("no variants matched — cannot compute a score")

    out_prefix = OUT_DIR / f"{args.pgs_id}_scores"
    subprocess.run([
        plink2, "--bfile", BED_PREFIX,
        "--score", scoring_file, "1", "2", "3", "header", "cols=+scoresums",
        "--out", out_prefix,
    ], check=True)

    super_pop = load_super_pop_labels()
    dest = OUT_DIR / f"{args.pgs_id}_scores_by_ancestry.tsv"
    with open(f"{out_prefix}.sscore") as fin, open(dest, "w") as fout:
        header = fin.readline().lstrip("#").split()
        idx = {name: i for i, name in enumerate(header)}
        fout.write("IID\tsuper_pop\tPGS\n")
        for line in fin:
            fields = line.split()
            iid = fields[idx["IID"]]
            fout.write(f"{iid}\t{super_pop.get(iid, 'NA')}\t{fields[idx['SCORE1_SUM']]}\n")
    logging.info("Stage B2 complete: %s", dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
