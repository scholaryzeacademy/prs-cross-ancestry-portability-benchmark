#!/usr/bin/env python3
"""
download_pgs.py — Stage B1: download real PGS Catalog scores (BUILD_PLAN.md §6 Stage B1).

Fetches real, published PGS Catalog scoring files via the official `pgscatalog-download` tool
(installed by scripts/verify_tools.py) — Track B is explicitly about applying REAL published
scores, not custom-built ones, as a descriptive cross-check against Track A's synthetic results.

Default: height only (PGS000297, "GRS3290_Height", Xie et al. 2020, 3,290 variants). BMI's most
commonly cited PGS Catalog entry (PGS000027) is a genome-wide score with ~2.1 million variants —
appropriate for a full-genome analysis, but wasteful to fully process against this project's
chr21+22-only smoke-test genotype panel, where the overwhelming majority of those variants
couldn't be scored anyway. Height alone satisfies BUILD_PLAN.md's "one or more" requirement with
a well-powered, appropriately-scoped score; add BMI (or another trait) via --pgs-id following the
same pattern once working against the full genome.

Scores are downloaded in GRCh37-harmonized form (`-b GRCh37`), matching 1000 Genomes Phase 3's
build, so Stage B2 can match on chromosome+position directly rather than needing an rsID lookup.

Usage:
    python3 download_pgs.py                        # height only (default)
    python3 download_pgs.py --pgs-id PGS000297 PGS000027
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = REPO_ROOT / "data" / "processed" / "pgs_catalog"

DEFAULT_PGS_IDS = ["PGS000297"]  # height (GRS3290_Height); see module docstring re: BMI


def require_pgscatalog_download() -> str:
    exe = shutil.which("pgscatalog-download")
    if exe is None:
        raise SystemExit("pgscatalog-download not found on PATH. Run scripts/verify_tools.py --install first.")
    return exe


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pgs-id", nargs="+", default=DEFAULT_PGS_IDS)
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    exe = require_pgscatalog_download()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "-i", *args.pgs_id, "-b", "GRCh37", "-o", str(OUT_DIR)]
    logging.info("+ %s", " ".join(cmd))
    subprocess.run(cmd, check=True)

    downloaded = sorted(OUT_DIR.glob("*.txt.gz"))
    logging.info("Stage B1 complete: %d scoring file(s) in %s", len(downloaded), OUT_DIR)
    for f in downloaded:
        logging.info("  %s", f.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
