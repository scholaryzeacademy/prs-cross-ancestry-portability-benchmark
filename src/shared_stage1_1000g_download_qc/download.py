#!/usr/bin/env python3
"""
download.py — Stage 1 downloader for 1000 Genomes Phase 3 (see docs/BUILD_PLAN.md §6, Stage 1).

Fetches the sample panel (population / super-population labels, required to build
per-super-population QC groups and the EUR-only discovery GWAS cohort in Stage 3) and
per-chromosome VCFs from the IGSR/EBI FTP-over-HTTPS mirror. Stdlib-only, no third-party
dependencies. Downloads are idempotent: a file whose size already matches the server's
Content-Length is skipped unless --force is given.

Default is a smoke-test subset (chr21, chr22) — the full autosomal set is tens of GB and
should be fetched deliberately with --chromosomes all.

Usage:
    python3 download.py --list
    python3 download.py                                   # panel + chr21, chr22
    python3 download.py --chromosomes 1 2 3
    python3 download.py --chromosomes all
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

KG_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

VCF_FILENAMES = {
    **{c: f"ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
       for c in [str(i) for i in range(1, 23)]},
    "X": "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz",
}

PANEL_FILES = [
    "integrated_call_samples_v3.20130502.ALL.panel",
]

DEFAULT_CHROMOSOMES = ["21", "22"]


def fetch(url: str, dest: Path, *, max_retries: int = 5, timeout: float = 60.0, force: bool = False) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)

    remote_size = None
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            cl = resp.headers.get("Content-Length")
            remote_size = int(cl) if cl is not None else None
    except urllib.error.URLError as exc:
        logging.warning("HEAD failed for %s (%s); will attempt GET anyway", url, exc)

    if not force and dest.exists() and remote_size is not None and dest.stat().st_size == remote_size:
        logging.info("skip (already downloaded): %s", dest)
        return

    for attempt in range(1, max_retries + 1):
        try:
            logging.info("downloading (attempt %d/%d): %s -> %s", attempt, max_retries, url, dest)
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as f:
                while chunk := resp.read(1024 * 1024):
                    f.write(chunk)
            if remote_size is not None and dest.stat().st_size != remote_size:
                raise IOError(f"size mismatch: got {dest.stat().st_size}, expected {remote_size}")
            logging.info("done: %s (%d bytes)", dest, dest.stat().st_size)
            return
        except (urllib.error.URLError, IOError, TimeoutError) as exc:
            logging.warning("attempt %d failed for %s: %s", attempt, url, exc)
            if attempt == max_retries:
                raise
            time.sleep(min(2 ** attempt, 30))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, default=Path("data/1000genomes"))
    p.add_argument("--chromosomes", nargs="+", default=DEFAULT_CHROMOSOMES,
                   help="Chromosome list, or 'all' for 1-22 (default: %(default)s)")
    p.add_argument("--force", action="store_true")
    p.add_argument("--list", action="store_true", help="Print resolved file list and exit")
    p.add_argument("--max-retries", type=int, default=5)
    p.add_argument("--timeout", type=float, default=60.0)
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    chroms = list(VCF_FILENAMES.keys()) if args.chromosomes == ["all"] else args.chromosomes
    unknown = [c for c in chroms if c not in VCF_FILENAMES]
    if unknown:
        logging.error("unknown chromosome(s): %s (known: %s)", unknown, sorted(VCF_FILENAMES, key=str))
        return 2

    jobs = [(f"{KG_BASE}/{f}", args.out_dir / "metadata" / f) for f in PANEL_FILES]
    for c in chroms:
        fname = VCF_FILENAMES[c]
        jobs.append((f"{KG_BASE}/{fname}", args.out_dir / "vcf" / fname))
        jobs.append((f"{KG_BASE}/{fname}.tbi", args.out_dir / "vcf" / f"{fname}.tbi"))

    if args.list:
        for url, dest in jobs:
            print(f"{url}  ->  {dest}")
        return 0

    failed = []
    for url, dest in jobs:
        try:
            fetch(url, dest, max_retries=args.max_retries, timeout=args.timeout, force=args.force)
        except Exception as exc:
            logging.error("giving up on %s: %s", url, exc)
            failed.append(url)

    if failed:
        logging.error("%d file(s) failed: %s", len(failed), failed)
        return 1
    logging.info("all files present under %s", args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
