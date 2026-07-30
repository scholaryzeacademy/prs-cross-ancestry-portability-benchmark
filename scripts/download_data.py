#!/usr/bin/env python3
"""
download_data.py — Automated data/tool fetcher for the PRS Cross-Ancestry
Portability & Recalibration Benchmark (see docs/BUILD_PLAN.md).

Every URL below was individually verified (HTTP HEAD / API lookup) before
being hard-coded here. Sources, per BUILD_PLAN.md §4:

  - 1000 Genomes Phase 3 genotypes + sample panel/pedigree     (IGSR/EBI FTP-over-HTTPS)
  - PGS Catalog scoring files for height (PGS000297) and BMI (PGS000027)
  - GCTA binary                                                 (Yang Lab)
  - PRSice-2 binary + source                                    (GitHub releases)
  - PRS-CSx source code + 1000G/UKBB LD reference panels         (GitHub + Dropbox)
  - bigsnpr/LDpred2 R package                                    (CRAN, via Rscript)

Design notes
------------
* "First download" / idempotency: a file already present with a size that
  matches the server's Content-Length is treated as already-downloaded and
  is skipped, so re-running the script is safe and cheap. Use --force to
  redownload anyway.
* Retry + resume: every download goes through a retry loop with exponential
  backoff + jitter. Partial files are resumed with an HTTP Range request
  where the server advertises support; otherwise the partial file is
  discarded and restarted. After each attempt the file size is checked
  against Content-Length (and against a sha256 when one is known) so a
  truncated/corrupted download is detected and retried rather than silently
  accepted as "broken".
* Size, by default: the full 1000 Genomes VCF set + full PRS-CSx LD panels
  are tens to ~100+ GB combined, so by default this script only fetches a
  lightweight smoke-test subset (chr21/chr22, no LD panels) plus all small
  metadata/tool files. Pass --chromosomes all and/or --prscsx-ld to opt into
  the full, heavy downloads.

Usage
-----
  python3 download_data.py --list
  python3 download_data.py                                   # default light set
  python3 download_data.py --categories all --chromosomes all --prscsx-ld eur
  python3 download_data.py --out-dir ./data --workers 4 --max-retries 6
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import hashlib
import json
import logging
import random
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Optional

DEFAULT_USER_AGENT = "prs-ancestry-benchmark-downloader/1.0"
CHUNK_SIZE = 1024 * 1024  # 1 MiB

# ---------------------------------------------------------------------------
# Catalog
# ---------------------------------------------------------------------------

ALL_CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
ALL_ANCESTRIES = ["afr", "amr", "eas", "eur", "sas"]

VCF_FILENAMES = {
    **{c: f"ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
       for c in [str(i) for i in range(1, 23)]},
    "X": "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz",
    "Y": "ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz",
    "MT": "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz",
}

KG_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
PGS_BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores"


@dataclasses.dataclass
class DownloadItem:
    key: str
    url: str
    dest: Path            # relative to --out-dir
    category: str
    sha256: Optional[str] = None
    optional: bool = False   # failure doesn't affect exit code
    post_process: Optional[str] = None  # "unzip" | "gunzip" | None


def dropbox_direct(url: str) -> str:
    """Force Dropbox share links to serve raw file bytes instead of the HTML preview page."""
    if "dropbox.com" not in url:
        return url
    if "dl=1" in url:
        return url
    if "dl=0" in url:
        return url.replace("dl=0", "dl=1")
    sep = "&" if "?" in url else "?"
    return f"{url}{sep}dl=1"


def build_catalog(args: argparse.Namespace) -> list[DownloadItem]:
    items: list[DownloadItem] = []
    categories = set(args.categories)
    want_all = "all" in categories

    # -- tools -------------------------------------------------------------
    if want_all or "tools" in categories:
        items += [
            DownloadItem(
                key="gcta",
                url="https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.3-linux-x86_64.zip",
                dest=Path("tools/gcta/gcta-1.95.3-linux-x86_64.zip"),
                category="tools",
                post_process="unzip",
            ),
            DownloadItem(
                key="prsice2_linux",
                url="https://github.com/choishingwan/PRSice/releases/latest/download/PRSice_linux.zip",
                dest=Path("tools/prsice2/PRSice_linux.zip"),
                category="tools",
                post_process="unzip",
            ),
            DownloadItem(
                key="prscsx_source",
                url="https://codeload.github.com/getian107/PRScsx/zip/refs/heads/master",
                dest=Path("tools/prscsx/PRScsx-master.zip"),
                category="tools",
                post_process="unzip",
            ),
        ]

    # -- 1000 Genomes metadata/panel (small, always cheap) ------------------
    if want_all or "panels" in categories:
        for fname in [
            "integrated_call_samples_v3.20130502.ALL.panel",
            "integrated_call_male_samples_v3.20130502.ALL.panel",
            "integrated_call_samples_v3.20200731.ALL.ped",
        ]:
            items.append(DownloadItem(
                key=f"1kg_panel_{fname}",
                url=f"{KG_BASE}/{fname}",
                dest=Path("1000genomes/metadata") / fname,
                category="panels",
            ))

    # -- 1000 Genomes VCFs (large; default subset only) ----------------------
    if want_all or "1000g-vcf" in categories:
        chroms = ALL_CHROMOSOMES if args.chromosomes == ["all"] else args.chromosomes
        for c in chroms:
            if c not in VCF_FILENAMES:
                logging.warning("Unknown chromosome %r, skipping", c)
                continue
            fname = VCF_FILENAMES[c]
            items.append(DownloadItem(
                key=f"1kg_vcf_chr{c}",
                url=f"{KG_BASE}/{fname}",
                dest=Path("1000genomes/vcf") / fname,
                category="1000g-vcf",
            ))
            items.append(DownloadItem(
                key=f"1kg_vcf_chr{c}_tbi",
                url=f"{KG_BASE}/{fname}.tbi",
                dest=Path("1000genomes/vcf") / f"{fname}.tbi",
                category="1000g-vcf",
                optional=True,
            ))

    # -- PGS Catalog scoring files -------------------------------------------
    if want_all or "pgs" in categories:
        for pgs_id, label in [("PGS000297", "height"), ("PGS000027", "bmi")]:
            items.append(DownloadItem(
                key=f"pgs_{pgs_id}_{label}",
                url=f"{PGS_BASE}/{pgs_id}/ScoringFiles/{pgs_id}.txt.gz",
                dest=Path("pgs_catalog") / f"{pgs_id}_{label}.txt.gz",
                category="pgs",
            ))

    # -- PRS-CSx LD reference panels (huge; opt-in via --prscsx-ld) ----------
    if args.prscsx_ld:
        ancestries = ALL_ANCESTRIES if args.prscsx_ld == ["all"] else args.prscsx_ld
        dropbox_1kg = {
            "afr": "https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=0",
            "amr": "https://www.dropbox.com/s/uv5ydr4uv528lca/ldblk_1kg_amr.tar.gz?dl=0",
            "eas": "https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=0",
            "eur": "https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0",
            "sas": "https://www.dropbox.com/s/hsm0qwgyixswdcv/ldblk_1kg_sas.tar.gz?dl=0",
        }
        for a in ancestries:
            if a not in dropbox_1kg:
                logging.warning("Unknown ancestry %r for PRS-CSx LD panel, skipping", a)
                continue
            items.append(DownloadItem(
                key=f"prscsx_ld_1kg_{a}",
                url=dropbox_direct(dropbox_1kg[a]),
                dest=Path("tools/prscsx/ld_reference/1kg") / f"ldblk_1kg_{a}.tar.gz",
                category="prscsx-ld",
                post_process="untar",
            ))
        items.append(DownloadItem(
            key="prscsx_snpinfo_1kg",
            url=dropbox_direct("https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3?dl=0"),
            dest=Path("tools/prscsx/ld_reference") / "snpinfo_mult_1kg_hm3",
            category="prscsx-ld",
        ))

    return items


# ---------------------------------------------------------------------------
# Download engine
# ---------------------------------------------------------------------------

class DownloadError(Exception):
    pass


def head_content_length(url: str, timeout: float) -> Optional[int]:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": DEFAULT_USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            length = resp.headers.get("Content-Length")
            return int(length) if length is not None else None
    except Exception:
        return None


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
            h.update(chunk)
    return h.hexdigest()


def already_downloaded(item: DownloadItem, dest_path: Path, timeout: float) -> bool:
    if not dest_path.exists():
        return False
    if item.sha256:
        return sha256_of(dest_path) == item.sha256
    remote_size = head_content_length(item.url, timeout)
    if remote_size is None:
        # Can't verify remotely; trust a non-empty existing file.
        return dest_path.stat().st_size > 0
    return dest_path.stat().st_size == remote_size


def stream_download(url: str, dest_path: Path, timeout: float,
                     progress_label: str) -> None:
    tmp_path = dest_path.with_suffix(dest_path.suffix + ".part")
    resume_from = tmp_path.stat().st_size if tmp_path.exists() else 0

    headers = {"User-Agent": DEFAULT_USER_AGENT}
    if resume_from:
        headers["Range"] = f"bytes={resume_from}-"

    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            mode = "ab" if resume_from and resp.status == 206 else "wb"
            if mode == "wb" and resume_from:
                resume_from = 0  # server ignored Range; restart clean
            total = resp.headers.get("Content-Length")
            total = int(total) + resume_from if total is not None else None

            dest_path.parent.mkdir(parents=True, exist_ok=True)
            downloaded = resume_from
            last_log = time.monotonic()
            with open(tmp_path, mode) as f:
                while True:
                    chunk = resp.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    now = time.monotonic()
                    if now - last_log > 5:
                        pct = f"{downloaded/total*100:5.1f}%" if total else f"{downloaded/1e6:.0f} MB"
                        logging.info("  [%s] %s", progress_label, pct)
                        last_log = now
    except urllib.error.HTTPError as e:
        if e.code == 416:  # Range not satisfiable -> file already complete
            pass
        else:
            raise DownloadError(f"HTTP {e.code} for {url}") from e
    except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
        raise DownloadError(f"{type(e).__name__}: {e}") from e

    if not tmp_path.exists():
        raise DownloadError(f"No data written for {url}")
    tmp_path.rename(dest_path)


def verify(item: DownloadItem, dest_path: Path, timeout: float) -> None:
    if item.sha256:
        actual = sha256_of(dest_path)
        if actual != item.sha256:
            raise DownloadError(f"sha256 mismatch: expected {item.sha256}, got {actual}")
        return
    remote_size = head_content_length(item.url, timeout)
    if remote_size is not None and dest_path.stat().st_size != remote_size:
        raise DownloadError(
            f"size mismatch: expected {remote_size} bytes, got {dest_path.stat().st_size}"
        )
    if dest_path.stat().st_size == 0:
        raise DownloadError("downloaded file is empty")


def post_process(item: DownloadItem, dest_path: Path) -> None:
    if item.post_process == "unzip":
        target_dir = dest_path.parent
        try:
            with zipfile.ZipFile(dest_path) as zf:
                zf.extractall(target_dir)
        except zipfile.BadZipFile as e:
            raise DownloadError(f"corrupt zip archive: {e}") from e
    elif item.post_process == "gunzip":
        import gzip
        out_path = dest_path.with_suffix("")
        with gzip.open(dest_path, "rb") as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    elif item.post_process == "untar":
        import tarfile
        target_dir = dest_path.parent
        try:
            with tarfile.open(dest_path) as tf:
                tf.extractall(target_dir)
        except tarfile.TarError as e:
            raise DownloadError(f"corrupt tar archive: {e}") from e


def download_one(item: DownloadItem, out_dir: Path, args: argparse.Namespace) -> dict:
    dest_path = out_dir / item.dest
    result = {"key": item.key, "url": item.url, "dest": str(dest_path)}

    if not args.force and already_downloaded(item, dest_path, args.timeout):
        logging.info("SKIP  %-35s already downloaded -> %s", item.key, dest_path)
        result["status"] = "skipped"
        return result

    last_error: Optional[Exception] = None
    for attempt in range(1, args.max_retries + 1):
        try:
            logging.info("GET   %-35s (attempt %d/%d) %s",
                         item.key, attempt, args.max_retries, item.url)
            stream_download(item.url, dest_path, args.timeout, item.key)
            verify(item, dest_path, args.timeout)
            if item.post_process:
                post_process(item, dest_path)
            logging.info("OK    %-35s -> %s", item.key, dest_path)
            result["status"] = "ok"
            return result
        except DownloadError as e:
            last_error = e
            logging.warning("FAIL  %-35s attempt %d/%d: %s",
                            item.key, attempt, args.max_retries, e)
            # Corrupted/mismatched full file: drop it so the next attempt restarts clean.
            if dest_path.exists() and "mismatch" in str(e):
                dest_path.unlink(missing_ok=True)
            if attempt < args.max_retries:
                backoff = min(60, 2 ** attempt) + random.uniform(0, 1)
                logging.info("      retrying in %.1fs...", backoff)
                time.sleep(backoff)

    severity = logging.WARNING if item.optional else logging.ERROR
    logging.log(severity, "GIVE UP %-35s after %d attempts: %s",
                item.key, args.max_retries, last_error)
    result["status"] = "failed_optional" if item.optional else "failed"
    result["error"] = str(last_error)
    return result


# ---------------------------------------------------------------------------
# Extra, non-URL install step: bigsnpr/LDpred2 via R
# ---------------------------------------------------------------------------

def install_bigsnpr(out_dir: Path, args: argparse.Namespace) -> dict:
    key = "bigsnpr_ldpred2_r_package"
    if not shutil.which("Rscript"):
        logging.warning("SKIP  %-35s Rscript not found on PATH; install R to fetch LDpred2 "
                        "(CRAN package 'bigsnpr')", key)
        return {"key": key, "status": "skipped_no_r"}

    log_path = out_dir / "tools" / "bigsnpr_install.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    r_expr = (
        'if (!requireNamespace("bigsnpr", quietly=TRUE)) '
        'install.packages("bigsnpr", repos="https://cloud.r-project.org")'
    )
    last_error = None
    for attempt in range(1, args.max_retries + 1):
        logging.info("GET   %-35s (attempt %d/%d) via Rscript/CRAN", key, attempt, args.max_retries)
        try:
            with open(log_path, "ab") as logf:
                subprocess.run(
                    ["Rscript", "-e", r_expr],
                    stdout=logf, stderr=subprocess.STDOUT,
                    timeout=args.timeout * 20, check=True,
                )
            logging.info("OK    %-35s (log: %s)", key, log_path)
            return {"key": key, "status": "ok", "log": str(log_path)}
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            last_error = e
            logging.warning("FAIL  %-35s attempt %d/%d: %s", key, attempt, args.max_retries, e)
            if attempt < args.max_retries:
                time.sleep(min(60, 2 ** attempt))
    logging.warning("GIVE UP %-35s after %d attempts: %s", key, args.max_retries, last_error)
    return {"key": key, "status": "failed_optional", "error": str(last_error)}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

CATEGORY_CHOICES = ["tools", "panels", "pgs", "1000g-vcf", "all"]


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Download all datasets/tools for the PRS cross-ancestry benchmark.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--out-dir", type=Path, default=Path("data"),
                   help="Root directory to download everything into")
    p.add_argument("--categories", nargs="+", default=["tools", "panels", "pgs"],
                   choices=CATEGORY_CHOICES,
                   help="Which light-weight categories to fetch. '1000g-vcf' adds genotype "
                        "VCFs (large); use --chromosomes to scope it. 'all' includes every "
                        "light-weight category (still excludes --prscsx-ld, which is opt-in).")
    p.add_argument("--chromosomes", nargs="+", default=["21", "22"],
                   help="Chromosomes to fetch when '1000g-vcf'/'all' is selected. "
                        "Pass 'all' for the full genome (large, ~20GB).")
    p.add_argument("--prscsx-ld", nargs="*", default=None,
                   choices=ALL_ANCESTRIES + ["all"],
                   help="Opt-in: also fetch PRS-CSx 1000G LD reference panel(s) for these "
                        "ancestries (each ~4-6GB). Omit entirely to skip (default).")
    p.add_argument("--with-r-packages", action="store_true",
                   help="Also install the bigsnpr/LDpred2 R package via Rscript+CRAN")
    p.add_argument("--workers", type=int, default=4, help="Parallel downloads")
    p.add_argument("--max-retries", type=int, default=5, help="Retries per file before giving up")
    p.add_argument("--timeout", type=float, default=60.0, help="Per-request timeout (seconds)")
    p.add_argument("--force", action="store_true", help="Redownload even if already present")
    p.add_argument("--list", action="store_true", help="List the resolved catalog and exit")
    p.add_argument("--log-file", type=Path, default=None)
    return p.parse_args(argv)


def setup_logging(log_file: Optional[Path]) -> None:
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s",
                        handlers=handlers)


def main(argv=None) -> int:
    args = parse_args(argv)
    setup_logging(args.log_file)

    items = build_catalog(args)
    if args.list:
        for it in items:
            flag = " (optional)" if it.optional else ""
            print(f"[{it.category:10s}] {it.key:30s} {it.dest}{flag}\n    {it.url}")
        if args.with_r_packages:
            print("[r-package ] bigsnpr_ldpred2_r_package (via Rscript/CRAN, no direct URL)")
        return 0

    if not items and not args.with_r_packages:
        logging.error("Nothing to download for the selected categories: %s", args.categories)
        return 2

    args.out_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []

    if items:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(download_one, it, args.out_dir, args): it for it in items}
            for fut in concurrent.futures.as_completed(futures):
                results.append(fut.result())

    if args.with_r_packages:
        results.append(install_bigsnpr(args.out_dir, args))

    manifest_path = args.out_dir / "download_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(results, f, indent=2)

    ok = sum(1 for r in results if r["status"] == "ok")
    skipped = sum(1 for r in results if r["status"] == "skipped")
    failed = [r for r in results if r["status"] == "failed"]
    failed_optional = [r for r in results if r["status"] in ("failed_optional", "skipped_no_r")]

    logging.info("=" * 70)
    logging.info("Done. ok=%d skipped=%d failed=%d failed_optional=%d",
                 ok, skipped, len(failed), len(failed_optional))
    logging.info("Manifest written to %s", manifest_path)

    if failed:
        logging.error("Required downloads failed: %s", ", ".join(r["key"] for r in failed))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
