#!/usr/bin/env python3
"""
verify_tools.py — Phase 0 tool check/install for the PRS Cross-Ancestry Portability
Benchmark (see docs/BUILD_PLAN.md §5, §10, Phase 0 deliverable: "GCTA, PRSice-2,
LDpred2/bigsnpr, PRS-CSx, and pgscatalog-utils all installed and verified").

Everything installs into the user's own account — no sudo, no system package manager —
so it works in sandboxed/CI environments: static binaries go under ./tools/, Python
packages via `pip install --user`, R packages via `Rscript install.packages()`.

Usage:
    python3 scripts/verify_tools.py                 # check only, report status
    python3 scripts/verify_tools.py --install        # also install what's missing
    python3 scripts/verify_tools.py --install --with-r   # also attempt bigsnpr (slow: compiles from source)
"""

from __future__ import annotations

import argparse
import logging
import shutil
import stat
import subprocess
import sys
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

TOOLS_DIR = Path("tools")

GCTA_URL = "https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.3-linux-x86_64.zip"
PRSICE2_URL = "https://github.com/choishingwan/PRSice/releases/latest/download/PRSice_linux.zip"
PLINK2_URL = "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip"
PRSCSX_URL = "https://codeload.github.com/getian107/PRScsx/zip/refs/heads/master"


USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36"


def download_with_resume(url: str, dest: Path, *, max_retries: int = 8, timeout: float = 120.0) -> None:
    """Chunked, resumable, retrying download — this environment's bandwidth to some hosts
    is slow enough (tens of KB/s) that a single-shot request can stall or truncate."""
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    remote_size = None
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            cl = resp.headers.get("Content-Length")
            remote_size = int(cl) if cl is not None else None
    except urllib.error.URLError as exc:
        logging.warning("HEAD failed for %s (%s); will attempt GET anyway", url, exc)

    for attempt in range(1, max_retries + 1):
        try:
            existing = dest.stat().st_size if dest.exists() else 0
            if remote_size is not None and existing >= remote_size:
                return
            headers = {"User-Agent": USER_AGENT}
            mode = "wb"
            if existing:
                headers["Range"] = f"bytes={existing}-"
                mode = "ab"
            req = urllib.request.Request(url, headers=headers)
            logging.info("downloading (attempt %d/%d, resuming at %d bytes): %s", attempt, max_retries, existing, url)
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, mode) as f:
                while chunk := resp.read(1024 * 1024):
                    f.write(chunk)
            if remote_size is not None and dest.stat().st_size != remote_size:
                raise IOError(f"size mismatch: got {dest.stat().st_size}, expected {remote_size}")
            return
        except (urllib.error.URLError, IOError, TimeoutError) as exc:
            logging.warning("attempt %d failed for %s: %s", attempt, url, exc)
            if attempt == max_retries:
                raise
            time.sleep(min(2 ** attempt, 30))


def download_zip(url: str, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    zip_path = dest_dir / "download.zip"
    download_with_resume(url, zip_path)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(dest_dir)
    zip_path.unlink()
    return dest_dir


def make_executable(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


def check(name: str, path_or_none: str | None) -> bool:
    status = f"found: {path_or_none}" if path_or_none else "MISSING"
    logging.info("%-28s %s", name, status)
    return path_or_none is not None


def find_in_tools(*patterns: str) -> str | None:
    """Search ./tools for an executable *file* (not directory) matching any of the given
    filename patterns. GCTA/PRSice-2/etc. unzip into their own subdirectory, whose name can
    itself match a naive glob, so directories must be filtered out explicitly."""
    if not TOOLS_DIR.exists():
        return None
    for pattern in patterns:
        for p in sorted(TOOLS_DIR.rglob(pattern)):
            if p.is_file() and (p.stat().st_mode & stat.S_IXUSR):
                return str(p)
    return None


def verify_plink2(install: bool) -> bool:
    exe = shutil.which("plink2") or find_in_tools("plink2")
    if exe or not install:
        return check("plink2", exe)
    dest = download_zip(PLINK2_URL, TOOLS_DIR / "plink2")
    exe_path = dest / "plink2"
    if exe_path.exists():
        make_executable(exe_path)
    return check("plink2", str(exe_path) if exe_path.exists() else None)


def verify_gcta(install: bool) -> bool:
    exe = shutil.which("gcta64") or shutil.which("gcta") or find_in_tools("gcta64", "gcta")
    if exe or not install:
        return check("gcta", exe)
    dest = download_zip(GCTA_URL, TOOLS_DIR / "gcta")
    candidates = [p for p in sorted(dest.rglob("gcta64")) + sorted(dest.rglob("gcta")) if p.is_file()]
    exe_path = candidates[0] if candidates else None
    if exe_path:
        make_executable(exe_path)
    return check("gcta", str(exe_path) if exe_path else None)


def verify_prsice2(install: bool) -> bool:
    exe = shutil.which("PRSice_linux") or find_in_tools("PRSice_linux")
    if exe or not install:
        return check("prsice2", exe)
    dest = download_zip(PRSICE2_URL, TOOLS_DIR / "prsice2")
    exe_path = dest / "PRSice_linux"
    if exe_path.exists():
        make_executable(exe_path)
    return check("prsice2", str(exe_path) if exe_path.exists() else None)


def find_script_in_tools(name: str) -> str | None:
    """Like find_in_tools, but for interpreted scripts (run via `python3 script.py`, not
    directly executed) so no executable-bit check applies."""
    if not TOOLS_DIR.exists():
        return None
    matches = [p for p in sorted(TOOLS_DIR.rglob(name)) if p.is_file()]
    return str(matches[0]) if matches else None


def verify_prscsx(install: bool) -> bool:
    script = find_script_in_tools("PRScsx.py")
    if script or not install:
        return check("prscsx", script)
    download_zip(PRSCSX_URL, TOOLS_DIR / "prscsx")
    script = find_script_in_tools("PRScsx.py")
    return check("prscsx", script)


def verify_pgscatalog_utils(install: bool) -> bool:
    exe = shutil.which("pgscatalog-download")
    if exe or not install:
        return check("pgscatalog-utils", exe)
    subprocess.run([sys.executable, "-m", "pip", "install", "--user", "pgscatalog-utils"], check=True)
    exe = shutil.which("pgscatalog-download")
    return check("pgscatalog-utils", exe)


def r_user_library() -> Path:
    """The system R library (/usr/local/lib/R/site-library) isn't writable without sudo, so
    packages must go into a per-user library directory that we create ourselves."""
    return Path.home() / "R" / "library"


def verify_bigsnpr(install: bool, with_r: bool) -> bool:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return check("bigsnpr (R)", None)
    user_lib = r_user_library()
    check_expr = f'.libPaths(c("{user_lib}", .libPaths())); cat(requireNamespace("bigsnpr", quietly=TRUE))'
    result = subprocess.run([rscript, "-e", check_expr], capture_output=True, text=True)
    installed = result.stdout.strip() == "TRUE"
    if installed or not (install and with_r):
        return check("bigsnpr (R)", "installed" if installed else None)
    logging.info("installing bigsnpr via CRAN into %s (this compiles from source and can take several minutes)", user_lib)
    user_lib.mkdir(parents=True, exist_ok=True)
    install_expr = f'install.packages("bigsnpr", repos="https://cloud.r-project.org", lib="{user_lib}")'
    subprocess.run([rscript, "-e", install_expr], check=True)
    result = subprocess.run([rscript, "-e", check_expr], capture_output=True, text=True)
    installed = result.stdout.strip() == "TRUE"
    return check("bigsnpr (R)", "installed" if installed else None)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--install", action="store_true", help="Install anything missing")
    p.add_argument("--with-r", action="store_true", help="Also attempt to install bigsnpr (slow)")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    checks = {
        "plink2": lambda: verify_plink2(args.install),
        "gcta": lambda: verify_gcta(args.install),
        "prsice2": lambda: verify_prsice2(args.install),
        "prscsx": lambda: verify_prscsx(args.install),
        "pgscatalog-utils": lambda: verify_pgscatalog_utils(args.install),
        "bigsnpr (R)": lambda: verify_bigsnpr(args.install, args.with_r),
    }
    results = {}
    for name, fn in checks.items():
        try:
            results[name] = fn()
        except Exception as exc:
            logging.error("%-28s ERROR: %s", name, exc)
            results[name] = False

    missing = [name for name, ok in results.items() if not ok]
    if missing:
        logging.info("\nMissing: %s", ", ".join(missing))
        logging.info("Re-run with --install (and --with-r for bigsnpr) to attempt automatic setup.")
        return 1
    logging.info("\nAll tools verified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
