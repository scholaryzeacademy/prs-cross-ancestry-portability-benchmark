#!/usr/bin/env python3
"""
qc.py — Stage 1 QC pipeline for 1000 Genomes Phase 3 (see docs/BUILD_PLAN.md §6, Stage 1).

Given VCFs downloaded by download.py and the IGSR sample panel, this:
  1. Converts each chromosome VCF to PLINK2 pgen/pvar/psam.
  2. Merges chromosomes into one dataset.
  3. Builds a --keep sample list per super-population (AFR, AMR, EAS, EUR, SAS) from the panel.
  4. Applies QC filters (--maf, --geno, --hwe) *within each super-population separately* and
     keeps only variants that pass in every super-population — a single global filter would
     hide population-specific HWE/frequency problems, which is exactly what the build plan's
     "QC filters ... within each super-population" requirement is guarding against.
  5. Applies a global per-sample missingness filter (--mind) and writes the final QC'd pgen.

Requires PLINK2 (https://www.cog-genomics.org/plink/2.0/) on PATH. Run
`scripts/verify_tools.py` first to check/install it.

Usage:
    python3 qc.py --vcf-dir data/1000genomes/vcf --panel data/1000genomes/metadata/integrated_call_samples_v3.20130502.ALL.panel \
        --chromosomes 21 22 --out-dir data/processed/1000g_qc
"""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import sys
from pathlib import Path

SUPER_POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def run(cmd: list[str]) -> None:
    logging.info("+ %s", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True)


def require_plink2() -> str:
    exe = shutil.which("plink2")
    if exe is None:
        raise SystemExit(
            "plink2 not found on PATH. Install it (see scripts/verify_tools.py) before running QC."
        )
    return exe


def load_panel(panel_path: Path) -> dict[str, list[str]]:
    """Return {super_pop: [sample_id, ...]} from the IGSR panel file."""
    groups: dict[str, list[str]] = {sp: [] for sp in SUPER_POPULATIONS}
    with open(panel_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sp = row.get("super_pop", "").strip()
            sample = row.get("sample", "").strip()
            if sp in groups and sample:
                groups[sp].append(sample)
    missing = [sp for sp, ids in groups.items() if not ids]
    if missing:
        raise ValueError(f"panel file had no samples for super-population(s): {missing}")
    return groups


def write_keep_file(sample_ids: list[str], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w") as f:
        for sid in sample_ids:
            # PLINK2's VCF-derived pgen has no #FID column (IID-only mode), so --keep must be
            # given a single-column ID file — a two-column FID+IID file matches zero samples.
            f.write(f"{sid}\n")


def convert_and_merge(plink2: str, vcf_dir: Path, chromosomes: list[str], work_dir: Path) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    per_chrom_prefixes = []
    for c in chromosomes:
        vcfs = sorted(vcf_dir.glob(f"ALL.chr{c}.*.vcf.gz"))
        if not vcfs:
            raise FileNotFoundError(f"no VCF found for chromosome {c} in {vcf_dir}")
        prefix = work_dir / f"chr{c}"
        run([
            plink2, "--vcf", str(vcfs[0]), "--make-pgen", "--out", str(prefix),
            "--set-all-var-ids", "@:#:$r:$a",
            # 1000 Genomes VCFs include long indels/structural variants whose allele codes
            # exceed plink2's default --set-all-var-ids length cap; variants still over 1000
            # chars (very rare SVs) get a '.' ID instead of failing the whole conversion.
            "--new-id-max-allele-len", "1000", "missing",
        ])
        per_chrom_prefixes.append(prefix)

    if len(per_chrom_prefixes) == 1:
        return per_chrom_prefixes[0]

    merge_list = work_dir / "merge_list.txt"
    with open(merge_list, "w") as f:
        for prefix in per_chrom_prefixes[1:]:
            f.write(f"{prefix}\n")
    merged = work_dir / "merged_raw"
    run([
        plink2, "--pfile", str(per_chrom_prefixes[0]),
        "--pmerge-list", str(merge_list),
        "--make-pgen", "--out", str(merged),
    ])
    return merged


def per_superpop_passing_variants(
    plink2: str, merged_prefix: Path, panel_groups: dict[str, list[str]],
    work_dir: Path, maf: float, geno: float, hwe: float,
) -> Path:
    passing_lists = []
    for sp in SUPER_POPULATIONS:
        keep_file = work_dir / f"keep_{sp}.txt"
        write_keep_file(panel_groups[sp], keep_file)
        sp_prefix = work_dir / f"qc_{sp}"
        run([
            plink2, "--pfile", str(merged_prefix),
            "--keep", str(keep_file),
            "--maf", str(maf),
            "--geno", str(geno),
            "--hwe", str(hwe),
            "--write-snplist",
            "--out", str(sp_prefix),
        ])
        passing_lists.append(sp_prefix.with_suffix(".snplist"))

    variant_sets = [set(p.read_text().split()) for p in passing_lists]
    shared = set.intersection(*variant_sets)
    logging.info(
        "variants passing QC per super-pop: %s; shared across all 5: %d",
        {sp: len(v) for sp, v in zip(SUPER_POPULATIONS, variant_sets)}, len(shared),
    )
    shared_path = work_dir / "variants_pass_all_superpops.txt"
    shared_path.write_text("\n".join(sorted(shared)) + "\n")
    return shared_path


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf-dir", type=Path, required=True)
    p.add_argument("--panel", type=Path, required=True)
    p.add_argument("--chromosomes", nargs="+", default=["21", "22"])
    p.add_argument("--out-dir", type=Path, default=Path("data/processed/1000g_qc"))
    p.add_argument("--maf", type=float, default=0.01, help="Minor allele frequency threshold (default: %(default)s)")
    p.add_argument("--geno", type=float, default=0.05, help="Per-variant missingness threshold (default: %(default)s)")
    p.add_argument("--hwe", type=float, default=1e-6, help="Hardy-Weinberg equilibrium p-value threshold (default: %(default)s)")
    p.add_argument("--mind", type=float, default=0.05, help="Per-sample missingness threshold (default: %(default)s)")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    plink2 = require_plink2()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = args.out_dir / "_work"

    panel_groups = load_panel(args.panel)
    logging.info("sample counts by super-population: %s", {sp: len(v) for sp, v in panel_groups.items()})

    merged_prefix = convert_and_merge(plink2, args.vcf_dir, args.chromosomes, work_dir)
    shared_variants = per_superpop_passing_variants(
        plink2, merged_prefix, panel_groups, work_dir, args.maf, args.geno, args.hwe,
    )

    final_prefix = args.out_dir / "1000g_qc"
    run([
        plink2, "--pfile", str(merged_prefix),
        "--extract", str(shared_variants),
        "--mind", str(args.mind),
        "--make-pgen", "--out", str(final_prefix),
    ])

    for sp, ids in panel_groups.items():
        write_keep_file(ids, args.out_dir / f"samples_{sp}.txt")

    logging.info("Stage 1 QC complete: %s.{pgen,pvar,psam}", final_prefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())
