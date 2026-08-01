#!/usr/bin/env python3
"""
compare.py — Stage B3: descriptive cross-check against published statistics (BUILD_PLAN.md §6
Stage B3).

****************************************************************************
* DESCRIPTIVE / ILLUSTRATIVE ONLY — NOT A VALIDATED ACCURACY MEASUREMENT.  *
* 1000 Genomes carries no phenotype data, so there is no ground-truth      *
* height to validate a prediction against here. This compares the         *
* computed PGS's *relative ordering* across super-populations against     *
* published, external, population-level height statistics as a plausible- *
* ordering check only — see BUILD_PLAN.md §6 Stage B3 and §9 item 4.      *
****************************************************************************

Compares Stage B2's per-super-population mean PGS to published national height statistics
(published_height_reference.tsv, sourced from national surveys — see that file's header for
full citations) for representative countries per super-population. This is necessarily an
approximate comparison: 1000 Genomes super-populations aggregate several specific source
populations that don't map to a single country, the reference statistics span different survey
years/methods, and the PGS itself uses only 52 of 3,290 scoring variants (chr21+22
smoke-test scope) with no ancestry-appropriate calibration applied (see Stage 6 for why that
matters). A mismatch between PGS ranking and published height ranking is not evidence the tool
chain is broken — it can be exactly the real, well-documented PRS cross-ancestry portability
problem this entire project (Track A) is built to quantify with a controlled ground truth.

Usage:
    python3 compare.py
"""

from __future__ import annotations

import csv
import logging
import statistics
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
PGS_SCORES_PATH = REPO_ROOT / "data" / "processed" / "pgs_catalog" / "PGS000297_scores_by_ancestry.tsv"
REFERENCE_PATH = Path(__file__).resolve().parent / "published_height_reference.tsv"
OUT_PATH = REPO_ROOT / "data" / "processed" / "pgs_catalog" / "stageB3_descriptive_comparison.tsv"

SUPER_POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def load_pgs_means() -> dict:
    groups = defaultdict(list)
    with open(PGS_SCORES_PATH) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            groups[row["super_pop"]].append(float(row["PGS"]))
    return {sp: (statistics.mean(v), statistics.stdev(v), len(v)) for sp, v in groups.items()}


def load_reference_means() -> dict:
    groups = defaultdict(list)
    with open(REFERENCE_PATH) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "super_pop":
                continue
            sp, country, sex, height_cm = fields[0], fields[1], fields[2], fields[3]
            groups[sp].append(float(height_cm))
    return {sp: statistics.mean(v) for sp, v in groups.items()}


def rank(values: dict) -> dict:
    ordered = sorted(values.items(), key=lambda kv: kv[1], reverse=True)
    return {sp: i + 1 for i, (sp, _) in enumerate(ordered)}


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if not PGS_SCORES_PATH.exists():
        raise SystemExit(f"{PGS_SCORES_PATH} not found — run stageB2_compute_scores/compute_scores.py first.")

    pgs_means = load_pgs_means()
    ref_means = load_reference_means()

    pgs_rank = rank({sp: m for sp, (m, _, _) in pgs_means.items()})
    ref_rank = rank(ref_means)

    rows = []
    for sp in SUPER_POPS:
        pgs_mean, pgs_sd, n = pgs_means[sp]
        rows.append({
            "super_pop": sp, "n": n,
            "pgs_mean": round(pgs_mean, 4), "pgs_sd": round(pgs_sd, 4), "pgs_rank": pgs_rank[sp],
            "published_height_cm_mean": round(ref_means[sp], 1), "published_height_rank": ref_rank[sp],
            "rank_matches": pgs_rank[sp] == ref_rank[sp],
        })

    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    n_match = sum(r["rank_matches"] for r in rows)
    logging.info("DESCRIPTIVE ONLY — not a validated accuracy measurement (see module docstring).")
    for r in rows:
        logging.info(
            "%s: PGS mean=%.3f (rank %d/5) vs. published height=%.1fcm (rank %d/5) — %s",
            r["super_pop"], r["pgs_mean"], r["pgs_rank"],
            r["published_height_cm_mean"], r["published_height_rank"],
            "MATCH" if r["rank_matches"] else "differs",
        )
    logging.info(
        "%d/5 super-populations' PGS rank matches published height rank. A low match count here "
        "illustrates real PRS cross-ancestry portability limitations (Track A's whole subject), "
        "not necessarily a pipeline error — see this script's docstring and METHODS.md.",
        n_match,
    )
    logging.info("Stage B3 complete: %s", OUT_PATH)
    return 0


if __name__ == "__main__":
    sys.exit(main())
