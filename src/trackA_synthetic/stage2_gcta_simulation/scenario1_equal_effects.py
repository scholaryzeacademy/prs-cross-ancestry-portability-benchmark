#!/usr/bin/env python3
"""
scenario1_equal_effects.py — Stage 2, Scenario 1 (BUILD_PLAN.md §6 Stage 2).

Simulates a quantitative phenotype with GCTA using the *same* causal-variant effect sizes
applied uniformly across all five super-populations. This isolates the portability problem
caused purely by real differences in allele frequency and LD structure across ancestries,
holding the underlying biology (effect sizes) constant.

Runs GCTA once per super-population (same causal loci + effect sizes each time, --keep
restricting to that ancestry's samples), then concatenates the five .phen outputs into one
combined phenotype file.

Usage:
    python3 scenario1_equal_effects.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    config = common.load_config()
    plink2 = common.require_plink2()
    gcta = common.require_gcta()

    bed_prefix = common.prepare_bed(plink2)
    keep_files = common.prepare_gcta_keep_files()
    causal_variants = common.select_causal_variants(plink2, bed_prefix, config)

    base_seed = config["causal_variants"]["random_seed"]
    effects = common.generate_base_effects(causal_variants, seed=base_seed)

    scenario_dir = common.SIM_ROOT / "scenario1_equal_effects"
    causal_loci_file = scenario_dir / "causal_loci.txt"
    common.write_causal_loci_file(effects, causal_loci_file)

    hsq = config["trait"]["heritability"]
    phen_by_pop = {}
    for sp in common.SUPER_POPULATIONS:
        out_prefix = scenario_dir / f"gcta_{sp}"
        phen_by_pop[sp] = common.run_gcta_simu(
            gcta, bed_prefix, keep_files[sp], causal_loci_file, hsq, out_prefix,
        )

    common.concatenate_phenotypes(phen_by_pop, scenario_dir / "phenotypes.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
