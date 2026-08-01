#!/usr/bin/env python3
"""
scenario2_ancestry_varying_effects.py — Stage 2, Scenario 2 (BUILD_PLAN.md §6 Stage 2).

Simulates a quantitative phenotype with GCTA using the *same causal variants as Scenario 1*
but with effect sizes perturbed independently per super-population before simulation, modeling
genuine gene-by-ancestry effect-size heterogeneity (documented in
configs/simulation_parameters.yaml under scenarios.scenario2_ancestry_varying_effects):

    effect_size[ancestry] = base_effect * N(1, perturbation_sd)

applied with a per-ancestry random seed offset so the perturbation is deterministic and
reproducible. Runs GCTA once per super-population (its own perturbed causal-loci file, --keep
restricting to that ancestry's samples), then concatenates the five .phen outputs.

Usage:
    python3 scenario2_ancestry_varying_effects.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common


def perturb_effects(base_effects: dict, perturbation_sd: float, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    noise = rng.normal(loc=1.0, scale=perturbation_sd, size=len(base_effects))
    return {snp: effect * n for (snp, effect), n in zip(base_effects.items(), noise)}


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    config = common.load_config()
    plink2 = common.require_plink2()
    gcta = common.require_gcta()

    bed_prefix = common.prepare_bed(plink2)
    keep_files = common.prepare_gcta_keep_files()
    causal_variants = common.select_causal_variants(plink2, bed_prefix, config)

    base_seed = config["causal_variants"]["random_seed"]
    base_effects = common.generate_base_effects(causal_variants, seed=base_seed)

    perturbation = config["scenarios"]["scenario2_ancestry_varying_effects"]["per_ancestry_perturbation"]
    perturbation_sd = perturbation["perturbation_sd"]
    seed_offsets = perturbation["per_ancestry_random_seed_offset"]

    scenario_dir = common.SIM_ROOT / "scenario2_ancestry_varying_effects"
    hsq = config["trait"]["heritability"]
    phen_by_pop = {}
    for sp in common.SUPER_POPULATIONS:
        pop_seed = base_seed + seed_offsets[sp]
        pop_effects = perturb_effects(base_effects, perturbation_sd, seed=pop_seed)
        causal_loci_file = scenario_dir / f"causal_loci_{sp}.txt"
        common.write_causal_loci_file(pop_effects, causal_loci_file)

        out_prefix = scenario_dir / f"gcta_{sp}"
        phen_by_pop[sp] = common.run_gcta_simu(
            gcta, bed_prefix, keep_files[sp], causal_loci_file, hsq, out_prefix,
        )

    common.concatenate_phenotypes(phen_by_pop, scenario_dir / "phenotypes.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
