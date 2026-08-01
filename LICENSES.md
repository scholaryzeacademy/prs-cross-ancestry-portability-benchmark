# LICENSES.md

Tracks licensing status per `docs/BUILD_PLAN.md` §11. Original code in this repository is
licensed under MIT (see `LICENSE`) — a portfolio-appropriate default; revisit if the project
owner prefers Apache-2.0. Third-party data and tools each carry their own license, checked
individually below rather than assumed uniform.

## Checklist (BUILD_PLAN.md §11)

- [ ] Confirm 1000 Genomes/IGSR data usage terms
- [ ] Confirm PGS Catalog's overall EBI terms of use, and the specific license of each
      individual scoring file used in Track B (some carry their own CC/non-commercial terms)
- [ ] Confirm GCTA, PRSice-2, LDpred2/`bigsnpr`, PRS-CSx, and `pgscatalog-utils`/`pgsc_calc`
      license terms (check each repository's LICENSE file directly)
- [ ] Cite Momin et al. 2026, Ruan et al. 2022 (PRS-CSx), and the PGS Catalog's own
      publications as the methodological foundation this project builds on
- [x] Choose and apply an explicit open-source license for original code — MIT (see `LICENSE`)

## Data & Tool Sources

| Resource | License (to confirm) | Notes |
|---|---|---|
| 1000 Genomes Phase 3 | Open, consent-based public resource | Confirm exact IGSR terms before publication |
| PGS Catalog scoring files | Varies per score | EBI general terms + per-score license; check each score used in Track B individually |
| GCTA | See upstream repo | https://cnsgenomics.com/software/gcta/ |
| PRSice-2 | See upstream repo | https://github.com/choishingwan/PRSice |
| `bigsnpr`/LDpred2 | See upstream repo (CRAN) | https://privefl.github.io/bigsnpr/ |
| PRS-CSx | See upstream repo | https://github.com/getian107/PRScsx |
| `pgscatalog-utils`/`pgsc_calc` | See upstream repo | https://github.com/PGScatalog/pgscatalog_utils |
