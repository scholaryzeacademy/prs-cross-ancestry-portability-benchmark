# LICENSES.md

Tracks licensing status per `docs/BUILD_PLAN.md` §11. Original code in this repository is
licensed under MIT (see `LICENSE`) — a portfolio-appropriate default; revisit if the project
owner prefers Apache-2.0. Third-party data and tools each carry their own license, checked
individually below rather than assumed uniform.

## Checklist (BUILD_PLAN.md §11)

- [x] Confirm 1000 Genomes/IGSR data usage terms — see "1000 Genomes / IGSR data" below;
      **not** a blanket public-domain grant, contrary to this file's earlier draft note
- [x] Confirm PGS Catalog's overall EBI terms of use, and the specific license of each
      individual scoring file used in Track B (some carry their own CC/non-commercial terms) —
      see "PGS Catalog" below
- [x] Confirm GCTA, PRSice-2, LDpred2/`bigsnpr`, PRS-CSx, and `pgscatalog-utils`/`pgsc_calc`
      license terms (check each repository's LICENSE file directly) — see table below
- [x] Cite Momin et al. 2026, Ruan et al. 2022 (PRS-CSx), and the PGS Catalog's own
      publications as the methodological foundation this project builds on — done, see
      `reports/final_writeup/TECHNICAL_WRITEUP.md` References
- [x] Choose and apply an explicit open-source license for original code — MIT (see `LICENSE`)

## Data Sources

### 1000 Genomes / IGSR data

**Not public domain, and not covered by a single blanket license.** Per IGSR's own data
disclaimer (internationalgenome.org/IGSR_disclaimer, checked 2026-08-02): "data made available
by IGSR comes from many different owners and consequently restrictions on different pieces of
data within IGSR and rights claimed on pieces of data vary." IGSR states the Phase 3 data itself
is available without embargo following final publication (the standard cited paper is the 1000
Genomes Project Consortium's 2015 *Nature* paper, already in this project's References), but
explicitly puts the burden on users: "It remains the responsibility of users to ensure that
their exploitation of the data does not infringe any of the rights of third parties, including
the data owners." This project's use — public, published, de-identified genotype data, for
non-commercial portfolio/methodology-demonstration research, with the source paper cited — sits
squarely within the intended reuse case IGSR describes, but this is a judgment call, not a
license grant in the OSI sense. Anyone repurposing this project's *code* against 1000 Genomes
data (or any other IGSR dataset) should re-check IGSR's disclaimer for that specific dataset
rather than assume this project's use case transfers.

### PGS Catalog / EBI

EBI's general Terms of Use (ebi.ac.uk/about/terms-of-use, checked 2026-08-02): "EMBL-EBI itself
places no additional restrictions on the use or redistribution of the data available via its
Data Resources and Tools other than those provided by the original data owners," and EBI
"expects attribution ... in accordance with good scientific practice." The PGS Catalog's own
about page adds: "Individual PGS obtained from the database should also be cited appropriately,
and used in accordance with any licensing restrictions set by the authors" — i.e. EBI does not
impose a uniform license across all scoring files; some may carry author-specific CC or
non-commercial terms layered on top of EBI's general terms.

**This project's specific score, PGS000297** ("GRS3290_Height," Xie et al. 2020): no
score-specific license override is listed on its PGS Catalog page beyond the general terms
above. This project cites the original publication (Xie et al. 2020, already in this project's
References) and uses the score descriptively/illustratively (Track B, never as a validated
accuracy claim — see `CLAUDE.md`'s Non-negotiable framing section), consistent with the general
terms. Anyone using a *different* PGS Catalog score should re-check that score's own page for an
author-specific license before assuming these same terms apply.

## Tool Licenses

All confirmed directly from the license file bundled with each tool (where present locally in
`tools/`, gitignored) or the tool's own repository, checked 2026-08-02.

| Tool | License | Source checked |
|---|---|---|
| GCTA | MIT | `tools/gcta/gcta-1.95.3-linux-x86_64/MIT_License.txt` (local) |
| PRSice-2 | GPL-3.0 | github.com/choishingwan/PRSice (repo license badge) |
| `bigsnpr` (LDpred2-auto) | GPL-3 | Installed package `DESCRIPTION` (`/home/jonaid/R/library/bigsnpr`) |
| `bigsparser` (bigsnpr dependency, used directly by this project's LDpred2 wrapper) | GPL-3 | Installed package `DESCRIPTION` |
| PRS-CSx | MIT | `tools/prscsx/PRScsx-master/LICENSE` (local; Copyright Tian Ge) |
| `pgscatalog-utils`/`pgscatalog.core`/`pgscatalog.calc`/`pgscatalog.match` | Apache-2.0 | Installed package `dist-info/licenses/LICENSE` (pip) |
| `pgsc_calc` | Apache-2.0 | github.com/PGScatalog/pgsc_calc (repo LICENSE) |
| plink2 (used throughout Stages 1/3/5/6 and Track B for QC/scoring) | GPLv3 core (`COPYING`), LGPLv3 for the `pgenlibr` library component (`COPYING.LESSER`) | github.com/chrchang/plink-ng `2.0/COPYING`, `2.0/COPYING.LESSER` |
| plink2's bundled Intel MKL component | Intel Simplified Software License (binary-only redistribution, no modification/reverse-engineering) | `tools/plink2/intel-simplified-software-license.txt` (local) — governs only the bundled MKL binary, not plink2 itself |

**Compliance note:** this project invokes every tool above as an external binary/CLI/R-package
dependency (subprocess calls or library imports), never by copying or modifying their source
into this repository — so none of GCTA/PRSice-2/bigsnpr/PRS-CSx/pgscatalog-utils/plink2's
licenses (including the GPL/LGPL ones) impose any licensing obligation on this project's own
MIT-licensed code. `tools/` and `data/` are gitignored (per `CLAUDE.md`), so none of these
third-party binaries or their license files are redistributed via this repository — anyone
reproducing this project fetches them independently via `scripts/verify_tools.py --install`.
