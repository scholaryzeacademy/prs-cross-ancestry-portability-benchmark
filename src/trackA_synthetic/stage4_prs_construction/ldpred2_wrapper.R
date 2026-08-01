#!/usr/bin/env Rscript
#
# ldpred2_wrapper.R — Stage 4, LDpred2-auto (BUILD_PLAN.md §6 Stage 4).
#
# Runs LDpred2-auto (bigsnpr::snp_ldpred2_auto) on Stage 3's EUR discovery-only GWAS, using the
# EUR discovery samples themselves as the LD reference (never the held-out EUR subset or the
# other four ancestries — construction must stay independent of any Stage 5 evaluation target).
# ldpred2-auto is used specifically because it estimates its own hyperparameters (h2, causal
# fraction p) from the summary statistics and LD reference alone via Gibbs sampling — no target
# phenotype needed, keeping construction cleanly separated from evaluation, same as the
# PRSice-2 wrapper's design.
#
# Known simplification: this project has no genetic map, so snp_cor's LD window is defined in
# physical base pairs (250kb radius) rather than centimorgans. This is a documented, common
# approximation when a genetic map isn't available — noted in METHODS.md.
#
# LD (genotype-only) is computed once and cached, then reused for both scenarios' df_beta — it
# doesn't depend on phenotype/scenario, and this LD step is the expensive part.
#
# Usage:
#   Rscript ldpred2_wrapper.R --scenario 1
#   Rscript ldpred2_wrapper.R --scenario all

Sys.setenv(R_LIBS_USER = file.path(Sys.getenv("HOME"), "R", "library"))  # so parallel workers (ncores>1) find bigsnpr's deps (e.g. RhpcBLASctl)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
  library(bigsnpr)
  library(bigsparser)
  library(data.table)
  library(Matrix)
})

repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), "..", "..", ".."))
shared_dir <- file.path(repo_root, "data", "processed", "simulated_phenotypes", "_shared")
bed_prefix <- file.path(shared_dir, "1000g_qc_biallelic")
eur_discovery_keep <- file.path(repo_root, "data", "processed", "gwas", "_shared", "keep_EUR_discovery.txt")
gwas_root <- file.path(repo_root, "data", "processed", "gwas")
prs_models_root <- file.path(repo_root, "data", "processed", "prs_models")
ldpred2_shared_dir <- file.path(prs_models_root, "_shared", "ldpred2")
plink2_bin <- if (nzchar(Sys.which("plink2"))) "plink2" else file.path(repo_root, "tools", "plink2", "plink2")

scenarios <- list(`1` = "scenario1_equal_effects", `2` = "scenario2_ancestry_varying_effects")

ld_window_bp <- 100000  # radius in bp; see module docstring re: no genetic map available
n_cores <- 4
# LD computation over the full chr21+chr22 panel (131K SNPs) didn't complete in 28+ minutes on
# this shared, heavily-loaded host (load average 150+ on 32 cores during testing). Restricted to
# chr22 alone for this smoke test, consistent with the same scope reduction CI already uses
# (.github/workflows/ci.yml, BUILD_PLAN.md §5 CI row) — documented in METHODS.md, not silent.
ld_chroms <- "22"

prepare_shared_ld <- function() {
  # Computed once per script invocation and reused in-memory for both scenarios (main() always
  # processes all requested scenarios within a single R session) since LD depends only on
  # genotypes, not on which scenario's phenotype/betas are being scored.
  dir.create(ldpred2_shared_dir, recursive = TRUE, showWarnings = FALSE)

  eur_bed_prefix <- file.path(ldpred2_shared_dir, "eur_discovery")
  system2(plink2_bin, c(
    "--bfile", bed_prefix, "--keep", eur_discovery_keep,
    "--make-bed", "--out", eur_bed_prefix
  ), stdout = file.path(ldpred2_shared_dir, "plink2_subset.log"), stderr = file.path(ldpred2_shared_dir, "plink2_subset.log"))

  rds_path <- snp_readBed(paste0(eur_bed_prefix, ".bed"), backingfile = file.path(ldpred2_shared_dir, "eur_discovery_bigsnp"))
  obj.bigSNP <- snp_attach(rds_path)
  G <- obj.bigSNP$genotypes
  POS <- obj.bigSNP$map$physical.pos
  map <- obj.bigSNP$map
  names(map) <- c("chr", "rsid", "genetic.dist", "pos", "a1", "a0")

  # Any scenario's sumstats covers the same 154,384 variant IDs/alleles (only BETA/SE differ
  # between scenarios), so matching against Scenario 1 here defines the canonical
  # ambiguous-SNP-filtered SNP set/order that both scenarios' df_beta will be reindexed to.
  sumstats <- fread(file.path(gwas_root, scenarios[["1"]], "eur_discovery_gwas.SYNTH_PHENO.glm.linear"))
  names(sumstats)[1] <- "CHROM"
  setnames(sumstats,
           c("CHROM", "POS", "A1", "OMITTED", "BETA", "SE", "OBS_CT"),
           c("chr", "pos", "a1", "a0", "beta", "beta_se", "n_eff"))
  info_snp <- snp_match(sumstats, map)
  info_snp <- info_snp[info_snp$chr %in% ld_chroms, ]
  cat(sprintf("canonical SNP set: %d matched (of %d in map), restricted to chr %s\n",
              nrow(info_snp), nrow(map), paste(ld_chroms, collapse = ",")))

  chrs <- sort(unique(info_snp$chr))
  corr_blocks <- vector("list", length(chrs))
  canonical_ids <- vector("list", length(chrs))
  for (i in seq_along(chrs)) {
    ind.chr <- which(info_snp$chr == chrs[i])
    ind.col <- info_snp$`_NUM_ID_`[ind.chr]
    t0 <- Sys.time()
    corr_blocks[[i]] <- snp_cor(G, ind.col = ind.col, size = ld_window_bp, infos.pos = POS[ind.col], ncores = n_cores)
    cat(sprintf("chr%s: %d SNPs, snp_cor took %.1fs\n", chrs[i], length(ind.col), as.numeric(Sys.time() - t0, units = "secs")))
    canonical_ids[[i]] <- info_snp$rsid[ind.chr]
  }
  corr_dense <- bdiag(corr_blocks)
  corr <- as_SFBM(corr_dense, backingfile = file.path(ldpred2_shared_dir, "corr_sfbm"), compact = TRUE)
  ld <- Matrix::colSums(corr_dense^2)
  list(corr = corr, ld = ld, canonical_rsids = unlist(canonical_ids), map_rds = rds_path)
}

run_ldpred2_for_scenario <- function(scenario_name, ld_meta) {
  gwas_file <- file.path(gwas_root, scenario_name, "eur_discovery_gwas.SYNTH_PHENO.glm.linear")
  if (!file.exists(gwas_file)) {
    stop(sprintf("%s not found - run src/trackA_synthetic/stage3_eur_gwas/run_gwas.py first.", gwas_file))
  }
  out_dir <- file.path(prs_models_root, scenario_name, "ldpred2")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  obj.bigSNP <- snp_attach(ld_meta$map_rds)
  map <- obj.bigSNP$map
  names(map) <- c("chr", "rsid", "genetic.dist", "pos", "a1", "a0")

  sumstats <- fread(gwas_file)
  names(sumstats)[1] <- "CHROM"
  setnames(sumstats,
           c("CHROM", "POS", "A1", "OMITTED", "BETA", "SE", "OBS_CT"),
           c("chr", "pos", "a1", "a0", "beta", "beta_se", "n_eff"))
  info_snp <- snp_match(sumstats, map)

  # Reindex to the exact canonical order the shared corr matrix's columns were built in.
  df_beta <- info_snp[match(ld_meta$canonical_rsids, info_snp$rsid), ]
  stopifnot(!anyNA(df_beta$rsid))
  cat(sprintf("[%s] %d SNPs aligned to shared LD matrix\n", scenario_name, nrow(df_beta)))

  chi2 <- (df_beta$beta / df_beta$beta_se)^2
  ldsc <- snp_ldsc(ld_meta$ld, ld_size = nrow(df_beta), chi2 = chi2, sample_size = df_beta$n_eff, blocks = NULL)
  h2_init <- max(ldsc[["h2"]], 0.05)
  cat(sprintf("[%s] LDSC h2 estimate (init for auto model): %.4f\n", scenario_name, h2_init))

  # burn_in/num_iter/chain count reduced from bigsnpr's usual defaults (500/200/10 chains) —
  # this host's contention meant even the chr22-only LD step alone took ~14 minutes; a smaller
  # MCMC budget keeps the whole run tractable while still giving a real (if noisier) estimate.
  auto <- snp_ldpred2_auto(
    ld_meta$corr, df_beta, h2_init = h2_init, vec_p_init = seq_log(1e-4, 0.5, length.out = 4),
    burn_in = 200, num_iter = 100, allow_jump_sign = FALSE, use_MLE = FALSE, ncores = n_cores,
  )
  beta_auto <- rowMeans(sapply(auto, function(a) a$beta_est))

  scoring <- data.table(ID = df_beta$rsid, A1 = df_beta$a1, BETA = beta_auto)
  scoring <- scoring[BETA != 0]
  dest <- file.path(out_dir, "ldpred2_auto.scoring.tsv")
  fwrite(scoring, dest, sep = "\t")
  cat(sprintf("[%s] LDpred2-auto scoring file: %d non-zero SNPs -> %s\n", scenario_name, nrow(scoring), dest))
  dest
}

main <- function() {
  option_list <- list(
    make_option("--scenario", type = "character", default = "all")
  )
  args <- parse_args(OptionParser(option_list = option_list))

  ld_meta <- prepare_shared_ld()

  scenario_ids <- if (args$scenario == "all") c("1", "2") else args$scenario
  for (sid in scenario_ids) {
    run_ldpred2_for_scenario(scenarios[[sid]], ld_meta)
  }
}

main()
