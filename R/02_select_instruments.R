#' 02_select_instruments.R
#'
#' Purpose: select genetic instruments for each exposure and mediator per
#'   ANALYSIS_PLAN.md §6 (tasks T2.1-T2.9).
#'   1. genome-wide significance, p < 5e-8
#'   2. independence: local PLINK clumping, r2 < 0.001 within 10,000 kb,
#'      1000 Genomes EUR reference (data/reference/1kg_v3/EUR)
#'   3. per-SNP F = (beta/se)^2 > 10
#'   5. palindromic SNPs (A/T, C/G) with MAF > 0.42 excluded
#'   Criteria 4 (presence in outcome GWAS / proxies) and 6 (MR-Steiger) need
#'   outcome data and are applied at harmonization (03_harmonize.R; T3.6).
#'
#' Inputs:  data/raw/* summary statistics (see SOURCES.md / MANIFEST.tsv);
#'          data/reference/1kg_v3/EUR.{bed,bim,fam}; PLINK via plinkbinr.
#' Outputs: data/instruments/{dataset}_instruments.tsv (one per exposure/mediator)
#'          data/instruments/INSTRUMENT_SUMMARY.tsv
#'          data/instruments/EXPOSURE_DIRECTIONS.tsv
#'          data/instruments/SANITY_CHECKS.tsv
#'          results/logs/{date}_02_select_instruments.log (+ sessionInfo)
#' Depends on: data.table, ieugwasr, plinkbinr, glue, here
#'
#' Sign convention: all SES exposures oriented so that higher = higher SES.
#' TDI (higher = more deprived) is multiplied by -1 before selection outputs.

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
ensure_pkgs_loaded(c("data.table", "ieugwasr", "plinkbinr", "glue"))

P_THRESH  <- 5e-8
CLUMP_KB  <- 10000
CLUMP_R2  <- 0.001
F_MIN     <- 10
PAL_MAF   <- 0.42
BFILE_EUR <- here::here("data", "reference", "1kg_v3", "EUR")
PLINK_BIN <- plinkbinr::get_plink_exe()
fs::dir_create(DIR_INSTRUMENTS)

log_con <- log_run("02_select_instruments")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)  # sessionInfo() is written by close_run_log() at the end
logl("PLINK: {PLINK_BIN}; reference: {BFILE_EUR}; p<{P_THRESH}, r2<{CLUMP_R2}, {CLUMP_KB} kb, F>{F_MIN}, palindrome MAF>{PAL_MAF} excluded")

# ---- Reference panel (EUR .bim) ----------------------------------------------
# Used to (a) map TDI chr:pos:alleles to rsIDs, (b) supply EUR allele
# frequencies where a file has none, and (c) know which rsIDs PLINK can clump.
bim <- fread(paste0(BFILE_EUR, ".bim"), header = FALSE,
             col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"),
             colClasses = c("character", "character", "numeric", "integer", "character", "character"))
bim <- bim[grepl("^rs[0-9]+$", rsid)]
logl("EUR bim: {nrow(bim)} rsID variants")

# ---- Dataset registry and readers (shared) ----------------------------------
source(here::here("R", "utils", "sumstats_registry.R"))

# ---- Main loop ---------------------------------------------------------------

summary_rows <- list(); counts <- list()
# Optional subset (added 2026-10-06): SELECT_ONLY="key1,key2" selects instruments for those registry keys only and
# updates just their rows in INSTRUMENT_SUMMARY.tsv; other instrument files are left untouched.
SELECT_ONLY <- Sys.getenv("SELECT_ONLY")
run_keys <- if (nzchar(SELECT_ONLY)) intersect(strsplit(SELECT_ONLY, ",")[[1]], names(registry)) else names(registry)
if (nzchar(SELECT_ONLY)) logl("SELECT_ONLY: {paste(run_keys, collapse = ', ')}")
for (key in run_keys) {
  cfg <- registry[[key]]
  t0 <- Sys.time()
  logl("=== {key}: {cfg$label}")
  d <- if (identical(cfg$special, "panukb")) read_panukb(cfg, p_max = P_THRESH, bim = bim) else read_standard(cfg)
  if (!identical(cfg$special, "panukb")) d <- d[!is.na(pval) & pval < P_THRESH]
  d <- std(d, cfg)
  n_gw <- nrow(d)
  d <- d[!is.na(beta) & !is.na(se) & se > 0 & grepl("^rs[0-9]+$", rsid)]
  d <- unique(d, by = "rsid")
  in_ref <- d$rsid %in% bim$rsid
  n_gw_rs <- nrow(d); n_not_in_ref <- sum(!in_ref)
  d <- d[in_ref]
  # clump
  cl <- ieugwasr::ld_clump(data.frame(rsid = d$rsid, pval = d$pval, id = key),
                           clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, clump_p = P_THRESH,
                           bfile = BFILE_EUR, plink_bin = PLINK_BIN)
  d <- d[rsid %in% cl$rsid]
  n_clumped <- nrow(d)
  # fill missing EAF from EUR reference, aligned to the effect allele
  if (anyNA(d$eaf)) {
    fr <- eur_freq(d[is.na(eaf), rsid])
    d <- merge(d, fr, by = "rsid", all.x = TRUE)
    d[is.na(eaf) & ea == ref_a1, eaf := ref_maf]
    d[is.na(eaf) & ea == ref_a2, eaf := 1 - ref_maf]
    d[, c("ref_a1", "ref_a2", "ref_maf") := NULL]
  }
  d[, F_stat := (beta / se)^2]
  n_weak <- sum(d$F_stat <= F_MIN)
  d <- d[F_stat > F_MIN]
  d[, maf := pmin(eaf, 1 - eaf)]
  pal <- is_palindromic(d$ea, d$oa)
  n_pal_amb <- sum(pal & (is.na(d$maf) | d$maf > PAL_MAF))
  d <- d[!(pal & (is.na(maf) | maf > PAL_MAF))]
  d[, palindromic := is_palindromic(ea, oa)]
  # variance explained, scale-free: r2_j = z^2 / (z^2 + N - 2)
  d[, r2 := F_stat / (F_stat + n - 2)]
  setorder(d, chr, pos)
  out <- d[, .(rsid, chr, pos, effect_allele = ea, other_allele = oa, eaf, beta, se, pval, n, F_stat, r2, palindromic)]
  fwrite(out, file.path(DIR_INSTRUMENTS, glue("{key}_instruments.tsv")), sep = "\t")
  K <- nrow(out); R2 <- sum(out$r2); Nm <- median(out$n, na.rm = TRUE)
  summary_rows[[key]] <- data.table(
    dataset = key, role = cfg$role, label = cfg$label, ancestry_ref = "EUR (1000G)",
    n_gw_rows = n_gw, n_gw_rsid = n_gw_rs, n_not_in_ref = n_not_in_ref, n_clumped = n_clumped,
    n_weak_removed = n_weak, n_palindromic_ambiguous_removed = n_pal_amb, n_snps = K,
    mean_F = round(mean(out$F_stat), 1), median_F = round(median(out$F_stat), 1),
    min_F = round(min(out$F_stat), 1), R2 = signif(R2, 4),
    overall_F = round(((Nm - K - 1) / K) * (R2 / (1 - R2)), 1),
    median_n = round(Nm),
    n_source = if ("n" %in% names(cfg$cols)) "per-SNP column" else if (!is.null(cfg$n_source)) cfg$n_source else "documented total (n_default)",
    units = cfg$units)
  logl("  GW-sig rows {n_gw}; rsID {n_gw_rs}; not in ref {n_not_in_ref}; clumped {n_clumped}; final {K}; mean F {round(mean(out$F_stat),1)}; R2 {signif(R2,3)} ({round(as.numeric(difftime(Sys.time(), t0, units='mins')),1)} min)")
  rm(d, cl, out); gc()
}

summary_dt <- rbindlist(summary_rows)
if (nzchar(SELECT_ONLY) && file.exists(file.path(DIR_INSTRUMENTS, "INSTRUMENT_SUMMARY.tsv"))) {
  old_dt <- fread(file.path(DIR_INSTRUMENTS, "INSTRUMENT_SUMMARY.tsv"))
  summary_dt <- rbind(old_dt[!dataset %in% run_keys], summary_dt, fill = TRUE)
}
fwrite(summary_dt, file.path(DIR_INSTRUMENTS, "INSTRUMENT_SUMMARY.tsv"), sep = "\t")

directions <- rbindlist(lapply(names(registry), function(k)
  data.table(dataset = k, role = registry[[k]]$role, sign_applied = registry[[k]]$sign,
             orientation = if (registry[[k]]$sign < 0) "flipped: higher = less deprived (higher SES)" else "as published",
             units = registry[[k]]$units)))
fwrite(directions, file.path(DIR_INSTRUMENTS, "EXPOSURE_DIRECTIONS.tsv"), sep = "\t")

# ---- T2.9 sanity checks: known EA loci (GRCh37) --------------------------------
# ANKRD55 was listed in TASKS.md T2.9 as an expected EA locus, but EA4 has no genome-wide
# significant SNP within 1 Mb of it (min p = 2e-5); it is kept here as a documented negative check.
loci <- data.table(locus = c("ANKRD55", "MAPT"), chr = c("5", "17"),
                   start = c(55395000L, 43971000L), end = c(55530000L, 44106000L))
checks <- rbindlist(lapply(c("ea_okbay2022", "ea_lee2018"), function(k) {
  ins <- fread(file.path(DIR_INSTRUMENTS, glue("{k}_instruments.tsv")), colClasses = list(character = "chr"))
  rbindlist(lapply(seq_len(nrow(loci)), function(i) {
    l <- loci[i]; s <- ins[chr == l$chr]
    dist <- if (nrow(s)) min(pmax(0, pmax(l$start - s$pos, s$pos - l$end))) else NA_real_
    data.table(dataset = k, locus = l$locus, nearest_instrument_kb = round(dist / 1000, 1),
               within_500kb = !is.na(dist) & dist <= 5e5)
  }))
}))
fwrite(checks, file.path(DIR_INSTRUMENTS, "SANITY_CHECKS.tsv"), sep = "\t")
logl("Sanity checks:\n{paste(capture.output(print(checks)), collapse = '\n')}")
logl("Done.")
close_run_log(log_con)
