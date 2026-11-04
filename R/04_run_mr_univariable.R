#' 04_run_mr_univariable.R
#'
#' Purpose: univariable MR (plan v1.0 §7, §8, §11, §12) for every pre-specified exposure x outcome cell: the 16-cell
#'   primary H1 grid (EA4, Kweon income, TDI, CAMSIS x CAD Aragam EUR, PAD MVP EUR, HF HERMES 2.0, any stroke
#'   GIGASTROKE multi), comparators (§7.2), exploratory cells (MVP AFR §7.3; van Zuydam rounded data §4.2) and the gSES
#'   secondary analysis (§7.4). Every cell is reported regardless of its result. ISEI is not analysed univariably (it
#'   enters only the comparable MVMR model, §3).
#'   Main estimator: IVW, multiplicative random effects. Battery (§8): MR-Egger, weighted median, weighted mode,
#'   MR-PRESSO, Cochran's Q / I2, leave-one-out, single-SNP, MR-Steiger (docs/analysis_specs/steiger_prevalence_spec.md),
#'   MR-RAPS. Sensitivity variants (IVW): proxies, MHC excluded, residual-LD pruned, MVP r2 >= 0.8, Steiger-filtered,
#'   van Zuydam conservative SE. Implementation details: results/logs/DEVIATIONS.md row "Phase 4 implementation".
#' Modes (Rscript R/04_run_mr_univariable.R <mode>):
#'   prep        (local only; needs data/raw) write per-SNP and study-total case/control counts to data/cell_inputs/
#'   cell K      run grid cell K (1..96; UGER array task) -> results/cells/univariable/{set}__{outcome}.rds (skipped if done)
#'   merge       combine the 96 cell files into the tables below (stops unless all 96 are present)
#'   check       re-run the reference cells in a temp dir and compare with results/reference/ (stops on any difference)
#'   local       (default) prep if needed, run all missing cells in parallel (MR_CORES), then merge
#'   Every cell sets its own seeds (PROJECT_SEED + cell index), so results do not depend on where or in which order
#'   cells run.
#' Outputs: results/tables/mr_univariable_{estimates,diagnostics,variants,loo,singlesnp,steiger_snp}.tsv,
#'   results/tables/mr_primary_grid.tsv
#' Depends on: data.table, MendelianRandomization, MRPRESSO, mr.raps, TwoSampleMR, psych, parallel

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args)) args[1] else "local"
stopifnot(MODE %in% c("prep", "cell", "merge", "check", "local"))
CELL_K <- if (MODE == "cell") as.integer(if (length(args) >= 2) args[2] else Sys.getenv("SGE_TASK_ID")) else NA_integer_
log_con <- log_run(if (MODE == "cell") sprintf("04_run_mr_univariable_cell%02d", CELL_K) else paste0("04_run_mr_univariable_", MODE))
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
source(here::here("R", "utils", "mr_methods.R"))
H <- file.path(DIR_DATA_HARM, "large"); RAW <- DIR_DATA_RAW
N_CORES <- as.integer(Sys.getenv("MR_CORES", "14")); N_PRESSO <- 10000L; N_BOOT <- 10000L
DIR_CELLS <- file.path(DIR_RESULTS, "cells", "univariable"); DIR_CIN <- here::here("data", "cell_inputs")
DIR_REF <- file.path(DIR_RESULTS, "reference")

sets <- c(EA4 = "ea_okbay2022", Kweon = "income_kweon2025", TDI = "tdi_panukb", CAMSIS = "occupation_akimova2024",
          EA3 = "ea_lee2018", Hill = "income_hill2019", TDI_EUR = "tdi_panukb_eur", gSES = "ses_gses_xia2025")
outcomes <- c("cad_aragam_eur", "cad_aragam_multi", "cad_mvp_phe411_eur", "pad_mvp_eur", "pad_finngen_athscle",
              "pad_vanzuydam", "pad_mvp_afr", "hf_hermes2", "hf_mvp_phe428_eur", "stroke_gigastroke_multi",
              "stroke_gigastroke_eur", "stroke_mvp_phe433_21_eur")
primary_out <- c("cad_aragam_eur", "pad_mvp_eur", "hf_hermes2", "stroke_gigastroke_multi")
outcome_class <- c(cad = "CAD", pad = "PAD", hf = "HF", stroke = "Any stroke")
prev <- c(CAD = 0.084, PAD = 0.071, HF = 0.024, `Any stroke` = 0.029)   # steiger_prevalence_spec.md (fixed 2026-10-08)
role_of <- function(s, o) {
  if (s %in% c("EA4", "Kweon", "TDI", "CAMSIS") && o %in% primary_out) return("primary")
  if (s == "gSES") return("secondary (gSES)")
  if (o == "pad_mvp_afr") return("exploratory (transportability)")
  if (o == "pad_vanzuydam") return("exploratory (rounded data)")
  "comparator / sensitivity"
}

# ---- outcome case/control counts (spec rule 2: per SNP where the file gives both, else study totals) -------------
prep_n <- function() {
man <- fread(PATH_MANIFEST)
tot <- function(ds) man[dataset == ds, .(ncase = as.numeric(cases), ncontrol = as.numeric(controls))]
totals <- list(cad_aragam_eur = tot("cad_aragam2022_eur"), cad_aragam_multi = tot("cad_aragam2022_multi"),
               pad_finngen_athscle = tot("pad_finngen_r12_athscle"), stroke_gigastroke_multi = tot("stroke_gigastroke_any_multi"),
               stroke_gigastroke_eur = tot("stroke_gigastroke_any_eur"), pad_mvp_eur = tot("pad_mvp_phe440_2_eur"),
               pad_mvp_afr = tot("pad_mvp_phe440_2_afr"), cad_mvp_phe411_eur = tot("cad_mvp_phe411_eur"),
               hf_hermes2 = tot("hf_hermes2_2025"), hf_mvp_phe428_eur = tot("hf_mvp_phe428_eur"),
               stroke_mvp_phe433_21_eur = tot("stroke_mvp_phe433_21_eur"), pad_vanzuydam = tot("pad_vanzuydam2021"))
stopifnot(all(sapply(totals, nrow) == 1), !anyNA(unlist(totals)))
snps_all <- unique(unlist(lapply(names(sets), function(s) fread(file.path(DIR_INSTRUMENTS, paste0(sets[[s]], "_instruments.tsv")))$rsid)))
inst_pos <- unique(rbindlist(lapply(sets, function(k) fread(file.path(DIR_INSTRUMENTS, paste0(k, "_instruments.tsv")),
                                                              select = c("rsid", "chr", "pos"), colClasses = list(character = "chr")))), by = "rsid")
per_snp_n <- function(o) {
  pos_join <- function(d) { d[, chr := sub("^chr", "", as.character(chr))]; merge(inst_pos, d, by = c("chr", "pos"))[, .(rsid, a = toupper(a), b = toupper(b), ncase, ncontrol)] }
  if (o %in% c("cad_mvp_phe411_eur", "pad_mvp_eur", "pad_mvp_afr")) {
    f <- c(cad_mvp_phe411_eur = "cad_mvp_phe411_EUR.txt.gz", pad_mvp_eur = "pad_mvp_phe440_2_EUR.txt.gz", pad_mvp_afr = "pad_mvp_phe440_2_AFR.txt.gz")[[o]]
    d <- fread(file.path(RAW, f), select = c("SNP_ID", "ref", "alt", "num_cases", "num_controls"))[SNP_ID %in% snps_all]
    return(d[, .(rsid = SNP_ID, a = toupper(alt), b = toupper(ref), ncase = as.numeric(num_cases), ncontrol = as.numeric(num_controls))])
  }
  if (o %in% c("hf_mvp_phe428_eur", "stroke_mvp_phe433_21_eur")) {
    f <- c(hf_mvp_phe428_eur = "hf_mvp_phe428_EUR.gwascat.tsv.gz", stroke_mvp_phe433_21_eur = "stroke_mvp_phe433_21_EUR.gwascat.tsv.gz")[[o]]
    d <- fread(file.path(RAW, f), select = c("rsid", "effect_allele", "other_allele", "num_cases", "num_controls"))[rsid %in% snps_all]
    return(d[, .(rsid, a = toupper(effect_allele), b = toupper(other_allele), ncase = as.numeric(num_cases), ncontrol = as.numeric(num_controls))])
  }
  if (o == "hf_hermes2") {
    d <- fread(file.path(RAW, "hf_hermes2_2025_hfall.tsv.gz"), select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "n_case", "N_total"))
    d <- d[, .(chr = chromosome, pos = as.integer(base_pair_location), a = effect_allele, b = other_allele, ncase = as.numeric(n_case), ncontrol = as.numeric(N_total) - as.numeric(n_case))]
    return(pos_join(d))
  }
  if (o == "pad_vanzuydam") {
    d <- fread(file.path(RAW, "pad_vanzuydam2021", "final_meta", "pad_primary.txt.gz"), select = c("chr", "bp", "or_allele", "nonor_allele", "ncases", "ncontrols"))
    d <- d[, .(chr = as.character(chr), pos = as.integer(bp), a = or_allele, b = nonor_allele, ncase = as.numeric(ncases), ncontrol = as.numeric(ncontrols))]
    return(pos_join(d))
  }
  NULL  # Aragam (controls not given per SNP), GIGASTROKE, FinnGen: study totals
}
nmap <- lapply(setNames(nm = outcomes), per_snp_n)
N <- rbindlist(lapply(names(nmap), function(o) {
  m <- nmap[[o]]; if (is.null(m)) return(NULL)
  m <- m[!is.na(ncase) & !is.na(ncontrol)]
  m <- m[, .SD[1], by = .(rsid, k = paste(pmin(a, b), pmax(a, b)))]
  m[, .(outcome = o, rsid, k, ncase, ncontrol)]
}))
fs::dir_create(DIR_CIN)
fwrite(N, file.path(DIR_CIN, "outcome_n_per_snp.tsv.gz"), sep = "\t")
fwrite(rbindlist(lapply(names(totals), function(o) cbind(outcome = o, totals[[o]]))), file.path(DIR_CIN, "outcome_n_totals.tsv"), sep = "\t")
logl("prep: {nrow(N)} per-SNP count rows for {uniqueN(N$outcome)} outcomes; totals for {length(totals)} outcomes")
}
load_n <- function() {
  N <- fread(file.path(DIR_CIN, "outcome_n_per_snp.tsv.gz"))
  Tt <- fread(file.path(DIR_CIN, "outcome_n_totals.tsv"))
  stopifnot(setequal(Tt$outcome, outcomes), !anyNA(Tt))
  list(N = N, totals = split(Tt[, .(ncase, ncontrol)], Tt$outcome))
}
attach_n <- function(d, o) {
  tt <- NT$totals[[o]]
  d[, `:=`(ncase = tt$ncase, ncontrol = tt$ncontrol, n_source = "study total")]
  m <- NT$N[outcome == o]
  if (nrow(m)) {
    d[, k := paste(pmin(ea_y, oa_y), pmax(ea_y, oa_y))]
    d[m, on = .(rsid, k), `:=`(ncase = i.ncase, ncontrol = i.ncontrol, n_source = "per SNP")]
    d[, k := NULL]
  }
  d
}

# ---- one cell ---------------------------------------------------------------------------------------------------
run_cell <- function(s, o, seed) {
  set.seed(seed)
  cls <- outcome_class[[sub("_.*", "", o)]]
  h <- fread(file.path(H, "univariable", paste0(s, "__", o, ".tsv.gz")), colClasses = list(character = c("ea_x", "oa_x", "ea_y", "oa_y", "chr")))
  inst <- fread(file.path(DIR_INSTRUMENTS, paste0(sets[[s]], "_instruments.tsv")), select = c("rsid", "n"))
  d <- h[mr_keep == TRUE]; d <- merge(d, inst, by = "rsid", all.x = TRUE)
  stopifnot(!anyNA(d$n))
  bx <- d$beta_x; sx <- d$se_x; by <- d$beta_y_aligned; sy <- d$se_y; L <- nrow(d)
  meta <- data.table(set = s, outcome = o, outcome_class = cls, role = role_of(s, o))
  # IVW (MRE) and cross-check against MendelianRandomization
  iv <- ivw_mre(bx, by, sy)
  mi <- MendelianRandomization::mr_input(bx = bx, bxse = sx, by = by, byse = sy, snps = d$rsid)
  mr_iv <- MendelianRandomization::mr_ivw(mi, model = "random", distribution = "normal")
  if (abs(mr_iv@Estimate - iv$b) > 1e-10 || abs(mr_iv@StdError - iv$se) > 1e-10) stop("IVW cross-check failed: ", s, " ", o)
  est <- list(est_row("IVW (MRE)", iv$b, iv$se, L))
  eg <- MendelianRandomization::mr_egger(mi, distribution = "normal")
  est[[length(est) + 1]] <- est_row("MR-Egger slope", eg@Estimate, eg@StdError.Est, L, eg@CILower.Est, eg@CIUpper.Est, eg@Pvalue.Est)
  wm <- MendelianRandomization::mr_median(mi, weighting = "weighted", iterations = N_BOOT, seed = seed)
  est[[length(est) + 1]] <- est_row("Weighted median", wm@Estimate, wm@StdError, L, wm@CILower, wm@CIUpper, wm@Pvalue)
  mb <- MendelianRandomization::mr_mbe(mi, weighting = "weighted", stderror = "delta", phi = 1, iterations = N_BOOT, seed = seed)
  est[[length(est) + 1]] <- est_row("Weighted mode", mb@Estimate, mb@StdError, L, mb@CILower, mb@CIUpper, mb@Pvalue)
  rp <- tryCatch(mr.raps::mr.raps(data.frame(beta.exposure = bx, beta.outcome = by, se.exposure = sx, se.outcome = sy),
                                  diagnostics = FALSE, over.dispersion = TRUE, loss.function = "huber"), error = function(e) NULL)
  if (!is.null(rp)) est[[length(est) + 1]] <- est_row("MR-RAPS (overdispersed, Huber)", rp$beta.hat, rp$beta.se, L)
  # MR-PRESSO
  pr <- tryCatch(suppressWarnings(MRPRESSO::mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sy", SdExposure = "sx",
                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = data.frame(bx, by, sx, sy), NbDistribution = N_PRESSO,
                   SignifThreshold = 0.05, seed = seed)), error = function(e) e)
  presso_ok <- !inherits(pr, "error")
  if (presso_ok) {
    me <- pr$`Main MR results`
    raw_r <- me[me$`MR Analysis` == "Raw", ]; cor_r <- me[me$`MR Analysis` == "Outlier-corrected", ]
    est[[length(est) + 1]] <- est_row("MR-PRESSO raw", raw_r$`Causal Estimate`, raw_r$Sd, L, p = raw_r$`P-value`)
    if (nrow(cor_r) && !is.na(cor_r$`Causal Estimate`)) {
      n_out <- length(pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
      est[[length(est) + 1]] <- est_row("MR-PRESSO outlier-corrected", cor_r$`Causal Estimate`, cor_r$Sd, L - n_out, p = cor_r$`P-value`)
    }
    glob_p <- pr$`MR-PRESSO results`$`Global Test`$Pvalue
    outl <- if (!is.null(pr$`MR-PRESSO results`$`Distortion Test`)) pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices` else integer()
    if (!is.numeric(outl)) outl <- integer()  # MRPRESSO returns the string "No significant outliers" (fixed 2026-10-08)
    dist_p <- if (!is.null(pr$`MR-PRESSO results`$`Distortion Test`)) pr$`MR-PRESSO results`$`Distortion Test`$Pvalue else NA_real_
  }
  E <- cbind(meta, rbindlist(est))
  # leave-one-out, single-SNP
  loo <- cbind(meta, loo_ivw(d$rsid, bx, by, sy))
  loo[, change := b - iv$b]
  ss <- cbind(meta, data.table(rsid = d$rsid, b = by / bx, se = abs(sy / bx)))
  # Steiger at 0.5x, 1x, 2x prevalence
  d <- attach_n(d, o)
  eaf <- fifelse(is.na(d$eaf_y_aligned), d$eaf_x, d$eaf_y_aligned)
  st <- lapply(c(half = 0.5, primary = 1, double = 2), function(f) steiger_one(bx, sx, d$n, by, eaf, d$ncase, d$ncontrol, prev[[cls]] * f))
  sp <- st$primary
  st_snp <- cbind(meta, data.table(rsid = d$rsid, r2_exp = sp$r2_exp, r2_out = sp$r2_out, steiger_p = sp$p_snp,
                                   reverse_p05 = (!sp$correct) & sp$p_snp < 0.05, eaf_filled = is.na(d$eaf_y_aligned), n_source = d$n_source))
  n_rev <- sapply(st, function(z) sum((!z$correct) & z$p_snp < 0.05))
  D <- cbind(meta, data.table(nsnp = L, Q = iv$Q, Q_df = iv$Q_df, Q_p = iv$Q_p, I2 = iv$I2,
    egger_intercept = eg@Intercept, egger_intercept_se = eg@StdError.Int, egger_intercept_lo = eg@CILower.Int,
    egger_intercept_hi = eg@CIUpper.Int, egger_intercept_p = eg@Pvalue.Int,
    presso_global_p = if (presso_ok) glob_p else NA_real_, presso_n_outliers = if (presso_ok) length(outl) else NA_integer_,
    presso_distortion_p = if (presso_ok) dist_p else NA_real_, presso_error = if (presso_ok) NA_character_ else conditionMessage(pr),
    loo_max_abs_change = max(abs(loo$change)), loo_max_snp = loo$left_out[which.max(abs(loo$change))],
    steiger_n_correct = sum(sp$correct), steiger_share_correct = mean(sp$correct),
    steiger_overall_correct = sp$overall_correct, steiger_overall_p = sp$overall_p,
    steiger_overall_correct_half = st$half$overall_correct, steiger_overall_correct_double = st$double$overall_correct,
    steiger_n_reverse_p05_half = n_rev[["half"]], steiger_n_reverse_p05 = n_rev[["primary"]], steiger_n_reverse_p05_double = n_rev[["double"]],
    n_snp_n_per_snp = sum(d$n_source == "per SNP"), n_outcome_eaf_missing = sum(is.na(d$eaf_y_aligned)),
    n_or_equal_1 = if (o == "pad_vanzuydam") sum(round(exp(d$beta_y), 2) == 1) else NA_integer_,
    n_mvp_r2_lt_0.8 = if ("r2" %in% names(d) && any(!is.na(d$r2))) sum(d$r2 < 0.8, na.rm = TRUE) else NA_integer_))
  # sensitivity variants (IVW, MRE)
  V <- list()
  addv <- function(name, keep_rows = NULL, dd = NULL, sy_alt = NULL) {
    dd <- if (is.null(dd)) d[keep_rows] else dd
    syy <- if (is.null(sy_alt)) dd$se_y else sy_alt
    if (nrow(dd) < 2) return(invisible())
    r <- ivw_mre(dd$beta_x, dd$beta_y_aligned, syy)
    V[[length(V) + 1]] <<- cbind(meta, est_row(paste0("IVW: ", name), r$b, r$se, nrow(dd)))
  }
  hp <- fread(file.path(H, "univariable_proxy", paste0(s, "__", o, ".tsv.gz")), colClasses = list(character = c("ea_x", "oa_x", "ea_y", "oa_y", "chr")))
  if (sum(hp$mr_keep) != L) addv("with proxies", dd = hp[mr_keep == TRUE])
  if (any(d$mhc)) addv("MHC excluded", which(!d$mhc))
  ldk <- fread(file.path(DIR_INSTRUMENTS, "ld_pruned", paste0(s, ".tsv")))$rsid
  if (!all(d$rsid %in% ldk)) addv("residual-LD pruned", which(d$rsid %in% ldk))
  if (!is.na(D$n_mvp_r2_lt_0.8) && D$n_mvp_r2_lt_0.8 > 0) addv("MVP imputation r2 >= 0.8", which(d$r2 >= 0.8))
  if (n_rev[["primary"]] > 0) addv("Steiger-filtered", which(!((!sp$correct) & sp$p_snp < 0.05)))
  if (o == "pad_vanzuydam") addv("van Zuydam conservative SE", seq_len(L), sy_alt = pmax(d$se_y, d$se_y_sens, na.rm = TRUE))
  list(E = E, D = D, V = rbindlist(V), loo = loo, ss = ss, st = st_snp)
}

grid <- CJ(s = names(sets), o = outcomes, sorted = FALSE); grid[, seed := PROJECT_SEED + .I]
stopifnot(nrow(grid) == 96)
cell_file <- function(i, dir = DIR_CELLS) file.path(dir, paste0(grid$s[i], "__", grid$o[i], ".rds"))
do_cell <- function(i, dir = DIR_CELLS) {
  f <- cell_file(i, dir)
  if (file.exists(f)) return(invisible("skipped"))
  r <- run_cell(grid$s[i], grid$o[i], grid$seed[i])
  r$meta <- list(cell = i, set = grid$s[i], outcome = grid$o[i], seed = grid$seed[i], R = R.version.string,
                 pkgs = sapply(c("MendelianRandomization", "MRPRESSO", "mr.raps", "TwoSampleMR", "psych", "data.table"),
                               function(p) as.character(packageVersion(p))), finished = format(Sys.time()))
  fs::dir_create(dir); tmp <- paste0(f, ".tmp"); saveRDS(r, tmp); file.rename(tmp, f)
  invisible("done")
}
merge_cells <- function(dir = DIR_CELLS, out = DIR_TABLES) {
  fl <- vapply(seq_len(nrow(grid)), cell_file, "", dir = dir)
  miss <- !file.exists(fl)
  if (any(miss)) stop(sum(miss), " of 96 cell files missing, e.g. ", basename(fl[which(miss)[1]]))
  res <- lapply(fl, readRDS)
  pick <- function(k) rbindlist(lapply(res, `[[`, k), fill = TRUE)
  E <- pick("E"); E[, `:=`(or = exp(b), or_lo = exp(lo), or_hi = exp(hi))]
  D <- pick("D")
  # Reporting fix (2026-10-08, after the UGER run): cells computed before the run_cell fix counted MRPRESSO's string
  # "No significant outliers" as 1 outlier. That string comes with Distortion p = NA and no outlier-corrected estimate.
  nso <- D$presso_n_outliers == 1 & is.na(D$presso_distortion_p) & is.na(D$presso_error)
  D[nso, presso_n_outliers := 0L]
  if (any(nso)) logl("merge: presso_n_outliers corrected 1 -> 0 ('No significant outliers') in {sum(nso)} cells")
  V <- pick("V"); V[, `:=`(or = exp(b), or_lo = exp(lo), or_hi = exp(hi))]
  fwrite(E, file.path(out, "mr_univariable_estimates.tsv"), sep = "\t")
  fwrite(D, file.path(out, "mr_univariable_diagnostics.tsv"), sep = "\t")
  fwrite(V, file.path(out, "mr_univariable_variants.tsv"), sep = "\t")
  fwrite(pick("loo"), file.path(out, "mr_univariable_loo.tsv"), sep = "\t")
  fwrite(pick("ss"), file.path(out, "mr_univariable_singlesnp.tsv"), sep = "\t")
  fwrite(pick("st"), file.path(out, "mr_univariable_steiger_snp.tsv"), sep = "\t")
  fwrite(rbindlist(lapply(res, function(x) as.data.table(c(x$meta[c("cell", "set", "outcome", "seed", "R", "finished")], as.list(x$meta$pkgs))))),
         file.path(out, "mr_univariable_cell_provenance.tsv"), sep = "\t")
  P <- E[role == "primary" & method == "IVW (MRE)"]
  P[, `:=`(alpha = 0.05 / 16, bonferroni_sig = p < 0.05 / 16,
           or_lo_bonf = exp(b - qnorm(1 - 0.05 / 32) * se), or_hi_bonf = exp(b + qnorm(1 - 0.05 / 32) * se))]
  stopifnot(nrow(P) == 16)
  fwrite(P, file.path(out, "mr_primary_grid.tsv"), sep = "\t")
  logl("merge: {nrow(E)} estimate rows, {nrow(D)} cells, {nrow(V)} variant rows; primary grid 16 cells")
}
# reference cells for --check: one small and one PRESSO-outlier-test cell, chosen for runtime, not for results
REF_CELLS <- which(grid$s == "TDI" & grid$o == "hf_mvp_phe428_eur" | grid$s == "Kweon" & grid$o == "pad_mvp_eur")
compare_cells <- function(new, ref) {
  num_equal <- function(a, b) {
    for (k in names(b)) if (is.numeric(b[[k]])) { d <- abs(a[[k]] - b[[k]]); d[is.na(a[[k]]) & is.na(b[[k]])] <- 0
      if (any(is.na(d)) || any(d > 1e-8 * pmax(1, abs(b[[k]]), na.rm = TRUE))) return(paste("column", k)) }
      else if (!identical(as.character(a[[k]]), as.character(b[[k]]))) return(paste("column", k))
    NA_character_ }
  for (k in c("E", "D", "V", "loo", "ss", "st")) {
    a <- new[[k]]; b <- ref[[k]]
    if (!identical(dim(a), dim(b))) return(paste(k, "dimensions differ"))
    w <- num_equal(a, b); if (!is.na(w)) return(paste(k, w))
  }
  NA_character_
}

if (MODE == "prep") { prep_n()
} else {
  if (MODE == "local" && !file.exists(file.path(DIR_CIN, "outcome_n_totals.tsv"))) prep_n()
  NT <- load_n()
  if (MODE == "cell") {
    stopifnot(!is.na(CELL_K), CELL_K >= 1, CELL_K <= 96)
    logl("cell {CELL_K}: {grid$s[CELL_K]} x {grid$o[CELL_K]}; PRESSO {N_PRESSO}; bootstrap {N_BOOT}")
    st <- do_cell(CELL_K); logl("cell {CELL_K}: {st}")
  } else if (MODE == "merge") { merge_cells()
  } else if (MODE == "check") {
    tmpd <- tempfile("check"); bad <- character()
    for (i in REF_CELLS) {
      ref <- readRDS(cell_file(i, DIR_REF)); do_cell(i, tmpd); new <- readRDS(cell_file(i, tmpd))
      w <- compare_cells(new, ref); logl("check cell {i} ({grid$s[i]} x {grid$o[i]}): {ifelse(is.na(w), 'identical (tolerance 1e-8)', paste('DIFFERS:', w))}")
      if (!is.na(w)) bad <- c(bad, w)
    }
    if (length(bad)) stop("CHECK FAILED") else logl("CHECK PASSED")
  } else {
    todo <- which(!file.exists(vapply(seq_len(nrow(grid)), cell_file, "")))
    logl("local: {length(todo)} of 96 cells to run; cores {N_CORES}; PRESSO {N_PRESSO}; bootstrap {N_BOOT}")
    res <- parallel::mclapply(todo, function(i) tryCatch(do_cell(i), error = function(e) structure(list(msg = conditionMessage(e)), class = "cell_error")),
                              mc.cores = N_CORES, mc.preschedule = FALSE)
    bad <- which(sapply(res, inherits, "cell_error"))
    if (length(bad)) stop("cells failed: ", paste(grid$s[todo[bad]], grid$o[todo[bad]], sapply(res[bad], `[[`, "msg"), sep = " / ", collapse = "; "))
    merge_cells()
  }
}
close_run_log(log_con)
