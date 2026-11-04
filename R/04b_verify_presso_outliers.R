#' 04b_verify_presso_outliers.R
#'
#' Purpose: verify the 12 MR-PRESSO outlier counts that `R/04 merge` corrected from 1 to 0 ("No significant outliers"),
#'   by re-running MR-PRESSO on exactly the same inputs and seed as R/04 (seed = PROJECT_SEED + grid index; MRPRESSO
#'   sets this seed itself) and recording its raw outlier-test output. A cell is confirmed when the re-run global-test
#'   p-value equals the stored one and `Outliers Indices` is the string "No significant outliers" (or the outlier test
#'   finds no SNP with Bonferroni p <= 0.05). No estimate is changed or added.
#' Usage: Rscript R/04b_verify_presso_outliers.R K   (K = 1..12; UGER array task) | merge
#' Outputs: results/cells/presso_verify/{set}__{outcome}.rds ; results/tables/mr_presso_outlier_verification.tsv
#' Depends on: data.table, MRPRESSO

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
args <- commandArgs(trailingOnly = TRUE); MODE <- if (length(args)) args[1] else stop("give K or merge")
log_con <- log_run(paste0("04b_verify_presso_", MODE))
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
H <- file.path(DIR_DATA_HARM, "large"); OUT <- file.path(DIR_RESULTS, "cells", "presso_verify"); fs::dir_create(OUT)
sets <- c("EA4", "Kweon", "TDI", "CAMSIS", "EA3", "Hill", "TDI_EUR", "gSES")
outcomes <- c("cad_aragam_eur", "cad_aragam_multi", "cad_mvp_phe411_eur", "pad_mvp_eur", "pad_finngen_athscle",
              "pad_vanzuydam", "pad_mvp_afr", "hf_hermes2", "hf_mvp_phe428_eur", "stroke_gigastroke_multi",
              "stroke_gigastroke_eur", "stroke_mvp_phe433_21_eur")
grid <- CJ(s = sets, o = outcomes, sorted = FALSE); grid[, seed := PROJECT_SEED + .I]   # identical to R/04
# the 12 cells corrected by R/04 merge (presso_n_outliers 1 -> 0)
targets <- data.table(s = c("EA4", "EA4", "EA4", "Kweon", "Kweon", "EA3", "EA3", "TDI_EUR", "TDI_EUR", "gSES", "gSES", "gSES"),
                      o = c("pad_mvp_eur", "stroke_gigastroke_eur", "stroke_mvp_phe433_21_eur", "pad_mvp_afr", "stroke_gigastroke_multi",
                            "stroke_gigastroke_eur", "stroke_mvp_phe433_21_eur", "stroke_gigastroke_multi", "stroke_mvp_phe433_21_eur",
                            "stroke_gigastroke_multi", "stroke_gigastroke_eur", "stroke_mvp_phe433_21_eur"))
targets <- merge(targets, grid[, .(s, o, seed, cell = .I)], by = c("s", "o"), sort = FALSE)
stopifnot(nrow(targets) == 12)

if (MODE == "merge") {
  D <- fread(file.path(DIR_TABLES, "mr_univariable_diagnostics.tsv"))
  res <- rbindlist(lapply(seq_len(nrow(targets)), function(k) {
    f <- file.path(OUT, paste0(targets$s[k], "__", targets$o[k], ".rds")); stopifnot(file.exists(f)); readRDS(f)
  }))
  res <- merge(res, D[, .(set, outcome, stored_global_p = as.character(presso_global_p), stored_n_outliers = presso_n_outliers)],
               by = c("set", "outcome"))
  res[, global_p_matches := rerun_global_p == stored_global_p]
  res[, confirmed := global_p_matches & n_numeric_outliers == 0 & stored_n_outliers == 0]
  fwrite(res, file.path(DIR_TABLES, "mr_presso_outlier_verification.tsv"), sep = "\t")
  logl("confirmed {sum(res$confirmed)} of {nrow(res)}")
  if (!all(res$confirmed)) stop("not all corrected PRESSO counts confirmed")
} else {
  k <- as.integer(MODE); t <- targets[k]
  h <- fread(file.path(H, "univariable", paste0(t$s, "__", t$o, ".tsv.gz")), colClasses = list(character = c("ea_x", "oa_x", "ea_y", "oa_y", "chr")))
  d <- h[mr_keep == TRUE]
  set.seed(t$seed)
  pr <- suppressWarnings(MRPRESSO::mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sy", SdExposure = "sx",
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = data.frame(bx = d$beta_x, by = d$beta_y_aligned, sx = d$se_x, sy = d$se_y),
          NbDistribution = 10000L, SignifThreshold = 0.05, seed = t$seed))
  ot <- pr$`MR-PRESSO results`$`Outlier Test`; dt <- pr$`MR-PRESSO results`$`Distortion Test`
  oi <- if (is.null(dt)) NA_character_ else paste(dt$`Outliers Indices`, collapse = ",")
  r <- data.table(set = t$s, outcome = t$o, cell = t$cell, seed = t$seed, nsnp = nrow(d),
                  rerun_global_p = as.character(pr$`MR-PRESSO results`$`Global Test`$Pvalue),
                  outliers_indices_raw = oi,
                  n_numeric_outliers = if (is.null(dt) || !is.numeric(dt$`Outliers Indices`)) 0L else length(dt$`Outliers Indices`),
                  smallest_outlier_test_p = if (is.null(ot)) NA_character_ else paste(head(ot$Pvalue[order(suppressWarnings(as.numeric(sub("<", "", ot$Pvalue))))], 3), collapse = ";"))
  saveRDS(r, file.path(OUT, paste0(t$s, "__", t$o, ".rds")))
  logl("{t$s} x {t$o}: global p {r$rerun_global_p}; outliers '{oi}'; smallest outlier-test p {r$smallest_outlier_test_p}")
}
close_run_log(log_con)
