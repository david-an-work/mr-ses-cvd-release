#' 02g_vanzuydam_se_scale.R
#'
#' Purpose: decide, with reviewable variant-level evidence, whether the van Zuydam 2021 PAD column `or_se` is the
#'   standard error of the OR (delta method: SE(OR) = OR x SE(log OR)) or of the log-OR. No MR estimate is computed;
#'   only the outcome file's own columns are used.
#' Rounding: `or` and `or_se` are deposited at 2 decimals; `p` at 2 significant figures (e.g. 2.4e-05, 3e-04, 0.54).
#'   For each variant, true OR lies in [or - 0.005, or + 0.005] and true p in [p - h, p + h], h = 0.5 x 10^(floor(log10 p) - 1).
#'   Implied SE(log OR) = |ln OR| / z(p), z(p) = qnorm(p/2, lower = FALSE). Each hypothesis predicts a range for the
#'   unrounded `or_se`: H_log: |ln OR|/z; H_or: OR x |ln OR|/z, evaluated over the OR x p rounding box (41 x 41 grid,
#'   corners included). The deposited `or_se` is compatible with a hypothesis if its rounding interval
#'   [or_se - 0.005, or_se + 0.005] intersects the predicted range.
#' Variant set (fixed, as in the first ad hoc check of 2026-10-06): p < 1e-6 and (or >= 1.20 or or <= 0.83), i.e. large
#'   effects where the two hypotheses predict different values. A second, pre-stated wider set (p < 1e-4 and
#'   |ln or| >= ln 1.10) is reported separately as a robustness check.
#' Interpretation (2026-10-06 run): or_se fits the OR scale far better than the log-OR scale, but p-values are
#'   systematically less extreme than or/or_se imply (median or_se / implied SE(OR) < 1), so p and or_se are not
#'   mutually consistent for a minority of variants under either hypothesis. Conversion SE(logOR) = or_se/OR is
#'   therefore "supported but provisional"; a p-derived SE is used as a sensitivity analysis (plan v0.4 §4.2).
#' Supersedes the ad hoc Python check of 2026-10-06 (14 informative, 9 SE(OR) matches, 0 SE(logOR) matches), whose
#'   point-rounding comparison differs slightly from this script's because of rounding implementation.
#' Outputs: results/tables/vanzuydam_se_scale_check.{tsv,xlsx} (variant level), vanzuydam_se_scale_summary.{tsv,xlsx};
#'   results/logs/{date}_02g_vanzuydam_se_scale.log
#' Depends on: data.table

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })

log_con <- log_run("02g_vanzuydam_se_scale")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)

f <- here::here("data", "raw", "pad_vanzuydam2021", "final_meta", "pad_primary.txt.gz")
logl("input: {f}; md5 {unname(tools::md5sum(f))}")
# the header has 11 names but rows have 12 fields (leading row index); fread names the extra first column V1
d <- fread(f, select = c("chr", "bp", "snp", "or", "or_se", "p"), colClasses = list(character = "p"))
d[, p_num := as.numeric(p)]

sets <- list(
  primary = d[p_num < 1e-6 & (or >= 1.20 | or <= 0.83)],
  wider   = d[p_num < 1e-4 & abs(log(or)) >= log(1.10)]
)

check_set <- function(v, set_name) {
  v <- copy(v)
  v[, p_half := 0.5 * 10^(floor(log10(p_num)) - 1)]
  g <- seq(0, 1, length.out = 41)
  res <- v[, {
    ors <- or + (g - 0.5) * 0.01
    ps  <- pmax(p_num + (g - 0.5) * 2 * p_half, .Machine$double.xmin)
    grid <- CJ(o = ors, pp = ps)
    z <- qnorm(grid$pp / 2, lower.tail = FALSE)
    se_log <- abs(log(grid$o)) / z
    se_or  <- grid$o * se_log
    se_log_point <- abs(log(or)) / qnorm(p_num / 2, lower.tail = FALSE)
    .(se_log_point = se_log_point,
      pred_log_round = round(se_log_point, 2), pred_or_round = round(or * se_log_point, 2),
      log_lo = min(se_log), log_hi = max(se_log), or_lo = min(se_or), or_hi = max(se_or))
  }, by = .(chr, bp, snp, or, or_se, p, p_num)]
  res[, `:=`(
    point_informative = pred_log_round != pred_or_round,
    point_match_log = pred_log_round == or_se,
    point_match_or  = pred_or_round == or_se,
    compat_log = (or_se + 0.005) >= log_lo & (or_se - 0.005) <= log_hi,
    compat_or  = (or_se + 0.005) >= or_lo  & (or_se - 0.005) <= or_hi
  )]
  res[, interval_class := fifelse(compat_or & !compat_log, "SE(OR) only",
                          fifelse(compat_log & !compat_or, "SE(logOR) only",
                          fifelse(compat_or & compat_log, "both", "neither")))]
  res[, set := set_name]
  res[]
}

out <- rbindlist(lapply(names(sets), function(s) check_set(sets[[s]], s)))
summ <- out[, .(n_variants = .N,
                n_point_informative = sum(point_informative),
                point_match_or = sum(point_informative & point_match_or),
                point_match_log = sum(point_informative & point_match_log),
                interval_or_only = sum(interval_class == "SE(OR) only"),
                interval_log_only = sum(interval_class == "SE(logOR) only"),
                interval_both = sum(interval_class == "both"),
                interval_neither = sum(interval_class == "neither"),
                median_ratio_to_SE_OR = median(or_se / (or * se_log_point)),
                median_ratio_to_SE_logOR = median(or_se / se_log_point),
                q25_ratio_to_SE_OR = quantile(or_se / (or * se_log_point), 0.25),
                q75_ratio_to_SE_OR = quantile(or_se / (or * se_log_point), 0.75)), by = set]
write_result_table(out, "vanzuydam_se_scale_check")
write_result_table(summ, "vanzuydam_se_scale_summary")
logl("summary:\n{paste(capture.output(print(summ)), collapse = '\n')}")
logl("Done.")
close_run_log(log_con)
