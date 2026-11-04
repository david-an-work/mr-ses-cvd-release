#' mr_methods.R
#'
#' Univariable MR estimators and diagnostics for plan v1.0 §7–§8 (Phase 4). Inputs are vectors on the harmonized scale:
#' bx, sx (exposure, per SD or per released unit), by, sy (outcome log-OR, aligned to the exposure-increasing allele).
#' All CIs are normal-approximation 95% intervals unless the method defines its own.

Z975 <- qnorm(0.975)

#' IVW, multiplicative random effects (first-order weights 1/sy^2; residual SE floored at 1).
ivw_mre <- function(bx, by, sy) {
  w <- 1 / sy^2; L <- length(bx)
  b <- sum(w * bx * by) / sum(w * bx^2)
  se_fe <- 1 / sqrt(sum(w * bx^2))
  Q <- sum(w * (by - b * bx)^2); phi <- if (L > 1) Q / (L - 1) else NA_real_
  se <- se_fe * max(1, sqrt(ifelse(is.na(phi), 1, phi)))
  list(b = b, se = se, se_fe = se_fe, Q = Q, Q_df = L - 1, Q_p = if (L > 1) pchisq(Q, L - 1, lower.tail = FALSE) else NA_real_,
       I2 = if (L > 1) max(0, (Q - (L - 1)) / Q) else NA_real_, nsnp = L)
}

est_row <- function(method, b, se, nsnp, lo = b - Z975 * se, hi = b + Z975 * se, p = 2 * pnorm(-abs(b / se))) {
  data.table::data.table(method = method, nsnp = nsnp, b = b, se = se, lo = lo, hi = hi, p = p)
}

#' Leave-one-out IVW (MRE): returns per-SNP estimates
loo_ivw <- function(rsid, bx, by, sy) {
  data.table::rbindlist(lapply(seq_along(rsid), function(i) {
    r <- ivw_mre(bx[-i], by[-i], sy[-i]); data.table::data.table(left_out = rsid[i], b = r$b, se = r$se)
  }))
}

#' Steiger (Hemani 2017 as implemented in TwoSampleMR): per-SNP and overall, at one prevalence.
#' r2 exposure from beta, se, n (continuous); r2 outcome from log-OR on the liability scale (logit model).
steiger_one <- function(bx, sx, nx, by, eaf, ncase, ncontrol, prevalence) {
  r_exp <- TwoSampleMR::get_r_from_bsen(bx, sx, nx)
  r_out <- TwoSampleMR::get_r_from_lor(by, eaf, ncase, ncontrol, prevalence, model = "logit")
  n_out <- ncase + ncontrol
  per <- psych::r.test(n = nx, n2 = n_out, r12 = abs(r_exp), r34 = abs(r_out))
  stopifnot(!anyNA(r_exp), !anyNA(r_out))  # r supplied for every SNP, so mr_steiger never falls back to p-values
  overall <- TwoSampleMR::mr_steiger(p_exp = NA, p_out = NA, n_exp = nx, n_out = n_out, r_exp = r_exp, r_out = r_out)
  list(r2_exp = r_exp^2, r2_out = r_out^2, p_snp = per$p, correct = r_exp^2 > r_out^2,
       overall_correct = overall$correct_causal_direction, overall_p = overall$steiger_test,
       r2_exp_sum = overall$r2_exp, r2_out_sum = overall$r2_out)
}
