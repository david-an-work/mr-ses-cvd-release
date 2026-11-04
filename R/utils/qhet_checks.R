#' qhet_checks.R — numerical checks around MVMR::qhet_mvmr (0.4.8)
#'
#' The package function (i) profiles a heterogeneity term tau2 with stats::optimize on (-10, 10), where the inner
#' minimisation is Nelder-Mead from 0 and the outer target is (Q - (L - 2))^2 whatever the number of exposures K;
#' (ii) re-minimises Q at that tau2 by Nelder-Mead from 0; (iii) reports no convergence diagnostic; (iv) bootstraps by
#' resampling SNPs (boot::boot) with BCa intervals rounded to 3 decimals and no handling of failed resamples.
#' Status (2026-10-06): under the usability rule below, no package fit was usable in simulated tests (R/02h), so
#' qhet_mvmr is deferred from analysis plan v0.5. Kept as documentation of that decision.
#' These helpers re-implement the objective to check the package's solution, and wrap the bootstrap so that failures
#' are counted. They do not change the estimator.
#'
#' qhet_reimpl(r_input, pcor)        -> list(b, tau2, Q, w, conv, n_inner, n_inner_fail): faithful re-implementation of the
#'                                       package point estimate, plus the optimizers' own convergence codes
#' qhet_valid(re)                     -> character reasons why a re-implemented fit is NOT usable (empty if usable)
#' qhet_Q(b, tau2, r_input, pcor)     -> Q at (b, tau2)
#' qhet_conv_check(r_input, pcor, b_candidate = NULL, n_starts = 10, seed = 1, rel_tol = 1e-6, tau2_override = NULL)
#'                                    -> one-row data.table with status "usable" / "not usable" / "not estimable";
#'                                       `reasons` lists every failed condition (tau2_override: TEST ONLY)
#'   Also reports tau2 re-solved at tol 1e-10 and the change in b at that tau2 (diagnostic; see plan §10.3).
#' qhet_boot_checked(r_input, pcor, R = 1000, ncores = 1, seed = 1, max_fail = 0.05, fail_hook = NULL)
#'                                    -> list(summary = data.table, boot = boot object)
#'   Each resample is fitted with qhet_reimpl and counted as FAILED unless qhet_valid() is empty (same rule as the point
#'   estimate). CI status: no failures -> "ok" (boot::boot.ci BCa, as the package); 0 < failures <= max_fail ->
#'   "failure-affected" (BCa on successful resamples, jackknife acceleration; exploratory only, since failed resamples may
#'   differ systematically); failures > max_fail -> "withheld" (no interval).
#'   fail_hook: TEST ONLY. A function(indices) returning TRUE to force a failure in that resample.

.qhet_parts <- function(r_input, pcor) {
  K <- (ncol(r_input) - 3) / 2
  pihat <- as.matrix(r_input[, 4:(3 + K)])
  se <- as.matrix(r_input[, (K + 4):ncol(r_input)])
  covlist <- lapply(seq_len(nrow(r_input)), function(l) pcor * outer(se[l, ], se[l, ]))
  list(K = K, L = nrow(r_input), pihat = pihat, gam = r_input$betaYG, segam = r_input$sebetaYG, covlist = covlist)
}

.qhet_w <- function(b, tau2, P) P$segam^2 + vapply(P$covlist, function(cv) drop(t(b) %*% cv %*% b), 0) + tau2

qhet_Q <- function(b, tau2, r_input, pcor, P = .qhet_parts(r_input, pcor)) {
  w <- .qhet_w(b, tau2, P); sum((1 / w) * (P$gam - P$pihat %*% b)^2)
}

qhet_reimpl <- function(r_input, pcor, tau2_override = NULL) {
  P <- .qhet_parts(r_input, pcor)
  inner_codes <- integer(0)
  PL <- function(tau2) {
    f <- stats::optim(rep(0, P$K), function(b) qhet_Q(b, tau2, P = P))
    inner_codes <<- c(inner_codes, f$convergence)
    (qhet_Q(f$par, tau2, P = P) - (P$L - 2))^2
  }
  tau2 <- if (is.null(tau2_override)) stats::optimize(PL, interval = c(-10, 10))$minimum else tau2_override
  fit <- stats::optim(rep(0, P$K), function(b) qhet_Q(b, tau2, P = P))
  list(b = fit$par, tau2 = tau2, Q = fit$value, w = .qhet_w(fit$par, tau2, P), conv = fit$convergence,
       n_inner = length(inner_codes), n_inner_fail = sum(inner_codes != 0), P = P)
}

# Conditions for a usable fit (applied to the point estimate AND to every bootstrap resample)
qhet_valid <- function(re) {
  r <- character(0)
  if (any(!is.finite(re$b)) || !is.finite(re$Q)) r <- c(r, "non-finite estimate")
  if (re$conv != 0) r <- c(r, paste0("final optim code ", re$conv))
  if (re$n_inner_fail > 0) r <- c(r, paste0(re$n_inner_fail, "/", re$n_inner, " profile optim codes != 0"))
  if (any(!is.finite(re$w)) || any(re$w <= 0)) r <- c(r, "non-positive weight denominator")
  if (abs(abs(re$tau2) - 10) < 1e-3) r <- c(r, "tau2 at search bound")
  r
}

qhet_conv_check <- function(r_input, pcor, b_candidate = NULL, n_starts = 10, seed = 1, rel_tol = 1e-6,
                            tau2_override = NULL) {
  pkg <- tryCatch(suppressWarnings(MVMR::qhet_mvmr(r_input, pcor = pcor, CI = FALSE, iterations = 0))[, 1],
                  error = function(e) structure(NA_real_, err = conditionMessage(e)))
  if (anyNA(pkg) || any(!is.finite(pkg)))
    return(data.table::data.table(status = "not estimable", reasons = attr(pkg, "err") %||% "non-finite estimate"))
  re <- qhet_reimpl(r_input, pcor, tau2_override = tau2_override); P <- re$P
  reimpl_maxdiff <- if (is.null(tau2_override)) max(abs(re$b - pkg)) else NA_real_
  b0 <- if (is.null(b_candidate)) re$b else b_candidate
  if (!is.null(b_candidate)) re$w <- .qhet_w(b0, re$tau2, P)
  Q0 <- qhet_Q(b0, re$tau2, P = P)
  ivw <- unname(stats::lm(P$gam ~ 0 + P$pihat, weights = 1 / P$segam^2)$coefficients)
  set.seed(seed)
  scale <- pmax(abs(ivw), 0.5)
  starts <- c(list(ivw, b0), lapply(seq_len(n_starts), function(i) ivw + stats::runif(P$K, -2, 2) * scale))
  fits <- lapply(starts, function(s) {
    f1 <- stats::optim(s, function(b) qhet_Q(b, re$tau2, P = P), control = list(maxit = 5000, reltol = 1e-12))
    f2 <- tryCatch(stats::optim(f1$par, function(b) qhet_Q(b, re$tau2, P = P), method = "BFGS",
                                control = list(maxit = 1000, reltol = 1e-14)), error = function(e) f1)
    if (f2$value <= f1$value) f2 else f1
  })
  Qs <- vapply(fits, `[[`, 0, "value"); best <- fits[[which.min(Qs)]]
  rel_gap <- (Q0 - min(Qs)) / max(abs(min(Qs)), 1e-12)
  # tau2 precision: the package uses optimize() with its default tol (~1.2e-4 on tau2), which can be as large as the
  # SNP-outcome variances. Re-solve tau2 with tol = 1e-10 (same inner routine) and re-minimise b at the refined tau2.
  PLr <- function(tau2) {
    bc <- stats::optim(rep(0, P$K), function(b) qhet_Q(b, tau2, P = P))$par
    (qhet_Q(bc, tau2, P = P) - (P$L - 2))^2
  }
  tau2_ref <- stats::optimize(PLr, interval = c(-10, 10), tol = 1e-10)$minimum
  b_ref <- stats::optim(b0, function(b) qhet_Q(b, tau2_ref, P = P), control = list(maxit = 5000, reltol = 1e-12))$par
  reasons <- qhet_valid(re)
  if (rel_gap > rel_tol) reasons <- c(reasons, sprintf("objective %.2g above best of %d starts", rel_gap, length(starts)))
  data.table::data.table(
    status = if (length(reasons)) "not usable" else "usable", reasons = paste(reasons, collapse = "; "),
    rel_Q_gap = rel_gap, max_abs_b_diff = max(abs(best$par - b0)), final_optim_code = re$conv,
    profile_optim_fail = re$n_inner_fail, profile_optim_n = re$n_inner,
    reimpl_maxdiff = reimpl_maxdiff, tau2 = re$tau2, tau2_refined = tau2_ref,
    max_abs_b_diff_refined_tau2 = max(abs(b_ref - b0)), tau2_negative = re$tau2 < 0,
    min_weight_denominator = min(re$w), df_used_by_pkg = P$L - 2, df_correct = P$L - P$K,
    detail = paste0("b_pkg=", paste(signif(pkg, 6), collapse = ","), "; b_best=", paste(signif(best$par, 6), collapse = ",")))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

qhet_boot_checked <- function(r_input, pcor, R = 1000, ncores = 1, seed = 1, max_fail = 0.05, fail_hook = NULL) {
  K <- (ncol(r_input) - 3) / 2
  n <- nrow(r_input)
  stat <- function(d, i, hook = TRUE) {
    if (hook && !is.null(fail_hook) && !identical(as.integer(i), seq_len(n)) && isTRUE(fail_hook(i))) return(rep(NA_real_, K))
    re <- tryCatch(qhet_reimpl(d[i, ], pcor), error = function(e) NULL)
    if (is.null(re) || length(qhet_valid(re))) rep(NA_real_, K) else re$b
  }
  set.seed(seed)
  b <- if (ncores > 1) boot::boot(r_input, stat, R = R, parallel = "multicore", ncpus = ncores) else boot::boot(r_input, stat, R = R)
  t0_ok <- all(is.finite(b$t0))
  ok <- stats::complete.cases(b$t); n_fail <- sum(!ok); frac <- n_fail / R
  status <- if (!t0_ok) "withheld (point estimate not usable)" else if (frac > max_fail) "withheld" else if (n_fail > 0) "failure-affected" else "ok"
  has_ci <- status %in% c("ok", "failure-affected") && sum(ok) >= 50
  jackm <- if (has_ci && n_fail > 0) t(vapply(seq_len(n), function(j) stat(r_input, seq_len(n)[-j], hook = FALSE), numeric(K))) else NULL
  if (!is.null(jackm) && K == 1) jackm <- t(jackm)
  ci <- lapply(seq_len(K), function(k) {
    if (!has_ci) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
    if (n_fail == 0) {  # identical to the package's own interval (before its 3-decimal rounding)
      bc <- boot::boot.ci(b, type = c("bca", "perc"), index = k)
      return(c(bc$bca[4:5], bc$percent[4:5]))
    }
    tk <- b$t[ok, k]
    pc <- stats::quantile(tk, c(0.025, 0.975), names = FALSE)
    # BCa on the successful resamples (exploratory): z0 from the bootstrap distribution, acceleration from the SNP jackknife
    jack <- jackm[, k]; jack <- jack[is.finite(jack)]
    z0 <- stats::qnorm(mean(tk < b$t0[k]))
    Lj <- mean(jack) - jack; a <- sum(Lj^3) / (6 * sum(Lj^2)^1.5)
    zq <- stats::qnorm(c(0.025, 0.975)); adj <- stats::pnorm(z0 + (z0 + zq) / (1 - a * (z0 + zq)))
    c(stats::quantile(tk, adj, names = FALSE, type = 6), pc)
  })
  ci <- do.call(rbind, ci)
  list(summary = data.table::data.table(exposure = seq_len(K), estimate = b$t0, bca_lo = ci[, 1], bca_hi = ci[, 2],
                                         pct_lo = ci[, 3], pct_hi = ci[, 4], R = R, n_fail = n_fail,
                                         fail_frac = frac, ci_status = status,
                                         ci_method = if (!has_ci) "none" else if (n_fail == 0) "boot.ci_bca" else "bca_subset_exploratory",
                                         n_jack_fail = if (is.null(jackm)) NA_integer_ else sum(!stats::complete.cases(jackm))),
       boot = b)
}
