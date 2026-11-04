#' 02h_qhet_checks_test.R
#'
#' Purpose: test the numerical checks in R/utils/qhet_checks.R on SIMULATED data before plan v0.4 describes them as
#'   operational. No project data are used, and no project estimate is computed. Coverage of the bootstrap CIs is NOT
#'   assessed here (that would need many CI runs per regime); this script tests that the checks run and flag correctly.
#' Tests (fixed before running):
#'   T1  Re-implementation and convergence: 10 datasets per regime (2 exposures, L = 400; weak tau 0.0068, strong
#'       0.026, the 02f design) and 10 datasets with 3 exposures (L = 400). Record status, relative Q gap, max |b diff|,
#'       re-implementation difference from the package, tau2 sign/bounds, tau2 re-solved at tol 1e-10 and the
#'       resulting change in b, weight positivity, df used vs correct.
#'   T2  Flag test: on 3 datasets, a deliberately perturbed candidate (package solution + 0.05 per exposure) must be
#'       "not usable"; the package solution itself must be "usable" (unless an optimizer code fails, which is recorded).
#'   T4  (added 2026-10-06, status-rule revision) Forced tau2 = -2e-4 (below -min SNP-outcome variance, 1e-4) on 3
#'       datasets: weight denominators <= 0, so the fit must be "not usable".
#'   T3  Bootstrap (L = 100 so that runs are short): (a) no failures, ncores = 1, R = 100: our CI must equal the
#'       package's own qhet_mvmr(CI = TRUE) interval to its 3-decimal rounding, same seed; (b) forced failure rate 3%
#'       (test hook), R = 200: status "failure-affected" (exploratory interval); (c) forced 10%: status "withheld".
#'       Forced failures test software handling only, not interval validity.
#' Outcome (2026-10-06, strict status rules): T1 0/30 usable. The only failing condition in all 30 was the inner
#'   profile fits (Nelder-Mead code 10, degenerate simplex) at 4-7 of 28 trial tau2 values. A diagnostic run shows these
#'   occur at NEGATIVE trial tau2 where weight denominators collapse to ~0 (objective ~1e38); the final fit converged
#'   (code 0) in 30/30, re-implementation difference 0, max relative Q gap 1.2e-7. Package tau2 = +/-2.28e-5 in 30/30;
#'   refined-tau2 change in b median 0.023 / 0.005 / 0.027 (weak / strong / three-exposure), max 0.075.
#'   T2: perturbed candidates "not usable" 3/3; package solutions also "not usable" (same profile reason), so T2 flags
#'   FALSE against its pre-run expectation. T4: forced tau2 = -2e-4 "not usable" 3/3 (non-positive weights).
#'   T3: every bootstrap resample failed the same rule (100/100, 200/200), so no interval was produced.
#'   Conclusion: the rule was not relaxed after seeing these results; qhet_mvmr is deferred from plan v0.5 (§10.3 B).
#' Outputs: results/tables/qhet_checks_test_{convergence,flag,nonpositive_w,bootstrap}.{tsv,xlsx}; results/logs/{date}_02h_qhet_checks_test.log
#' Depends on: MVMR, boot, data.table, MASS

source(here::here("R", "00_setup.R"))
source(file.path(DIR_R_UTILS, "qhet_checks.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
ensure_pkgs_loaded(c("MVMR", "boot", "data.table", "MASS"))
log_con <- log_run("02h_qhet_checks_test")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
set.seed(20261006)
NC <- max(1L, parallelly::availableCores(omit = 1) - 1L)

RHO <- 0.23; SEX <- 0.005; SEY <- 0.01
sim <- function(L, tau, K = 2) {
  g1 <- rnorm(L, 0, 0.02)
  G <- cbind(g1, 0.9 * g1 + rnorm(L, 0, tau))
  theta <- c(0.3, -0.1)
  if (K == 3) { G <- cbind(G, 0.8 * g1 + rnorm(L, 0, 0.015)); theta <- c(0.3, -0.1, 0.05) }
  PC <- matrix(RHO, K, K); diag(PC) <- 1
  bx <- G + MASS::mvrnorm(L, rep(0, K), PC * SEX^2)
  by <- as.vector(G %*% theta) + rnorm(L, 0, SEY)
  list(r = MVMR::format_mvmr(BXGs = bx, BYG = by, seBXGs = matrix(SEX, L, K), seBYG = rep(SEY, L),
                             RSID = paste0("s", seq_len(L))), pc = PC)
}

# T1
t1 <- rbindlist(lapply(list(list("weak", 0.0068, 2), list("strong", 0.026, 2), list("three_exposure", 0.0068, 3)),
  function(rg) rbindlist(lapply(1:10, function(i) {
    s <- sim(400, rg[[2]], rg[[3]]); t0 <- Sys.time()
    cbind(regime = rg[[1]], rep = i, qhet_conv_check(s$r, s$pc, seed = i),
          seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  }))))
logl("T1 status counts:\n{paste(capture.output(print(t1[, .N, by = .(regime, status)])), collapse = '\n')}")

# T2
t2 <- rbindlist(lapply(1:3, function(i) {
  s <- sim(400, 0.0068)
  pkg <- suppressWarnings(MVMR::qhet_mvmr(s$r, pcor = s$pc, CI = FALSE, iterations = 0))[, 1]
  rbind(cbind(rep = i, candidate = "package solution", qhet_conv_check(s$r, s$pc, seed = i)),
        cbind(rep = i, candidate = "package + 0.05", qhet_conv_check(s$r, s$pc, b_candidate = pkg + 0.05, seed = i)))
}), fill = TRUE)
t2[, expected := fifelse(candidate == "package + 0.05", "not usable", "usable")]
logl("T2 all as expected: {all(t2$status == t2$expected)}")

# T4
t4 <- rbindlist(lapply(1:3, function(i) {
  s <- sim(400, 0.0068)
  cbind(rep = i, qhet_conv_check(s$r, s$pc, seed = i, tau2_override = -2e-4))
}))
t4[, expected := "not usable"]
logl("T4 all as expected: {all(t4$status == t4$expected)}")

# T3
s100 <- sim(100, 0.0068)
set.seed(777); pk <- suppressWarnings(MVMR::qhet_mvmr(s100$r, pcor = s100$pc, CI = TRUE, iterations = 100, ncores = 1))
a <- qhet_boot_checked(s100$r, s100$pc, R = 100, ncores = 1, seed = 777)$summary
pk_ci <- do.call(rbind, lapply(strsplit(pk[["95% CI"]], "(?<=[0-9])-", perl = TRUE), as.numeric))
a[, `:=`(test = "a_no_failure_vs_package", pkg_lo = pk_ci[, 1], pkg_hi = pk_ci[, 2],
         matches_pkg_3dp = round(bca_lo, 3) == pk_ci[, 1] & round(bca_hi, 3) == pk_ci[, 2])]
b3 <- qhet_boot_checked(s100$r, s100$pc, R = 200, ncores = NC, seed = 778, fail_hook = function(i) runif(1) < 0.03)$summary
b3[, test := "b_forced_3pct"]
b10 <- qhet_boot_checked(s100$r, s100$pc, R = 200, ncores = NC, seed = 779, fail_hook = function(i) runif(1) < 0.10)$summary
b10[, test := "c_forced_10pct"]
t3 <- rbind(a, b3, b10, fill = TRUE)
logl("T3:\n{paste(capture.output(print(t3)), collapse = '\n')}")

write_result_table(t1, "qhet_checks_test_convergence")
write_result_table(t2, "qhet_checks_test_flag")
write_result_table(t4, "qhet_checks_test_nonpositive_w")
write_result_table(t3, "qhet_checks_test_bootstrap")
logl("Done.")
close_run_log(log_con)
