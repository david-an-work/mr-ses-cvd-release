#' 02f_qhet_validation.R
#'
#' Purpose: behaviour of the installed MVMR::qhet_mvmr (0.4.8; it warns that it is "currently undergoing
#'   development") versus IVW-MVMR on SIMULATED two-exposure data. This is a check before plan v0.3 makes qhet_mvmr
#'   part of Q2. No project data are used, and no project estimate is computed.
#' Design (fixed before running): L = 400 independent SNPs; true SNP-exposure effects g1 ~ N(0, 0.02^2),
#'   g2 = 0.9 g1 + N(0, tau^2); exposure SEs 0.005 with error correlation rho = 0.23 (as in our EA-income data);
#'   outcome effects by = 0.3 g1 - 0.1 g2 + N(0, 0.01^2) with no exposure-outcome overlap and no pleiotropy.
#'   tau is set to give a "weak" regime (conditional F about 2-3) and a "strong" regime (about 20).
#'   Calibration note (2026-10-06): the first run used tau = 0.0035 / 0.016, which gave median conditional F of
#'   1.4 / 7.5-8.6 (off target), and the CI runs failed (ncores = 4 exceeded the cores available). tau was
#'   recalibrated on exposure strength only (F - 1 is proportional to tau^2) to 0.0068 / 0.026, and ncores set to 1.
#'   First-run summary: weak IVW 0.157 / 0.050, qhet median 0.287 / -0.081 (SD 0.92 / 1.00); strong IVW 0.264 / -0.076,
#'   qhet median 0.303 / -0.099 (SD 0.04).
#'   Re-run outcome (2026-10-06, tau 0.0068 / 0.026): median conditional F 2.41 / 2.39 (weak) and 12.2 / 21.1 (strong;
#'   exposure 1 below the ~20 target). Weak: IVW mean 0.210 / -0.012; qhet median 0.295 / -0.089 (SD 0.10 / 0.11).
#'   Strong: IVW 0.276 / -0.090; qhet 0.300 / -0.100 (SD 0.03 / 0.02). 0 errors, 200/200 finite per regime.
#'   CI runs (200 iterations, 1 core) completed in 1925 s / 1744 s; weak 95% CI 0.161-0.534 and -0.309-0.034.
#'   R = 100 replicates per regime, point estimates only (CI = FALSE). Recorded per replicate: conditional F,
#'   IVW-MVMR and qhet estimates, warnings other than the fixed development/CI/iteration notices, errors, and
#'   whether qhet returned finite values. One bootstrap-CI run (200 iterations) per regime records run time and
#'   failures.
#' Outputs: results/tables/qhet_validation_{replicates,summary}.{tsv,xlsx}; results/logs/{date}_02f_qhet_validation.log
#' Depends on: MVMR, data.table, MASS

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(data.table))
ensure_pkgs_loaded(c("MVMR", "data.table", "MASS"))
log_con <- log_run("02f_qhet_validation")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
set.seed(20261006)

L <- 400; RHO <- 0.23; SEX <- 0.005; SEY <- 0.01; THETA <- c(0.3, -0.1); R <- 100
PC <- matrix(c(1, RHO, RHO, 1), 2)
fixed_warn <- "undergoing development|confidence interval not calculated|Iterations for bootstrap|Multi-core"

sim <- function(tau) {
  g1 <- rnorm(L, 0, 0.02); g2 <- 0.9 * g1 + rnorm(L, 0, tau)
  E <- MASS::mvrnorm(L, c(0, 0), PC * SEX^2)
  bx <- cbind(g1, g2) + E; sx <- matrix(SEX, L, 2)
  by <- as.vector(cbind(g1, g2) %*% THETA) + rnorm(L, 0, SEY)
  MVMR::format_mvmr(BXGs = bx, BYG = by, seBXGs = sx, seBYG = rep(SEY, L), RSID = paste0("s", seq_len(L)))
}
run_one <- function(r_in) {
  w <- character(0)
  F_ <- withCallingHandlers(suppressMessages(capture.output(f <- MVMR::strength_mvmr(r_in, gencov = MVMR::phenocov_mvmr(PC, as.matrix(r_in[, c("sebetaX1", "sebetaX2")]))))),
                            warning = function(x) invokeRestart("muffleWarning"))
  ivw <- suppressWarnings(MVMR::ivw_mvmr(r_in, gencov = MVMR::phenocov_mvmr(PC, as.matrix(r_in[, c("sebetaX1", "sebetaX2")]))))
  q <- tryCatch(withCallingHandlers(MVMR::qhet_mvmr(r_in, pcor = PC, CI = FALSE, iterations = 0),
                  warning = function(x) { if (!grepl(fixed_warn, conditionMessage(x))) w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }),
                error = function(e) structure(NA, err = conditionMessage(e)))
  qv <- if (is.data.frame(q) || is.matrix(q)) as.numeric(unlist(q[, 1])) else c(NA_real_, NA_real_)
  data.table(F1 = as.numeric(f[1, 1]), F2 = as.numeric(f[1, 2]), ivw1 = ivw[1, 1], ivw2 = ivw[2, 1],
             q1 = qv[1], q2 = qv[2], q_error = if (is.null(attr(q, "err"))) NA_character_ else attr(q, "err"),
             q_warnings = paste(unique(w), collapse = " | "))
}
regimes <- c(weak = 0.0068, strong = 0.026)
reps <- rbindlist(lapply(names(regimes), function(nm) {
  logl("regime {nm} (tau {regimes[[nm]]})")
  rbindlist(lapply(seq_len(R), function(i) cbind(regime = nm, rep = i, run_one(sim(regimes[[nm]])))))
}))
summ <- reps[, .(medF1 = median(F1), medF2 = median(F2),
                 ivw1_mean = mean(ivw1), ivw1_sd = sd(ivw1), ivw2_mean = mean(ivw2), ivw2_sd = sd(ivw2),
                 q_finite = sum(is.finite(q1) & is.finite(q2)), q_errors = sum(!is.na(q_error)),
                 q_extra_warnings = sum(nzchar(q_warnings)),
                 q1_mean = mean(q1, na.rm = TRUE), q1_median = median(q1, na.rm = TRUE), q1_sd = sd(q1, na.rm = TRUE),
                 q2_mean = mean(q2, na.rm = TRUE), q2_median = median(q2, na.rm = TRUE), q2_sd = sd(q2, na.rm = TRUE)),
             by = regime]
summ[, `:=`(true1 = THETA[1], true2 = THETA[2])]
# one bootstrap-CI run per regime: run time and outcome
ci <- rbindlist(lapply(names(regimes), function(nm) {
  r_in <- sim(regimes[[nm]]); t0 <- Sys.time()
  out <- tryCatch(suppressWarnings(MVMR::qhet_mvmr(r_in, pcor = PC, CI = TRUE, iterations = 200, ncores = 1)), error = function(e) conditionMessage(e))
  data.table(regime = nm, seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
             result = paste(capture.output(print(out)), collapse = " / "))
}))
write_result_table(reps, "qhet_validation_replicates")
write_result_table(cbind(summ), "qhet_validation_summary")
write_result_table(ci, "qhet_validation_ci_runs")
logl("summary:\n{paste(capture.output(print(summ)), collapse = '\n')}")
logl("CI runs:\n{paste(capture.output(print(ci)), collapse = '\n')}")
logl("Done.")
close_run_log(log_con)
