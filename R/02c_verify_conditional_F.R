#' 02c_verify_conditional_F.R
#'
#' Purpose: independent re-computation of the conditional F-statistics reported by
#'   R/02b_ses_feasibility.R (TASKS.md T2.12), on the identical inputs it saved.
#'   Base-R matrix algebra only; MVMR is not called.
#' Scope: this verifies the conditional-F ARITHMETIC given the saved inputs: the MVMR
#'   call, the formatting of B/SE into it, and the per-SNP covariance built from the saved
#'   correlation matrix and its scaling. It does NOT verify anything upstream of the saved
#'   files: SNP selection and joint clumping, allele alignment, palindrome handling, or the
#'   conversion of LDSC cross-trait intercepts into the error-correlation matrix. Both
#'   calculations share those inputs, so an error there would be reproduced, not caught.
#' Failure: the script stops (non-zero exit) if any recomputed value differs from the
#'   unrounded package value by more than TOL, or from the rounded table by more than 0.005.
#'
#' Method (Sanderson et al. 2021, Stat Med 40:2541): for exposure j with SNP-exposure
#'   estimates y = B[, j] and X = B[, -j], delta = unweighted OLS of y on X (no intercept);
#'   for each SNP l, sigma2_l = d' Sigma_l d with d = (-1 at j, delta elsewhere) and
#'   Sigma_l = pcor_scaled * (se_l se_l'); Q = sum_l (y_l - X_l delta)^2 / sigma2_l.
#'   F_pkg = Q / L matches MVMR::strength_mvmr (0.4.8); F_paper = Q / (L - (K - 1)) is
#'   the paper's degrees-of-freedom version. Both are reported.
#'
#' Inputs:  results/tables/condF_inputs/{model}_B.tsv, {model}_pcor.tsv, {model}_condF_unrounded.tsv
#'          results/tables/feasibility_conditional_F.tsv
#' Outputs: results/tables/feasibility_conditional_F_verification.{tsv,xlsx}
#'          results/logs/{date}_02c_verify_conditional_F.log (+ sessionInfo)
#' Depends on: data.table

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(data.table))

log_con <- log_run("02c_verify_conditional_F")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)

DIR_CONDF <- file.path(DIR_TABLES, "condF_inputs")
TOL <- 1e-8  # relative tolerance against the unrounded package values
pkg <- fread(file.path(DIR_TABLES, "feasibility_conditional_F.tsv"))
scale_cols <- grep("^condF_", names(pkg), value = TRUE)
col_scale <- c(condF_independent = 0, condF_overlap_adjusted = 1)
scale_of <- function(cn) if (cn %in% names(col_scale)) col_scale[[cn]] else as.numeric(sub("condF_pcor_x", "", cn))

cond_F <- function(B, pcor, scale) {
  keys <- sub("^b_", "", grep("^b_", names(B), value = TRUE))
  bx <- as.matrix(B[, paste0("b_", keys), with = FALSE])
  sx <- as.matrix(B[, paste0("se_", keys), with = FALSE])
  L <- nrow(bx); K <- ncol(bx)
  pc <- pcor * scale; diag(pc) <- 1
  out <- matrix(NA_real_, K, 2, dimnames = list(NULL, c("F_pkg", "F_paper")))
  for (j in seq_len(K)) {
    X <- bx[, -j, drop = FALSE]; y <- bx[, j]
    delta <- solve(crossprod(X), crossprod(X, y))
    d <- numeric(K); d[j] <- -1; d[-j] <- delta
    # sigma2_l = sum_ab d_a d_b pc_ab se_la se_lb, vectorised over SNPs
    ds <- sweep(sx, 2, d, `*`)
    sigma2 <- rowSums((ds %*% pc) * ds)
    Q <- sum((y - X %*% delta)^2 / sigma2)
    out[j, ] <- c(Q / L, Q / (L - (K - 1)))
  }
  out
}

rows <- list()
for (mdl in unique(pkg$model)) {
  model_id <- gsub(" \\+ ", "+", mdl)
  B <- fread(file.path(DIR_CONDF, paste0(model_id, "_B.tsv")))
  pcor_dt <- fread(file.path(DIR_CONDF, paste0(model_id, "_pcor.tsv")))
  pcor <- as.matrix(pcor_dt[, -1]); dimnames(pcor) <- list(pcor_dt$trait, pcor_dt$trait)
  p <- pkg[model == mdl]
  unr <- fread(file.path(DIR_CONDF, paste0(model_id, "_condF_unrounded.tsv")))
  stopifnot(identical(p$exposure, colnames(pcor)), nrow(B) == p$n_snps[1], identical(unr$exposure, p$exposure))
  for (cn in scale_cols) {
    sc <- scale_of(cn)
    f <- cond_F(B, pcor, sc)
    rows[[length(rows) + 1]] <- data.table(model = mdl, exposure = p$exposure, n_snps = nrow(B),
      pcor_scale = sc, condF_package_rounded = p[[cn]], condF_package_unrounded = unr[[paste0("condF_pcor_x", sc)]],
      condF_recomputed = f[, "F_pkg"], condF_recomputed_df_corrected = f[, "F_paper"])
  }
}
ver <- rbindlist(rows)
ver[, rel_diff_unrounded := abs(condF_package_unrounded - condF_recomputed) / abs(condF_recomputed)]
ver[, abs_diff_vs_rounded := abs(condF_package_rounded - condF_recomputed)]
ver[, pass := !is.na(rel_diff_unrounded) & rel_diff_unrounded <= TOL & abs_diff_vs_rounded <= 0.005 + 1e-12]
write_result_table(ver, "feasibility_conditional_F_verification")
logl("rows: {nrow(ver)}; pass: {sum(ver$pass)}; max relative diff vs unrounded: {signif(max(ver$rel_diff_unrounded), 3)}; max abs diff vs 2-dp table: {signif(max(ver$abs_diff_vs_rounded), 3)}")
logl("max conditional F at any covariance scale: {signif(max(ver$condF_recomputed), 5)} (Q/L); {signif(max(ver$condF_recomputed_df_corrected), 5)} (Q/(L-K+1))")
if (!all(ver$pass)) stop(sum(!ver$pass), " of ", nrow(ver), " conditional-F values failed verification; see feasibility_conditional_F_verification.tsv")
logl("Done.")
close_run_log(log_con)
