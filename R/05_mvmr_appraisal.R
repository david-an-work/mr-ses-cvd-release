#' 05_mvmr_appraisal.R
#'
#' Purpose: appraisal component (plan v1.0 §10; Q1, Q2) on the outcome-specific harmonized MVMR sets (R/03).
#'   For every model (§10.2) x outcome file:
#'   - Sanderson-Windmeijer conditional F at 0 / 0.5 / 1 / 1.5 x the LDSC-intercept correlation (§10.1), computed as in
#'     R/02b (MVMR::strength_mvmr + phenocov_mvmr) on the SNPs retained for that outcome, and re-computed independently
#'     with the R/02c formula; the script stops on a relative difference > 1e-8.
#'   - IVW-MVMR, multiplicative random effects (MendelianRandomization::mr_mvivw, model = "random", normal CIs;
#'     residual SE floored at 1), cross-checked against a closed-form implementation (stops on difference > 1e-10).
#'   - The conventional fixed-effect SE sqrt(diag((B'WB)^-1)) as a precision benchmark (§11).
#'   - Label (§10.4): a coefficient with conditional F < 10 at any of the four scales is "weak-instrument-limited".
#'   No other MVMR estimator is used (§10.3).
#' Models: project A1 four-SES, A2 EA+Kweon, A3 EA+TDI, A4 EA+CAMSIS, A5 four-SES with Hill; EUR-only TDI sensitivity
#'   versions of A1 and A3; comparable B1 Li D (407 SNPs; B1-pruned LD-pruned sensitivity), B2/B3 Le single-adjustment
#'   (EA+Hill, which uses the same SNP set as B1, and EA+ISEI). Ye 2025 (B3) is not reconstructed (files not obtained).
#' Outcomes: all 12 harmonized outcome files; the role column marks primary, comparator and exploratory outcomes.
#' Outputs: results/tables/mvmr_appraisal.tsv (one row per model x outcome x exposure),
#'   results/tables/mvmr_condF_verification.tsv
#' Depends on: data.table, MVMR, MendelianRandomization

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
log_con <- log_run("05_mvmr_appraisal")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
H <- file.path(DIR_DATA_HARM, "large", "mvmr"); DIR_CONDF <- file.path(DIR_TABLES, "condF_inputs")
SCALES <- c(0, 0.5, 1, 1.5)

models <- data.table(
  model = c("EA+Income+TDI+Occupation", "EA+Income", "EA+TDI", "EA+Occupation", "EA+Income_Hill+TDI+Occupation",
            "EA+TDI_EUR", "EA+Income+TDI_EUR+Occupation",
            "LiD_S1", "LiD_S1-LDpruned", "EA+Income_Hill", "EA+Occupation_ISEI"),
  label = c("A1 four-measure (EA4, Kweon income, TDI, CAMSIS)", "A2 EA4 + Kweon income", "A3 EA4 + deprivation (TDI)", "A4 EA4 + occupation (CAMSIS)",
            "A5 four-measure with Hill income", "A3-EUR EA4 + EUR-only TDI", "A1-EUR A1 with EUR-only TDI",
            "B1 Li D-comparable (EA4 + Hill income)", "B1-pruned B1 with LD-pruned instrument set",
            "B2 Le-comparable (EA4 + Hill income; same set as B1)", "B3 Le-comparable (EA4 + occupation, ISEI)"),
  pcor = c("EA+Income+TDI+Occupation", "EA+Income", "EA+TDI", "EA+Occupation", "EA+Income_Hill+TDI+Occupation",
           "EA+TDI_EUR", "EA+Income+TDI_EUR+Occupation", "EA+Income_Hill", "EA+Income_Hill", "EA+Income_Hill",
           "EA+Occupation_ISEI"),
  group = c(rep("project", 5), rep("project (sensitivity)", 2), "comparable", "comparable (sensitivity)", "comparable", "comparable"))
outcomes <- c(cad_aragam_eur = "primary", cad_aragam_multi = "comparator", cad_mvp_phe411_eur = "comparator (no UKB)",
              pad_mvp_eur = "primary", pad_finngen_athscle = "comparator", pad_vanzuydam = "exploratory (rounded data)",
              pad_mvp_afr = "exploratory (transportability)", hf_hermes2 = "primary", hf_mvp_phe428_eur = "comparator (no UKB)",
              stroke_gigastroke_multi = "primary", stroke_gigastroke_eur = "comparator",
              stroke_mvp_phe433_21_eur = "comparator (no UKB)")

# independent conditional F (R/02c formula)
cond_F_indep <- function(bx, sx, pcor, scale) {
  K <- ncol(bx); L <- nrow(bx); pc <- pcor * scale; diag(pc) <- 1
  vapply(seq_len(K), function(j) {
    X <- bx[, -j, drop = FALSE]; y <- bx[, j]
    delta <- solve(crossprod(X), crossprod(X, y)); d <- numeric(K); d[j] <- -1; d[-j] <- delta
    ds <- sweep(sx, 2, d, `*`); sigma2 <- rowSums((ds %*% pc) * ds)
    sum((y - X %*% delta)^2 / sigma2) / L
  }, 0)
}
# closed-form IVW-MVMR (first-order weights, multiplicative RE floored at 1)
ivw_mv <- function(bx, by, sy) {
  W <- diag(1 / sy^2); A <- solve(t(bx) %*% W %*% bx); th <- A %*% t(bx) %*% W %*% by
  r <- by - bx %*% th; L <- nrow(bx); K <- ncol(bx)
  phi <- as.numeric(t(r) %*% W %*% r) / (L - K)
  list(b = as.numeric(th), se_fe = sqrt(diag(A)), se = sqrt(diag(A)) * max(1, sqrt(phi)), phi = phi,
       Q = as.numeric(t(r) %*% W %*% r), Q_df = L - K)
}

rows <- list(); ver <- list()
for (i in seq_len(nrow(models))) for (o in names(outcomes)) {
  m <- models[i]
  h <- fread(file.path(H, paste0(m$model, "__", o, ".tsv.gz")))
  d <- h[mr_keep == TRUE]
  bc <- grep("^b_", names(d), value = TRUE); sc <- sub("^b_", "se_", bc)
  bx <- as.matrix(d[, ..bc]); sx <- as.matrix(d[, ..sc]); by <- d$beta_y_aligned; sy <- d$se_y
  pcf <- fread(file.path(DIR_CONDF, paste0(m$pcor, "_pcor.tsv")))
  pcor <- as.matrix(pcf[, -1]); dimnames(pcor) <- list(pcf$trait, pcf$trait)
  stopifnot(ncol(pcor) == length(bc))
  exps <- pcf$trait
  # conditional F, package route (as R/02b) and independent route (as R/02c)
  r_input <- MVMR::format_mvmr(BXGs = bx, BYG = rep(0, nrow(d)), seBXGs = sx, seBYG = rep(1, nrow(d)), RSID = d$rsid)
  F <- sapply(SCALES, function(s) {
    pc <- pcor * s; diag(pc) <- 1
    utils::capture.output(f <- suppressWarnings(if (s == 0) MVMR::strength_mvmr(r_input, gencov = 0)
                                                 else MVMR::strength_mvmr(r_input, gencov = MVMR::phenocov_mvmr(pc, sx))))
    as.numeric(f[1, ])
  })
  Fi <- sapply(SCALES, function(s) cond_F_indep(bx, sx, pcor, s))
  F <- matrix(F, nrow = length(exps)); Fi <- matrix(Fi, nrow = length(exps))
  rel <- max(abs(F - Fi) / pmax(1, abs(Fi)))
  ver[[length(ver) + 1]] <- data.table(model = m$model, outcome = o, n_snps = nrow(d), max_rel_diff = rel)
  if (rel > 1e-8) stop("conditional F verification failed: ", m$model, " x ", o, " (", rel, ")")
  # IVW-MVMR
  mvi <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = bx, bxse = sx, by = by, byse = sy, snps = d$rsid),
                                         model = "random", distribution = "normal")
  cf <- ivw_mv(bx, by, sy)
  if (max(abs(mvi@Estimate - cf$b)) > 1e-10 || max(abs(mvi@StdError - cf$se)) > 1e-10) stop("IVW-MVMR cross-check failed: ", m$model, " x ", o)
  Fmin <- apply(F, 1, min); Fmax <- apply(F, 1, max)
  rows[[length(rows) + 1]] <- data.table(model = m$model, model_label = m$label, group = m$group, outcome = o,
    outcome_role = outcomes[[o]], n_snps = nrow(d), exposure = exps,
    condF_x0 = F[, 1], condF_x0.5 = F[, 2], condF_x1 = F[, 3], condF_x1.5 = F[, 4], condF_min = Fmin, condF_max = Fmax,
    weak_instrument_limited = Fmin < 10,
    b = mvi@Estimate, se = mvi@StdError, lo = mvi@CILower, hi = mvi@CIUpper, p = mvi@Pvalue, se_conventional = cf$se_fe,
    or = exp(mvi@Estimate), or_lo = exp(mvi@CILower), or_hi = exp(mvi@CIUpper),
    residual_se = sqrt(cf$phi), Q = cf$Q, Q_df = cf$Q_df)
}
R <- rbindlist(rows); Vt <- rbindlist(ver)
fwrite(R, file.path(DIR_TABLES, "mvmr_appraisal.tsv"), sep = "\t")
fwrite(Vt, file.path(DIR_TABLES, "mvmr_condF_verification.tsv"), sep = "\t")
logl("models x outcomes: {nrow(Vt)}; coefficient rows: {nrow(R)}; conditional F verification max rel diff {signif(max(Vt$max_rel_diff), 3)}")
logl("coefficients with conditional F >= 10 at all four scales: {sum(!R$weak_instrument_limited)} of {nrow(R)}")
close_run_log(log_con)
