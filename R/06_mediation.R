#' 06_mediation.R
#'
#' Purpose: exploratory education mediation (plan v1.0 §9; H2). EA4 only; seven single-mediator models (BMI, smoking,
#'   LDL, HDL, TG, SBP, T2D) x the four primary outcomes (CAD Aragam EUR, PAD MVP EUR, HF HERMES 2.0, any stroke
#'   GIGASTROKE multi). Run regardless of the univariable results.
#'   - Strength gate (§9.3): conditional F for EA and the mediator on the outcome-specific harmonized set, at 0 / 0.5 /
#'     1 / 1.5 x the LDSC-intercept correlation (pcor from R/02d); eligible only if both are >= 10 at all four scales.
#'   - Main (§9.2): product of coefficients a x b, a = IVW EA -> M (EA4 instruments; R/03 mediation cells), b = mediator
#'     coefficient of the single-mediator IVW-MVMR (multiplicative random effects). Secondary: difference method,
#'     total (univariable IVW EA4 -> Y) minus direct (EA coefficient of the MVMR). Proportions only if the total effect's
#'     95% CI excludes 0; never summed across mediators.
#'   - Main CI (§9.4): joint SNP bootstrap, 2,000 resamples, percentile 95% intervals. Pool = union of the SNPs of the
#'     three components for this outcome (MVMR set, EA4 -> Y set, EA4 -> M set); each resample draws from the pool with
#'     replacement and every component is re-estimated on the instruments of its own type in the resample (with
#'     multiplicity). Failure: any component with < 3 instruments, an error or a non-finite value. Any failure ->
#'     "failure-affected"; > 5% failures -> interval withheld. Delta-method (independence) intervals for comparison.
#'   - Labels: weak-instrument-limited (gate failed); mediator-outcome overlap for the MVP outcome (§9.5).
#'   Intervals do not account for participant overlap between the GWAS (§9.4).
#' Outputs: results/tables/mediation_education.tsv ; results/tables/mediation_education_condF.tsv
#' Depends on: data.table, MVMR, MendelianRandomization

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
log_con <- log_run("06_mediation")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
HL <- file.path(DIR_DATA_HARM, "large"); DIR_PC <- file.path(DIR_TABLES, "condF_inputs_mediators")
N_BOOT <- as.integer(Sys.getenv("MED_BOOT", "2000")); SCALES <- c(0, 0.5, 1, 1.5)
mediators <- c("BMI", "Smoking", "LDL", "HDL", "TG", "SBP", "T2D")
outcomes <- c("cad_aragam_eur", "pad_mvp_eur", "hf_hermes2", "stroke_gigastroke_multi")
overlap_label <- function(m, o) {
  if (o != "pad_mvp_eur") return(NA_character_)
  if (m %in% c("SBP", "LDL", "HDL", "TG")) return("mediator-outcome overlap: MVP")
  if (m == "T2D") return("mediator-outcome overlap: MVP possible (no shared cohort found)")
  NA_character_
}

ivw <- function(bx, by, sy) {  # IVW, first-order weights; returns estimate and multiplicative-RE SE
  w <- 1 / sy^2; b <- sum(w * bx * by) / sum(w * bx^2); L <- length(bx)
  phi <- if (L > 1) sum(w * (by - b * bx)^2) / (L - 1) else 1
  c(b = b, se = sqrt(1 / sum(w * bx^2)) * max(1, sqrt(phi)))
}
mvivw <- function(X, by, sy) {  # IVW-MVMR coefficients and multiplicative-RE SEs
  W <- 1 / sy^2; A <- solve(crossprod(X * sqrt(W))); th <- A %*% crossprod(X, W * by)
  r <- by - X %*% th; phi <- sum(W * r^2) / (nrow(X) - ncol(X))
  list(b = as.numeric(th), se = sqrt(diag(A)) * max(1, sqrt(phi)))
}
cond_F <- function(bx, sx, pcor, scale) {  # Sanderson-Windmeijer (R/02c formula; = MVMR::strength_mvmr)
  K <- ncol(bx); L <- nrow(bx); pc <- pcor * scale; diag(pc) <- 1
  vapply(seq_len(K), function(j) {
    X <- bx[, -j, drop = FALSE]; y <- bx[, j]
    delta <- solve(crossprod(X), crossprod(X, y)); d <- numeric(K); d[j] <- -1; d[-j] <- delta
    ds <- sweep(sx, 2, d, `*`); sum((y - X %*% delta)^2 / rowSums((ds %*% pc) * ds)) / L
  }, 0)
}
pct <- function(x) if (all(is.na(x))) c(NA_real_, NA_real_) else quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

grid <- CJ(m = mediators, o = outcomes, sorted = FALSE); grid[, seed := PROJECT_SEED + .I]
rows <- list(); frows <- list()
for (i in seq_len(nrow(grid))) {
  M <- grid$m[i]; O <- grid$o[i]; set.seed(grid$seed[i])
  mv <- fread(file.path(HL, "mvmr", paste0("EA+", M, "__", O, ".tsv.gz")))[mr_keep == TRUE]
  bc <- grep("^b_", names(mv), value = TRUE); stopifnot(length(bc) == 2, bc[1] == "b_ea_okbay2022")
  X <- as.matrix(mv[, ..bc]); SX <- as.matrix(mv[, sub("^b_", "se_", bc), with = FALSE])
  am <- fread(file.path(HL, "mediation", paste0("EA4__", M, ".tsv.gz")))[mr_keep == TRUE]
  tt <- fread(file.path(HL, "univariable", paste0("EA4__", O, ".tsv.gz")))[mr_keep == TRUE]
  # strength gate on the outcome-specific MVMR set
  pcf <- fread(file.path(DIR_PC, paste0("EA+", M, "_pcor.tsv"))); pcor <- as.matrix(pcf[, -1])
  Fm <- sapply(SCALES, function(s) cond_F(X, SX, pcor, s))
  # cross-check one scale against the package route (as R/02b/02d)
  ri <- MVMR::format_mvmr(BXGs = X, BYG = rep(0, nrow(X)), seBXGs = SX, seBYG = rep(1, nrow(X)), RSID = mv$rsid)
  utils::capture.output(fp <- suppressWarnings(MVMR::strength_mvmr(ri, gencov = MVMR::phenocov_mvmr(pcor, SX))))
  if (max(abs(as.numeric(fp[1, ]) - Fm[, 3]) / Fm[, 3]) > 1e-8) stop("conditional F cross-check failed: ", M, " x ", O)
  gate <- all(Fm >= 10)
  frows[[i]] <- data.table(mediator = M, outcome = O, n_snps_mvmr = nrow(X), exposure = c("EA", M),
                           condF_x0 = Fm[, 1], condF_x0.5 = Fm[, 2], condF_x1 = Fm[, 3], condF_x1.5 = Fm[, 4])
  # point estimates (closed form, cross-checked against MendelianRandomization)
  a <- ivw(am$beta_x, am$beta_y_aligned, am$se_y); tot <- ivw(tt$beta_x, tt$beta_y_aligned, tt$se_y)
  mvf <- mvivw(X, mv$beta_y_aligned, mv$se_y)
  chk <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = X, bxse = SX, by = mv$beta_y_aligned, byse = mv$se_y),
                                          model = "random", distribution = "normal")
  chk2 <- MendelianRandomization::mr_ivw(MendelianRandomization::mr_input(bx = am$beta_x, bxse = am$se_x, by = am$beta_y_aligned, byse = am$se_y),
                                         model = "random", distribution = "normal")
  if (max(abs(chk@Estimate - mvf$b), abs(chk@StdError - mvf$se), abs(chk2@Estimate - a[["b"]]), abs(chk2@StdError - a[["se"]])) > 1e-10)
    stop("point-estimate cross-check failed: ", M, " x ", O)
  direct <- mvf$b[1]; bM <- mvf$b[2]; ind <- a[["b"]] * bM; dif <- tot[["b"]] - direct
  # joint SNP bootstrap
  pool <- unique(c(mv$rsid, tt$rsid, am$rsid))
  setkey(mv, rsid); setkey(tt, rsid); setkey(am, rsid)
  B <- matrix(NA_real_, N_BOOT, 4, dimnames = list(NULL, c("indirect", "difference", "total", "direct")))
  for (r in seq_len(N_BOOT)) {
    s <- sample(pool, length(pool), replace = TRUE)
    ok <- tryCatch({
      i1 <- mv[J(s), nomatch = NULL]; i2 <- tt[J(s), nomatch = NULL]; i3 <- am[J(s), nomatch = NULL]
      if (nrow(i1) < 3 || nrow(i2) < 3 || nrow(i3) < 3) stop("fewer than 3 instruments")
      fb <- mvivw(as.matrix(i1[, ..bc]), i1$beta_y_aligned, i1$se_y)
      ab <- ivw(i3$beta_x, i3$beta_y_aligned, i3$se_y)[["b"]]; tb <- ivw(i2$beta_x, i2$beta_y_aligned, i2$se_y)[["b"]]
      v <- c(ab * fb$b[2], tb - fb$b[1], tb, fb$b[1]); if (!all(is.finite(v))) stop("non-finite"); v
    }, error = function(e) NULL)
    if (!is.null(ok)) B[r, ] <- ok
  }
  nfail <- sum(!complete.cases(B)); ffrac <- nfail / N_BOOT
  ci_label <- if (nfail == 0) "joint SNP bootstrap" else if (ffrac <= 0.05) "joint SNP bootstrap (failure-affected, exploratory)" else "withheld (> 5% failed resamples)"
  withheld <- ffrac > 0.05
  ci_ind <- if (withheld) c(NA, NA) else pct(B[, "indirect"]); ci_dif <- if (withheld) c(NA, NA) else pct(B[, "difference"])
  tot_ci <- tot[["b"]] + c(-1, 1) * qnorm(0.975) * tot[["se"]]; prop_ok <- tot_ci[1] > 0 || tot_ci[2] < 0
  pr_ind <- if (prop_ok && !withheld) pct(B[, "indirect"] / B[, "total"]) else c(NA, NA)
  pr_dif <- if (prop_ok && !withheld) pct(B[, "difference"] / B[, "total"]) else c(NA, NA)
  se_d_ind <- sqrt(a[["b"]]^2 * mvf$se[2]^2 + bM^2 * a[["se"]]^2); se_d_dif <- sqrt(tot[["se"]]^2 + mvf$se[1]^2)
  rows[[i]] <- data.table(mediator = M, outcome = O, seed = grid$seed[i],
    n_snps_mvmr = nrow(X), n_snps_total = nrow(tt), n_snps_ea_to_m = nrow(am),
    condF_EA_min = min(Fm[1, ]), condF_EA_max = max(Fm[1, ]), condF_M_min = min(Fm[2, ]), condF_M_max = max(Fm[2, ]),
    gate_passed = gate, label = if (gate) "eligible for interpretation (exploratory)" else "weak-instrument-limited",
    overlap_label = overlap_label(M, O),
    a_ea_to_m = a[["b"]], a_se = a[["se"]], b_m_to_y = bM, b_se = mvf$se[2],
    total = tot[["b"]], total_se = tot[["se"]], total_lo = tot_ci[1], total_hi = tot_ci[2],
    direct = direct, direct_se = mvf$se[1],
    indirect = ind, indirect_lo = ci_ind[1], indirect_hi = ci_ind[2],
    difference = dif, difference_lo = ci_dif[1], difference_hi = ci_dif[2],
    prop_indirect = if (prop_ok) ind / tot[["b"]] else NA_real_, prop_indirect_lo = pr_ind[1], prop_indirect_hi = pr_ind[2],
    prop_difference = if (prop_ok) dif / tot[["b"]] else NA_real_, prop_difference_lo = pr_dif[1], prop_difference_hi = pr_dif[2],
    indirect_delta_lo = ind - qnorm(0.975) * se_d_ind, indirect_delta_hi = ind + qnorm(0.975) * se_d_ind,
    difference_delta_lo = dif - qnorm(0.975) * se_d_dif, difference_delta_hi = dif + qnorm(0.975) * se_d_dif,
    n_boot = N_BOOT, n_failed = nfail, ci_label = ci_label)
  logl("{M} x {O}: gate {gate}; indirect {signif(ind, 3)}; failures {nfail}")
}
R <- rbindlist(rows); Fm <- rbindlist(frows)
fwrite(R, file.path(DIR_TABLES, "mediation_education.tsv"), sep = "\t")
fwrite(Fm, file.path(DIR_TABLES, "mediation_education_condF.tsv"), sep = "\t")
logl("models: {nrow(R)}; gate passed: {sum(R$gate_passed)}; any failed resamples: {sum(R$n_failed > 0)}")
close_run_log(log_con)
