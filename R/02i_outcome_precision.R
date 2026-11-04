#' 02i_outcome_precision.R
#'
#' Purpose: comparison of candidate headline outcomes (CAD vs PAD, with HF and stroke for reference) before
#'   registration, without any causal estimate. For each outcome, reads the outcome standard errors (and alleles/positions
#'   for matching) at the instrument SNPs. Outcome betas and p-values are not read. Exception: for van Zuydam the
#'   deposited OR is read, only to convert or_se to SE(log OR) = or_se / OR. No project causal estimate is calculated,
#'   and none was used to choose the headline outcome.
#' Quantities:
#'   - Univariable expected SE of the IVW estimate (fixed effect, first order): 1 / sqrt(sum(bx^2 / seY^2)) over the
#'     exposure's instruments available in the outcome; minimum detectable OR per SD at 80% power,
#'     alpha = 0.05/16 (two-sided): exp(+/-(qnorm(1 - 0.003125/2) + qnorm(0.8)) * SE).
#'   - MVMR expected SE: sqrt(diag((B' W B)^-1)), W = diag(1/seY^2), B = exposure betas of the model's SNP set
#'     (results/tables/condF_inputs/*_B.tsv; Li D-comparable S1 from pilot_lid_cad_inputs), restricted to SNPs available
#'     in the outcome. A CONVENTIONAL PRECISION BENCHMARK: it treats exposure associations as known, so it does not
#'     reflect weak-instrument bias (whose direction is model-dependent with correlated exposures), overlap or pleiotropy.
#'     Useful only for comparing outcome precision.
#'   - R2 / bx come from the exposure discovery GWAS (winner's curse): detectable effects are optimistic.
#' Matching: rsID for MVP, van Zuydam, HERMES, FinnGen (rsids column); GRCh37 chr:pos + allele pair for Aragam and
#'   GIGASTROKE. van Zuydam SE = or_se / OR (provisional, R/02g). MVP SE from the 95% CI.
#' Outputs: results/tables/outcome_precision_{univariable,mvmr,outcomes}.{tsv,xlsx}; results/logs/{date}_02i_outcome_precision.log
#' Depends on: data.table

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
log_con <- log_run("02i_outcome_precision")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
RAW <- here::here("data", "raw")
# MVP case/control counts from the dbGaP metadata files (filled in when the files are downloaded)
MVP_N <- list(hf = c(58749, 370923), st = c(20766, 420219))  # GWAS Catalog GCST90475981 / GCST90477991 (EUR, GIA)

# ---- exposures ----
expo <- c(EA4 = "ea_okbay2022", Kweon = "income_kweon2025", TDI = "tdi_panukb", CAMSIS = "occupation_akimova2024",
          Hill = "income_hill2019", ISEI = "occupation_akimova2024_isei", gSES = "ses_gses_xia2025")
inst <- rbindlist(lapply(names(expo), function(e) {
  d <- fread(file.path(DIR_INSTRUMENTS, paste0(expo[[e]], "_instruments.tsv")))
  d[, .(exposure = e, rsid, chr = as.character(chr), pos, a1 = toupper(effect_allele), a2 = toupper(other_allele), bx = beta)]
}))
Bfiles <- c(`EA4+Kweon` = "EA+Income", `EA4+TDI` = "EA+TDI", `EA4+CAMSIS` = "EA+Occupation",
            `4-SES` = "EA+Income+TDI+Occupation", `4-SES (Hill)` = "EA+Income_Hill+TDI+Occupation")
# Le-comparable pairs (added 2026-10-06), used if R/02b has written their SNP sets
le_models <- c(`EA4+Hill (Le)` = "EA+Income_Hill", `EA4+ISEI (Le)` = "EA+Occupation_ISEI")
le_ok <- file.exists(file.path(DIR_TABLES, "condF_inputs", paste0(le_models, "_B.tsv")))
Bfiles <- c(Bfiles, le_models[le_ok])
Bs <- lapply(Bfiles, function(f) fread(file.path(DIR_TABLES, "condF_inputs", paste0(f, "_B.tsv"))))
lid <- list.files(file.path(DIR_TABLES, "pilot_lid_cad_inputs"), pattern = "^S1_B\\.tsv$", full.names = TRUE)
if (length(lid) == 1) Bs[["Li D-comparable S1 (EA4+Hill)"]] <- fread(lid)
all_rsid <- unique(c(inst$rsid, unlist(lapply(Bs, `[[`, "rsid"))))
# rsIDs of adopted exposures only: the outcome-table coverage column excludes the gSES candidate (added 2026-10-07)
core_rsid <- unique(c(inst[exposure != "gSES", rsid], unlist(lapply(Bs, `[[`, "rsid"))))
pos_map <- unique(inst[, .(rsid, chr, pos, a1, a2)], by = "rsid")
logl("instrument rsIDs: {length(all_rsid)}; with GRCh37 positions: {nrow(pos_map)}; MVMR models: {length(Bs)}")

# ---- outcomes: return data.table(rsid, seY) ----
by_pos <- function(f, chrc, posc, a1c, a2c, sec) {
  d <- fread(f, select = c(chrc, posc, a1c, a2c, sec), showProgress = FALSE)
  setnames(d, c("chr", "pos", "o1", "o2", "seY")); d[, chr := as.character(chr)]
  d <- d[pos_map, on = .(chr, pos), nomatch = 0]
  d <- d[(toupper(o1) == a1 & toupper(o2) == a2) | (toupper(o1) == a2 & toupper(o2) == a1)]
  unique(d[, .(rsid, seY)], by = "rsid")
}
by_rsid <- function(d) unique(d[rsid %in% all_rsid & is.finite(seY) & seY > 0, .(rsid, seY)], by = "rsid")

outcomes <- list(
  `CAD Aragam EUR` = list(cases = 181522, controls = 984168, ukb = "yes (18.7% of cases, 44.5% of controls)",
    get = function() by_pos(file.path(RAW, "cad_aragam2022_eur.tsv"), "chromosome", "base_pair_location", "effect_allele", "other_allele", "standard_error")),
  `CAD Aragam multi` = list(cases = 210842, controls = 1167328, ukb = "yes (16.1% / 37.6%)",
    get = function() by_pos(file.path(RAW, "cad_aragam2022_multi.tsv"), "chromosome", "base_pair_location", "effect_allele", "other_allele", "standard_error")),
  `PAD MVP EUR` = list(cases = 15036, controls = 423003, ukb = "no (UKB-free; other overlap class c)",
    get = function() {
      d <- fread(file.path(RAW, "pad_mvp_phe440_2_EUR.txt.gz"), select = c("SNP_ID", "ci"), showProgress = FALSE)
      d <- d[SNP_ID %in% all_rsid]
      lims <- tstrsplit(d$ci, ",", fixed = TRUE)
      by_rsid(data.table(rsid = d$SNP_ID, seY = (log(as.numeric(lims[[2]])) - log(as.numeric(lims[[1]]))) / (2 * qnorm(0.975))))
    }),
  `CAD MVP Phe_411 EUR (IHD)` = list(cases = 141465, controls = 275809, ukb = "no (UKB-free; other overlap class c)",
    get = function() {
      d <- fread(file.path(RAW, "cad_mvp_phe411_EUR.txt.gz"), select = c("SNP_ID", "ci"), showProgress = FALSE)
      d <- d[SNP_ID %in% all_rsid]
      lims <- tstrsplit(d$ci, ",", fixed = TRUE)
      by_rsid(data.table(rsid = d$SNP_ID, seY = (log(as.numeric(lims[[2]])) - log(as.numeric(lims[[1]]))) / (2 * qnorm(0.975))))
    }),
  `PAD van Zuydam (rounded)` = list(cases = 12086, controls = 449548, ukb = "yes (30.4% / 91.2%)",
    get = function() {
      d <- fread(file.path(RAW, "pad_vanzuydam2021", "final_meta", "pad_primary.txt.gz"), select = c("snp", "or", "or_se"), showProgress = FALSE)
      by_rsid(d[snp %in% all_rsid, .(rsid = snp, seY = or_se / or)])
    }),
  `PAD-adjacent FinnGen I9_ATHSCLE` = list(cases = 19807, controls = 463106, ukb = "no (shares legacy cohorts with EA4/Kweon)",
    get = function() {
      d <- fread(file.path(RAW, "pad_finngen_r12_athscle.tsv.gz"), select = c("rsids", "sebeta"), showProgress = FALSE)
      d <- d[, .(rsid = unlist(strsplit(rsids, ",", fixed = TRUE)), seY = rep(sebeta, lengths(strsplit(rsids, ",", fixed = TRUE))))]
      by_rsid(d)
    }),
  `HF HERMES 2.0` = list(cases = 153174, controls = 1793175, ukb = "yes (8.0% / 22.6%)",
    get = function() {
      d <- fread(file.path(RAW, "hf_hermes2_2025_hfall.tsv.gz"), select = c("rsid", "standard_error"), showProgress = FALSE)
      by_rsid(d[, .(rsid, seY = as.numeric(standard_error))])
    }),
  `Any stroke GIGASTROKE EUR` = list(cases = 73652, controls = 1234808, ukb = "yes (4.0% / 29.9%)",
    get = function() by_pos(file.path(RAW, "stroke_gigastroke_any_eur.tsv.gz"), "chromosome", "base_pair_location", "effect_allele", "other_allele", "standard_error")),
  `Any stroke GIGASTROKE multi` = list(cases = 110182, controls = 1503898, ukb = "yes (2.7% / 24.6%)",
    get = function() by_pos(file.path(RAW, "stroke_gigastroke_any_multi.tsv.gz"), "chromosome", "base_pair_location", "effect_allele", "other_allele", "standard_error"))
)

mvp_get <- function(f) function() {
      d <- fread(file.path(RAW, f), select = c("SNP_ID", "ci"), showProgress = FALSE)
      d <- d[SNP_ID %in% all_rsid]
      lims <- tstrsplit(d$ci, ",", fixed = TRUE)
      by_rsid(data.table(rsid = d$SNP_ID, seY = (log(as.numeric(lims[[2]])) - log(as.numeric(lims[[1]]))) / (2 * qnorm(0.975))))
    }
# GWAS Catalog copies of the MVP gwPheWAS files (GWAS-SSF layout; identical to the dbGaP GIA files: Phe_411 checked
# at 3,556 instrument SNPs and 19,708,268 rows, 2026-10-07). Used because NCBI was unreachable from the sandbox.
mvp_gwascat_get <- function(f) function() {
  d <- fread(file.path(RAW, f), select = c("rsid", "ci_lower", "ci_upper"), showProgress = FALSE)
  d <- d[rsid %in% all_rsid]
  by_rsid(data.table(rsid = d$rsid, seY = (log(as.numeric(d$ci_upper)) - log(as.numeric(d$ci_lower))) / (2 * qnorm(0.975))))
}
if (file.exists(file.path(RAW, "hf_mvp_phe428_EUR.gwascat.tsv.gz")))
  outcomes[["HF MVP Phe_428 EUR"]] <- list(cases = MVP_N$hf[1], controls = MVP_N$hf[2], ukb = "no (UKB-free; class c)", get = mvp_gwascat_get("hf_mvp_phe428_EUR.gwascat.tsv.gz"))
if (file.exists(file.path(RAW, "stroke_mvp_phe433_21_EUR.gwascat.tsv.gz")))
  outcomes[["Stroke MVP Phe_433_21 EUR"]] <- list(cases = MVP_N$st[1], controls = MVP_N$st[2], ukb = "no (UKB-free; class c)", get = mvp_gwascat_get("stroke_mvp_phe433_21_EUR.gwascat.tsv.gz"))
Z <- qnorm(1 - 0.003125 / 2) + qnorm(0.8)
uni <- list(); mv <- list(); oc <- list()
for (o in names(outcomes)) {
  t0 <- Sys.time(); se <- outcomes[[o]]$get()
  logl("{o}: {nrow(se)} instrument rsIDs with outcome SE ({round(as.numeric(difftime(Sys.time(), t0, units = 'secs')))} s)")
  oc[[o]] <- data.table(outcome = o, cases = outcomes[[o]]$cases, controls = outcomes[[o]]$controls,
                        case_fraction = outcomes[[o]]$cases / (outcomes[[o]]$cases + outcomes[[o]]$controls),
                        ukb = outcomes[[o]]$ukb, rsids_with_se = nrow(se[rsid %in% core_rsid]))
  u <- inst[se, on = "rsid", nomatch = 0][, .(n_snps = .N, exp_se = 1 / sqrt(sum(bx^2 / seY^2))), by = exposure]
  u[, `:=`(outcome = o, n_instruments = inst[, .N, by = exposure][u, on = "exposure", N],
           mdor_protective = exp(-Z * exp_se), mdor_harmful = exp(Z * exp_se))]
  uni[[o]] <- u
  for (m in names(Bs)) {
    B <- Bs[[m]][se, on = "rsid", nomatch = 0]
    bcols <- grep("^b_", names(B), value = TRUE)
    if (nrow(B) <= length(bcols)) next
    X <- as.matrix(B[, ..bcols]); W <- 1 / B$seY^2
    V <- solve(crossprod(X * sqrt(W)))
    mv[[paste(o, m)]] <- data.table(outcome = o, model = m, n_snps_model = nrow(Bs[[m]]), n_snps_available = nrow(B),
                                    exposure = sub("^b_", "", bcols), exp_se_conventional = sqrt(diag(V)))
  }
}
uni <- rbindlist(uni); mv <- rbindlist(mv); oc <- rbindlist(oc)
setcolorder(uni, c("outcome", "exposure", "n_instruments", "n_snps"))
write_result_table(oc, "outcome_precision_outcomes")
write_result_table(uni, "outcome_precision_univariable")
write_result_table(mv, "outcome_precision_mvmr")
logl("univariable expected SE:\n{paste(capture.output(print(dcast(uni, exposure ~ outcome, value.var = 'exp_se'), digits = 3)), collapse = '\n')}")
logl("Done.")
close_run_log(log_con)
