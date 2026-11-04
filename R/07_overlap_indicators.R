#' 07_overlap_indicators.R
#'
#' Purpose: the registered sample-overlap indicators (plan v1.0 §2.1, Steps 2-3). Descriptive only; nothing here
#'   corrects any estimate.
#'   Step 3 - illustrative Burgess-type scenario per univariable cell: phi_UKB x (1/F) x 100 (% of the observational
#'     association), with F recomputed on the outcome-specific harmonized set as ((N - K - 1)/K) x R2/(1 - R2), using the
#'     instruments' per-SNP r2 and median N (same formula as R/02 INSTRUMENT_SUMMARY). phi_UKB = min(UKB N in exposure
#'     GWAS, UKB N in outcome GWAS) / total outcome N, from docs/prereg/d7_evidence.tsv. MVP and FinnGen outcomes:
#'     phi_UKB = 0 (UKB-mediated only; total overlap unknown). Not calculated for Kweon (UKB N of the income factor not
#'     documented) or gSES (§7.4).
#'   Step 2 - cross-trait LDSC intercepts, exposure x outcome (GenomicSEM::ldsc, HapMap3, 1000G EUR LD scores), using the
#'     EUR outcome file where one exists. HERMES 2.0 (multi-ancestry; EUR-only not obtained) and van Zuydam (rounded)
#'     are labelled approximate. Outcome N for munging = effective N 4/(1/cases + 1/controls) from MANIFEST; the
#'     intercept does not depend on this choice. Observed scale (no prevalence conversion).
#' Outputs: results/tables/overlap_scenarios.tsv ; results/tables/overlap_ldsc_exposure_outcome.tsv
#' Depends on: data.table, GenomicSEM

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
log_con <- log_run("07_overlap_indicators")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
source(here::here("R", "utils", "outcome_readers.R"))
H <- file.path(DIR_DATA_HARM, "large"); RAW <- DIR_DATA_RAW
MODE <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "all" }

# ---- Step 3: illustrative scenarios -------------------------------------------------------------------------------
ukb_exp <- c(EA4 = 441121, TDI = 439198, CAMSIS = 273157, Hill = 286301, EA3 = 442183)   # d7_evidence item 1.UKB-bound
out_ukb <- data.table(outcome = c("cad_aragam_multi", "cad_aragam_eur", "hf_hermes2", "stroke_gigastroke_multi",
                                  "stroke_gigastroke_eur", "pad_vanzuydam"),
  ukb_n = c(472335, 471900, 421803, 372351, 372351, 413616),
  total_n = c(210842 + 1167328, 181522 + 984168, 153174 + 1793175, 110182 + 1503898, 73652 + 1234808, 12086 + 449548),
  source = c("Aragam ST1", "Aragam EUR shares 18.7%/44.5% (plan §2.1)", "HERMES 2.0 MOESM6", "GIGASTROKE ST1",
             "GIGASTROKE ST1 (UKB only in EUR)", "van Zuydam Table I"))
zero_out <- c("cad_mvp_phe411_eur", "pad_mvp_eur", "pad_mvp_afr", "hf_mvp_phe428_eur", "stroke_mvp_phe433_21_eur", "pad_finngen_athscle")
keys <- c(EA4 = "ea_okbay2022", Kweon = "income_kweon2025", TDI = "tdi_panukb", CAMSIS = "occupation_akimova2024",
          EA3 = "ea_lee2018", Hill = "income_hill2019", TDI_EUR = "tdi_panukb_eur", gSES = "ses_gses_xia2025")
if (MODE %in% c("all", "scenarios")) {
  rows <- list()
  for (s in names(keys)) {
    inst <- fread(file.path(DIR_INSTRUMENTS, paste0(keys[[s]], "_instruments.tsv")), select = c("rsid", "r2", "n"))
    Nm <- median(inst$n, na.rm = TRUE)
    for (f in list.files(file.path(H, "univariable"), paste0("^", s, "__.*\\.tsv\\.gz$"), full.names = TRUE)) {
      o <- sub("^.*__(.*)\\.tsv\\.gz$", "\\1", f)
      d <- fread(f, select = c("rsid", "mr_keep"))[mr_keep == TRUE]
      r <- inst[rsid %in% d$rsid]; K <- nrow(r); R2 <- sum(r$r2)
      Fk <- ((Nm - K - 1) / K) * (R2 / (1 - R2))
      phi <- if (o %in% zero_out) 0 else if (s %in% names(ukb_exp) && o %in% out_ukb$outcome) {
        u <- out_ukb[outcome == o]; min(ukb_exp[[s]], u$ukb_n) / u$total_n } else NA_real_
      why <- if (o %in% zero_out) "UKB-free outcome (UKB-mediated only; total overlap unknown)" else
             if (s == "Kweon") "not calculated (UKB N of the income factor not documented)" else
             if (s == "gSES") "not calculated (§7.4)" else if (s == "TDI_EUR") "not calculated (UKB N of EUR stratum not in d7_evidence)" else ""
      if (s == "TDI_EUR") phi <- if (o %in% zero_out) 0 else NA_real_
      rows[[length(rows) + 1]] <- data.table(set = s, outcome = o, n_snps = K, F_outcome_specific = Fk,
        phi_ukb = phi, scenario_pct = 100 * phi / Fk, note = why)
    }
  }
  S <- rbindlist(rows); fwrite(S, file.path(DIR_TABLES, "overlap_scenarios.tsv"), sep = "\t")
  logl("scenarios: {nrow(S)} cells; calculated {sum(!is.na(S$scenario_pct))}")
}

# ---- Step 2: exposure x outcome cross-trait LDSC intercepts -------------------------------------------------------
if (MODE %in% c("all", "ldsc")) {
  LD_DIR <- paste0(here::here("data", "reference", "eur_w_ld_chr"), "/"); HM3 <- file.path(LD_DIR, "w_hm3.snplist")
  OUT_LDSC <- file.path(H, "ldsc", "outcomes"); fs::dir_create(OUT_LDSC)
  hm3 <- fread(HM3)
  bim <- fread(here::here("data", "reference", "1kg_v3", "EUR.bim"), header = FALSE, select = c(1, 2, 4),
               col.names = c("chr", "rsid", "pos"), colClasses = list(character = c(1, 2)))[rsid %in% hm3$SNP]
  snp_ref <- unique(merge(bim, hm3[, .(rsid = SNP, a1 = A1, a2 = A2)], by = "rsid"), by = "rsid")
  rd <- make_outcome_readers(snp_ref, RAW, logl)
  man <- fread(PATH_MANIFEST)
  ld_outcomes <- c(cad_aragam_eur = "cad_aragam2022_eur", pad_mvp_eur = "pad_mvp_phe440_2_eur", hf_hermes2 = "hf_hermes2_2025",
                   stroke_gigastroke_eur = "stroke_gigastroke_any_eur", cad_mvp_phe411_eur = "cad_mvp_phe411_eur",
                   hf_mvp_phe428_eur = "hf_mvp_phe428_eur", stroke_mvp_phe433_21_eur = "stroke_mvp_phe433_21_eur",
                   pad_finngen_athscle = "pad_finngen_r12_athscle", pad_vanzuydam = "pad_vanzuydam2021")
  files <- character()
  for (o in names(ld_outcomes)) {
    f <- file.path(OUT_LDSC, paste0(o, "_input.tsv"))
    if (!file.exists(paste0(file.path(OUT_LDSC, o), ".sumstats.gz"))) {
      d <- rd$outcomes[[o]]()[status == "ok"]
      mm <- man[dataset == ld_outcomes[[o]]]; neff <- 4 / (1 / as.numeric(mm$cases) + 1 / as.numeric(mm$controls))
      fwrite(d[, .(SNP = rsid, A1 = ea_y, A2 = oa_y, BETA = beta_y, SE = se_y, P = 2 * pnorm(-abs(beta_y / se_y)), N = neff)], f, sep = "\t")
      owd <- setwd(OUT_LDSC)
      GenomicSEM::munge(files = f, hm3 = HM3, trait.names = o, info.filter = 0.9, maf.filter = 0.01, log.name = paste0("munge_", o))
      setwd(owd); file.remove(f)
      logl("munged {o}: {nrow(d)} HapMap3 rows; N_eff {round(neff)}")
    }
    files <- c(files, paste0(file.path(OUT_LDSC, o), ".sumstats.gz"))
  }
  exp_files <- file.path(H, "ldsc", c("EA", "Income", "TDI", "Occupation"), fsep = "/")
  exp_files <- paste0(exp_files, ".sumstats.gz"); stopifnot(all(file.exists(exp_files)))
  nm <- c("EA4", "Kweon", "TDI", "CAMSIS", names(ld_outcomes))
  owd <- setwd(OUT_LDSC)
  L <- GenomicSEM::ldsc(traits = c(exp_files, files), sample.prev = rep(NA, length(nm)), population.prev = rep(NA, length(nm)),
                        ld = LD_DIR, wld = LD_DIR, trait.names = nm, stand = TRUE, ldsc.log = "ldsc_exposure_outcome")
  setwd(owd)
  I <- L$I; dimnames(I) <- list(nm, nm); rg <- L$S_Stand; dimnames(rg) <- list(nm, nm)
  res <- CJ(exposure = nm[1:4], outcome = names(ld_outcomes), sorted = FALSE)
  res[, `:=`(cross_trait_intercept = mapply(function(e, o) I[e, o], exposure, outcome),
             rg = mapply(function(e, o) rg[e, o], exposure, outcome),
             outcome_intercept_h2 = mapply(function(o) I[o, o], outcome),
             label = fifelse(outcome %in% c("hf_hermes2", "pad_vanzuydam"), "approximate", ""))]
  fwrite(res, file.path(DIR_TABLES, "overlap_ldsc_exposure_outcome.tsv"), sep = "\t")
  logl("LDSC exposure x outcome: {nrow(res)} pairs")
}

# Intercept standard errors are not returned by GenomicSEM::ldsc(); they are parsed from its log
# ("Cross trait Intercept: est (se)"), keyed by the two sumstats file names.
if (MODE %in% c("all", "ldsc", "ldsc_se")) {
  OUT_LDSC <- file.path(H, "ldsc", "outcomes")
  lg <- readLines(file.path(OUT_LDSC, "ldsc_exposure_outcome_ldsc.log"))
  key <- c(EA = "EA4", Income = "Kweon", TDI = "TDI", Occupation = "CAMSIS")
  hdr <- grep("^Results for genetic covariance between: ", lg)
  pr <- rbindlist(lapply(hdr, function(i) {
    f <- strsplit(sub("^Results for genetic covariance between: ", "", lg[i]), " and ", fixed = TRUE)[[1]]
    b <- sub("\\.sumstats\\.gz$", "", basename(f))
    ic <- grep("^Cross trait Intercept: ", lg[i:(i + 4)], value = TRUE)[1]
    v <- as.numeric(regmatches(ic, gregexpr("-?[0-9.]+(e-?[0-9]+)?", ic))[[1]])
    data.table(t1 = b[1], t2 = b[2], est = v[1], se = v[2])
  }))
  pr[, exposure := key[t1]]
  pr <- pr[!is.na(exposure) & !(t2 %in% names(key)), .(exposure, outcome = t2, est_log = est, cross_trait_intercept_se = se)]
  res <- fread(file.path(DIR_TABLES, "overlap_ldsc_exposure_outcome.tsv"))
  if ("cross_trait_intercept_se" %in% names(res)) res[, cross_trait_intercept_se := NULL]
  res <- merge(res, pr, by = c("exposure", "outcome"), all.x = TRUE, sort = FALSE)
  stopifnot(!anyNA(res$cross_trait_intercept_se), max(abs(res$est_log - res$cross_trait_intercept)) < 6e-4)
  res[, est_log := NULL]
  setcolorder(res, c("exposure", "outcome", "cross_trait_intercept", "cross_trait_intercept_se"))
  fwrite(res, file.path(DIR_TABLES, "overlap_ldsc_exposure_outcome.tsv"), sep = "\t")
  logl("intercept SEs parsed for {nrow(res)} pairs")
}
close_run_log(log_con)
