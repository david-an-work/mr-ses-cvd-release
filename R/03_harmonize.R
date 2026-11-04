#' 03_harmonize.R
#'
#' Purpose: Phase 3 harmonization under the frozen plan v1.0 (§6; HARMONIZATION_RULES_draft.md draft v3, rules 1-11).
#'   Registered on OSF (k7dcm; plan files in usz4x) before this script was run. No MR estimate is computed here.
#' Steps:
#'   1. Outcome extraction: for every outcome file, rows for the union of instrument SNPs (all exposure, mediator and
#'      MVMR SNP sets), matched by rule 1 (GRCh37 files: chr + GRCh37 pos + unordered allele pair, rsID cross-check;
#'      GRCh38 files [MVP, FinnGen]: rsID + unordered allele pair). Rule 3: exactly one allele-matching row, else drop.
#'      MVP SE from the 95% CI (rule 7); van Zuydam beta = ln(OR), SE = or_se/OR, sensitivity SE kept (rule 11).
#'   2. Harmonization (rules 2, 4): align outcome to the exposure effect allele; no strand-complement rescue;
#'      palindromes by TwoSampleMR action = 2 logic (tolerance 0.08; van Zuydam 0.085; MVP AFR: all palindromes
#'      dropped). Implemented here and cross-checked against TwoSampleMR::harmonise_data on every standard cell.
#'   3. Cells: univariable exposure x outcome; MVMR SNP sets (R/02b, R/02d, R/02e) x outcome, aligned to each set's
#'      reference allele; EA4 instruments x each mediator GWAS (beta(EA -> M) for the exploratory mediation).
#'   4. Flags carried, not applied (sensitivity analyses): MHC (chr6:25-34 Mb, GRCh37), MVP imputation r2 < 0.8.
#'   5. Retention log per cell (rule 10) and per-SNP drop log.
#' Outputs:
#'   data/harmonized/large/outcome_extracts/{outcome}.tsv.gz
#'   data/harmonized/large/{univariable,mvmr,mediation}/{set}__{outcome}.tsv.gz
#'   results/tables/harmonization_retention.{tsv,xlsx}; results/tables/harmonization_drop_log.tsv
#'   results/tables/harmonization_tsmr_check.tsv; results/logs/{date}_03_harmonize.log
#' Depends on: data.table, TwoSampleMR (cross-check only)

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
log_con <- log_run("03_harmonize")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
source(here::here("R", "utils", "sumstats_registry.R"))
RAW <- DIR_DATA_RAW
OUT <- file.path(DIR_DATA_HARM, "large")
for (d in c("outcome_extracts", "univariable", "mvmr", "mediation")) fs::dir_create(file.path(OUT, d))
TOL <- 0.08; TOL_VZ <- 0.085
MHC <- c(chr = "6", from = 25e6, to = 34e6)

# ---- 1. instrument sets ---------------------------------------------------------------------------------------
uni_keys <- c(EA4 = "ea_okbay2022", Kweon = "income_kweon2025", TDI = "tdi_panukb", CAMSIS = "occupation_akimova2024",
              EA3 = "ea_lee2018", Hill = "income_hill2019", TDI_EUR = "tdi_panukb_eur", ISEI = "occupation_akimova2024_isei",
              gSES = "ses_gses_xia2025")
med_keys <- c(BMI = "bmi_giant_yengo2018", Smoking = "smoking_wootton2020", LDL = "ldl_glgc_graham2021",
              HDL = "hdl_glgc_graham2021", TG = "tg_glgc_graham2021", SBP = "sbp_keaton2024", T2D = "t2d_diamante_mahajan2022")
read_inst <- function(key) {
  d <- fread(file.path(DIR_INSTRUMENTS, paste0(key, "_instruments.tsv")))
  d[, .(rsid, chr = as.character(chr), pos = as.integer(pos), ea_x = toupper(effect_allele), oa_x = toupper(other_allele),
        eaf_x = as.numeric(eaf), beta_x = beta, se_x = se, p_x = pval)]
}
inst <- lapply(c(uni_keys, med_keys), read_inst)
names(inst) <- c(names(uni_keys), names(med_keys))
bfiles <- c(list.files(file.path(DIR_TABLES, "condF_inputs"), "_B\\.tsv$", full.names = TRUE),
            list.files(file.path(DIR_TABLES, "condF_inputs_mediators"), "_B\\.tsv$", full.names = TRUE),
            file.path(DIR_TABLES, "pilot_lid_cad_inputs", c("S1_B.tsv", "S1-LDpruned_B.tsv")))
stopifnot(all(file.exists(bfiles)))
Bs <- lapply(bfiles, fread); names(Bs) <- sub("_B\\.tsv$", "", basename(bfiles))
names(Bs)[names(Bs) %in% c("S1", "S1-LDpruned")] <- paste0("LiD_", names(Bs)[names(Bs) %in% c("S1", "S1-LDpruned")])
snp_ref <- unique(rbindlist(inst)[, .(rsid, chr, pos, a1 = ea_x, a2 = oa_x)], by = "rsid")
b_rs <- unique(unlist(lapply(Bs, `[[`, "rsid")))
stopifnot(all(b_rs %in% snp_ref$rsid))
logl("union SNPs: {nrow(snp_ref)} (instrument sets {length(inst)}, MVMR sets {length(Bs)})")
source(here::here("R", "utils", "outcome_readers.R"))
rd <- make_outcome_readers(snp_ref, RAW, logl)
outcomes <- rd$outcomes; by_rsid <- rd$by_rsid; allele_pair_ok <- rd$allele_pair_ok

ext <- list()
for (o in names(outcomes)) {
  f <- file.path(OUT, "outcome_extracts", paste0(o, ".tsv.gz"))
  t0 <- Sys.time(); e <- outcomes[[o]]()
  fwrite(e, f, sep = "\t")
  ext[[o]] <- e
  logl("{o}: {sum(e$status == 'ok')} ok / {nrow(e)} union SNPs; missing {sum(e$status == 'missing')}, allele_mismatch {sum(e$status == 'allele_mismatch')}, duplicate {sum(e$status == 'duplicate')} ({round(as.numeric(difftime(Sys.time(), t0, units = 'secs')))} s)")
}
# mediators as 'outcomes' for beta(EA4 -> M): rsID + alleles (registry readers, sign applied)
med_ext <- list()
for (m in names(med_keys)) {
  d <- read_registry(med_keys[[m]], p_max = 1, rsids = inst$EA4$rsid)
  med_ext[[m]] <- by_rsid(d[, .(rsid, ea_y = ea, oa_y = oa, eaf_y = eaf, beta_y = beta, se_y = se, p_y = pval)])
  logl("mediator {m}: {sum(med_ext[[m]]$status == 'ok')} of {nrow(inst$EA4)} EA4 instruments ok")
}

# ---- 3. harmonization ------------------------------------------------------------------------------------------
source(here::here("R", "utils", "harmonize_utils.R"))
ret <- list(); drops <- list(); chk <- list()
log_drops <- function(h, set, outcome) if (any(!h$mr_keep)) h[!(mr_keep), .(set, outcome, rsid, drop_reason)] else NULL
tsmr_check <- function(x, y, set, outcome) {
  # cross-check against TwoSampleMR action = 2 on the SNPs with an allele-matched outcome row
  yy <- y[status == "ok"]; xx <- x[rsid %in% yy$rsid]
  if (nrow(xx) == 0) return(NULL)
  ed <- data.frame(SNP = xx$rsid, beta.exposure = xx$beta_x, se.exposure = xx$se_x, effect_allele.exposure = xx$ea_x,
                   other_allele.exposure = xx$oa_x, eaf.exposure = xx$eaf_x, exposure = "x", id.exposure = "x", stringsAsFactors = FALSE)
  od <- data.frame(SNP = yy$rsid, beta.outcome = yy$beta_y, se.outcome = yy$se_y, effect_allele.outcome = yy$ea_y,
                   other_allele.outcome = yy$oa_y, eaf.outcome = yy$eaf_y, outcome = "y", id.outcome = "y", stringsAsFactors = FALSE)
  t <- suppressMessages(TwoSampleMR::harmonise_data(ed, od, action = 2))
  data.table(set = set, outcome = outcome, tsmr_keep = sum(t$mr_keep), tsmr_beta_y = list(setNames(t$beta.outcome[t$mr_keep], t$SNP[t$mr_keep])))
}
std_outcomes <- setdiff(names(ext), c("pad_vanzuydam", "pad_mvp_afr"))
tol_for <- function(o) if (o == "pad_vanzuydam") TOL_VZ else TOL

# 3a. univariable cells
for (s in names(uni_keys)) for (o in names(ext)) {
  h <- harmonize(inst[[s]], ext[[o]], tol = tol_for(o), drop_all_pal = o == "pad_mvp_afr")
  fwrite(h, file.path(OUT, "univariable", paste0(s, "__", o, ".tsv.gz")), sep = "\t")
  ret[[length(ret) + 1]] <- retention(h, s, o, "univariable"); drops[[length(drops) + 1]] <- log_drops(h, s, o)
  if (o %in% std_outcomes) {
    t <- tsmr_check(inst[[s]], ext[[o]], s, o)
    mine <- setNames(h$beta_y_aligned[h$mr_keep], h$rsid[h$mr_keep])
    tb <- t$tsmr_beta_y[[1]]
    chk[[length(chk) + 1]] <- data.table(set = s, outcome = o, n_mine = length(mine), n_tsmr = length(tb),
      same_snps = setequal(names(mine), names(tb)),
      max_abs_beta_diff = if (length(intersect(names(mine), names(tb)))) max(abs(mine[intersect(names(mine), names(tb))] - tb[intersect(names(mine), names(tb))])) else NA_real_)
  }
}
logl("univariable cells: {length(uni_keys) * length(ext)}")

# 3b. MVMR SNP sets: align outcome to each set's reference allele
for (b in names(Bs)) for (o in names(ext)) {
  B <- Bs[[b]]
  x <- B[, c(list(rsid = rsid, ea_x = toupper(ref_ea), oa_x = toupper(ref_oa), eaf_x = as.numeric(ref_eaf)), .SD), .SDcols = patterns("^(b|se)_")]
  x <- merge(x, snp_ref[, .(rsid, chr, pos)], by = "rsid")
  h <- harmonize(x, ext[[o]], tol = tol_for(o), drop_all_pal = o == "pad_mvp_afr")
  fwrite(h, file.path(OUT, "mvmr", paste0(b, "__", o, ".tsv.gz")), sep = "\t")
  ret[[length(ret) + 1]] <- retention(h, b, o, "mvmr"); drops[[length(drops) + 1]] <- log_drops(h, b, o)
}
logl("MVMR cells: {length(Bs) * length(ext)}")

# 3c. EA4 -> mediator (exposure-side step of the product of coefficients)
for (m in names(med_ext)) {
  h <- harmonize(inst$EA4, med_ext[[m]])
  fwrite(h, file.path(OUT, "mediation", paste0("EA4__", m, ".tsv.gz")), sep = "\t")
  ret[[length(ret) + 1]] <- retention(h, "EA4", m, "EA4_to_mediator"); drops[[length(drops) + 1]] <- log_drops(h, "EA4", m)
}

# ---- 4. outputs -----------------------------------------------------------------------------------------------
R <- rbindlist(ret); C <- rbindlist(chk); Dl <- rbindlist(drops)
write_result_table(R, "harmonization_retention")
fwrite(Dl, file.path(DIR_TABLES, "harmonization_drop_log.tsv"), sep = "\t")
fwrite(C, file.path(DIR_TABLES, "harmonization_tsmr_check.tsv"), sep = "\t")
logl("TwoSampleMR cross-check: {sum(C$same_snps)}/{nrow(C)} cells with identical retained SNPs; max |beta diff| {signif(max(C$max_abs_beta_diff, na.rm = TRUE), 3)}")
bad <- C[!same_snps | is.na(max_abs_beta_diff) | max_abs_beta_diff > 1e-12]
if (nrow(bad) > 0) { print(bad); stop("TwoSampleMR cross-check failed in ", nrow(bad), " cell(s); outputs above must not be used") }
u <- R[type == "univariable"]
logl("univariable retention (final/selected), by outcome:\n{paste(capture.output(print(dcast(u, set ~ outcome, value.var = 'n_final'))), collapse = '\n')}")
logl("Done. No MR estimate computed.")
close_run_log(log_con)
