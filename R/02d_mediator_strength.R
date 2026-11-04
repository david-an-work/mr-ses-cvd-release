#' 02d_mediator_strength.R
#'
#' Purpose: conditional instrument strength (Sanderson et al. 2021) for the planned
#'   mediation MVMR models (ANALYSIS_PLAN.md §9), before any MR estimate is computed.
#'   Strong univariable instruments do not imply strong conditional instruments, so
#'   this is checked directly. Feasibility only: no outcome data are read.
#'   Models, for each primary SES exposure (EA4, Kweon income, TDI, occupation):
#'     (i)  SES + full mediator panel (BMI, smoking, LDL, HDL, TG, SBP, T2D)  [plan §9]
#'     (ii) SES + each mediator singly                                         [diagnostic]
#'   Instruments: union of the model's exposures' instruments (data/instruments/),
#'   re-clumped jointly on the minimum p (r2 < 0.001, 10 Mb, 1000G EUR), aligned to the
#'   SES exposure's alleles, ambiguous palindromes (MAF > 0.42) dropped. Same procedure
#'   as R/02b_ses_feasibility.R.
#'   Covariance of SNP-effect estimation errors: correlation matrix from the LDSC
#'   cross-trait intercepts estimated jointly for the 4 SES exposures + 7 mediators,
#'   with off-diagonals scaled by 0 (independent samples) / 0.5 / 1 / 1.5.
#'   Caveat: LDL/HDL/TG (GLGC) and T2D (DIAMANTE) are trans-ancestry meta-analyses;
#'   LDSC with European LD scores is approximate for them.
#'
#' Inputs:  data/raw/* (via R/utils/sumstats_registry.R); data/instruments/*;
#'          data/reference/eur_w_ld_chr/; data/reference/1kg_v3/EUR.*;
#'          data/harmonized/large/ldsc/{EA,Income,TDI,Occupation}.sumstats.gz (from 02b)
#' Outputs: results/tables/feasibility_mediator_conditional_F.{tsv,xlsx}
#'          results/tables/feasibility_ldsc_intercepts_mediators.{tsv,xlsx}
#'          results/tables/condF_inputs_mediators/{model}_B.tsv, _pcor.tsv, _condF_unrounded.tsv
#'          data/harmonized/large/ldsc/ (munged mediator sumstats; gitignored)
#'          results/logs/{date}_02d_mediator_strength.log (+ sessionInfo)
#' Depends on: data.table, GenomicSEM, MVMR, ieugwasr, plinkbinr

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
ensure_pkgs_loaded(c("data.table", "GenomicSEM", "MVMR", "ieugwasr", "plinkbinr"))

CLUMP_KB  <- 10000
CLUMP_R2  <- 0.001
PAL_MAF   <- 0.42
BFILE_EUR <- here::here("data", "reference", "1kg_v3", "EUR")
PLINK_BIN <- plinkbinr::get_plink_exe()
LD_DIR    <- paste0(here::here("data", "reference", "eur_w_ld_chr"), "/")
HM3       <- file.path(LD_DIR, "w_hm3.snplist")
OUT_LDSC  <- here::here("data", "harmonized", "large", "ldsc")
DIR_CONDF <- file.path(DIR_TABLES, "condF_inputs_mediators"); fs::dir_create(DIR_CONDF)
PCOR_SCALES <- c(0, 0.5, 1, 1.5)

source(here::here("R", "utils", "sumstats_registry.R"))

log_con <- log_run("02d_mediator_strength")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)

ses <- c("ea_okbay2022", "income_kweon2025", "tdi_panukb", "occupation_akimova2024")
med <- c("bmi_giant_yengo2018", "smoking_wootton2020", "ldl_glgc_graham2021", "hdl_glgc_graham2021",
         "tg_glgc_graham2021", "sbp_keaton2024", "t2d_diamante_mahajan2022")
short <- c(ea_okbay2022 = "EA", income_kweon2025 = "Income", tdi_panukb = "TDI", occupation_akimova2024 = "Occupation",
           bmi_giant_yengo2018 = "BMI", smoking_wootton2020 = "Smoking", ldl_glgc_graham2021 = "LDL",
           hdl_glgc_graham2021 = "HDL", tg_glgc_graham2021 = "TG", sbp_keaton2024 = "SBP",
           t2d_diamante_mahajan2022 = "T2D")
traits <- c(ses, med)

bim <- fread(paste0(BFILE_EUR, ".bim"), header = FALSE,
             col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"),
             colClasses = c("character", "character", "numeric", "integer", "character", "character"))
bim <- bim[grepl("^rs[0-9]+$", rsid)]
hm3 <- fread(HM3)

# ---- LDSC: munge mediators (fingerprinted cache, as in 02b), joint 11-trait ldsc ----

munged_med <- file.path(OUT_LDSC, paste0(short[med], ".sumstats.gz"))
munged_ses <- file.path(OUT_LDSC, paste0(short[ses], ".sumstats.gz"))
stopifnot("run R/02b_ses_feasibility.R first (SES munged files missing)" = all(file.exists(munged_ses)))
bim_fi <- file.info(paste0(BFILE_EUR, ".bim"))
fingerprint <- c(
  paste0("reader_md5=", unname(tools::md5sum(here::here("R", "utils", "sumstats_registry.R")))),
  "munge=info.filter:0.9|maf.filter:0.01",
  paste0("hm3_md5=", unname(tools::md5sum(HM3))),
  paste0("bim=", bim_fi$size, "|", format(bim_fi$mtime, "%Y-%m-%dT%H:%M:%S")),
  paste0("GenomicSEM=", packageVersion("GenomicSEM"), "|data.table=", packageVersion("data.table")),
  paste0("input=", vapply(med, function(k) { fi <- file.info(registry[[k]]$file)
    paste(k, basename(registry[[k]]$file), fi$size, format(fi$mtime, "%Y-%m-%dT%H:%M:%S"), sep = "|") }, character(1))))
fp_file <- file.path(OUT_LDSC, "munge_fingerprint_mediators.txt")
cache_ok <- all(file.exists(munged_med)) && file.exists(fp_file) && identical(readLines(fp_file), fingerprint)
logl("mediator munge cache: {if (cache_ok) 'valid, reusing' else 'missing or stale, rebuilding'}")
if (!cache_ok) {
  file.remove(c(munged_med[file.exists(munged_med)], fp_file[file.exists(fp_file)]))
  inputs <- character(0)
  for (k in med) {
    d <- read_registry(k, p_max = 1, rsids = hm3$SNP, bim = bim[rsid %in% hm3$SNP])
    out <- d[, .(SNP = rsid, A1 = ea, A2 = oa, BETA = beta, P = pval, N = n,
                 MAF = ifelse(is.na(eaf), NA_real_, pmin(eaf, 1 - eaf)))]
    if (all(is.na(out$MAF))) out[, MAF := NULL]
    f <- file.path(OUT_LDSC, paste0(short[k], "_input.tsv"))
    fwrite(out, f, sep = "\t", na = "NA"); inputs <- c(inputs, f)
    logl("ldsc input {k}: {nrow(out)} HapMap3 SNPs; N missing for {sum(is.na(out$N))}")
  }
  owd <- setwd(OUT_LDSC)
  GenomicSEM::munge(files = inputs, hm3 = HM3, trait.names = unname(short[med]),
                    N = rep(NA, length(med)), info.filter = 0.9, maf.filter = 0.01, log.name = "munge_mediators")
  setwd(owd); file.remove(inputs)
  writeLines(fingerprint, fp_file)
}

owd <- setwd(OUT_LDSC)
ld_res <- GenomicSEM::ldsc(traits = c(munged_ses, munged_med), sample.prev = rep(NA, length(traits)),
                           population.prev = rep(NA, length(traits)), ld = LD_DIR, wld = LD_DIR,
                           trait.names = unname(short[traits]), stand = TRUE, ldsc.log = "ldsc_ses_mediators")
setwd(owd)
I <- ld_res$I; dimnames(I) <- list(unname(short[traits]), unname(short[traits]))
write_result_table(cbind(data.table(trait = rownames(I)), as.data.table(round(I, 4))), "feasibility_ldsc_intercepts_mediators")
logl("LDSC intercepts (SES + mediators):\n{paste(capture.output(print(round(I, 3))), collapse = '\n')}")

# ---- Conditional F ------------------------------------------------------------------

inst <- rbindlist(lapply(traits, function(key) {
  x <- fread(file.path(DIR_INSTRUMENTS, glue("{key}_instruments.tsv")), colClasses = list(character = "chr"))
  x[, .(rsid, pval, dataset = key)]
}))
is_pal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C")

# read each trait once, restricted to every SNP that can enter any model
cand <- unique(inst$rsid)
full <- lapply(setNames(traits, traits), function(k) {
  logl("reading {k} at {length(cand)} candidate SNPs")
  read_registry(k, p_max = 1, rsids = cand, bim = bim[rsid %in% cand])[, .(rsid, ea, oa, eaf, beta, se)]
})

run_model <- function(keys) {
  u <- inst[dataset %in% keys, .(pval = min(pval)), by = rsid]
  cl <- ieugwasr::ld_clump(data.frame(rsid = u$rsid, pval = u$pval, id = "joint"),
                           clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, clump_p = 1,
                           bfile = BFILE_EUR, plink_bin = PLINK_BIN)
  snps <- cl$rsid
  ref <- full[[keys[1]]][rsid %in% snps, .(rsid, ref_ea = ea, ref_oa = oa, ref_eaf = eaf)]
  B <- copy(ref)
  for (key in keys) {
    m <- merge(ref, full[[key]][rsid %in% snps, .(rsid, ea, oa, beta, se)], by = "rsid")
    m[, aligned := fifelse(ea == ref_ea & oa == ref_oa, beta,
                    fifelse(ea == ref_oa & oa == ref_ea, -beta, NA_real_))]
    B <- merge(B, m[!is.na(aligned), .(rsid, b = aligned, s = se)], by = "rsid")
    setnames(B, c("b", "s"), c(paste0("b_", key), paste0("se_", key)))
  }
  maf <- pmin(B$ref_eaf, 1 - B$ref_eaf)
  B <- B[!(is_pal(ref_ea, ref_oa) & (is.na(maf) | maf > PAL_MAF))]
  bx <- as.matrix(B[, paste0("b_", keys), with = FALSE])
  sx <- as.matrix(B[, paste0("se_", keys), with = FALSE])
  model_id <- paste(short[keys], collapse = "+")
  pcor <- cov2cor(I[short[keys], short[keys], drop = FALSE])
  fwrite(B, file.path(DIR_CONDF, paste0(model_id, "_B.tsv")), sep = "\t")
  fwrite(cbind(data.table(trait = short[keys]), as.data.table(pcor)),
         file.path(DIR_CONDF, paste0(model_id, "_pcor.tsv")), sep = "\t")
  r_input <- MVMR::format_mvmr(BXGs = bx, BYG = rep(0, nrow(B)), seBXGs = sx,
                               seBYG = rep(1, nrow(B)), RSID = B$rsid)
  res <- data.table(model_type = if (length(keys) > 2) "full panel" else "single mediator",
                    model = paste(short[keys], collapse = " + "), ses_exposure = short[keys[1]],
                    exposure = short[keys], role = c("SES", rep("mediator", length(keys) - 1)),
                    n_snps_clumped = length(snps), n_snps = nrow(B))
  for (sc in PCOR_SCALES) {
    pc <- pcor * sc; diag(pc) <- 1
    col <- paste0("condF_pcor_x", sc)
    if (min(eigen(pc, symmetric = TRUE, only.values = TRUE)$values) <= 0) { res[, (col) := NA_real_]; next }
    f <- suppressWarnings(if (sc == 0) MVMR::strength_mvmr(r_input, gencov = 0)
                          else MVMR::strength_mvmr(r_input, gencov = MVMR::phenocov_mvmr(pc, sx)))
    res[, (col) := as.numeric(f[1, ])]
  }
  fwrite(res, file.path(DIR_CONDF, paste0(model_id, "_condF_unrounded.tsv")), sep = "\t")
  res
}

models <- c(lapply(ses, function(s) c(s, med)),
            unlist(lapply(ses, function(s) lapply(med, function(m) c(s, m))), recursive = FALSE))
cf <- rbindlist(lapply(models, function(m) { logl("conditional F: {paste(short[m], collapse = ' + ')}"); run_model(m) }))
fc <- grep("^condF_", names(cf), value = TRUE)
cf[, (fc) := lapply(.SD, round, 2), .SDcols = fc]
cf[, min_condF_any_scale := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = fc]
write_result_table(cf, "feasibility_mediator_conditional_F")
logl("conditional F (mediator models):\n{paste(capture.output(print(cf[, c('model', 'exposure', 'n_snps', fc), with = FALSE])), collapse = '\n')}")
logl("Done.")
close_run_log(log_con)
