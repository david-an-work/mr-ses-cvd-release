#' 02b_ses_feasibility.R
#'
#' Purpose: feasibility checks for the proposed (not adopted) reframing in
#'   literature/SCOPE_OPTIONS.md (TASKS.md T2.11-T2.12). These measure; they do
#'   not decide the framing (decision D5).
#'   T2.11  Pairwise genetic correlations (LD score regression via GenomicSEM)
#'          among the SES exposures, on HapMap3 SNPs with 1000G EUR LD scores.
#'   T2.12  Conditional F-statistics (Sanderson et al. 2021, MVMR::strength_mvmr)
#'          for the SES exposures modelled jointly: the 4-exposure model and EA
#'          paired with each other exposure. Two versions: independent-sample
#'          assumption (gencov = 0) and overlap-adjusted (covariance of SNP-
#'          exposure estimates from the LDSC cross-trait intercepts).
#'          Covariance sensitivity (added 2026-10-03): the off-diagonal phenotypic
#'          correlations are scaled by 0, 0.5, 1 and 1.5 (gencov = 0 is scale 0).
#'          EUR-only TDI sensitivity (Pan-UKB beta_EUR/se_EUR) runs alongside the meta.
#'          The exact SNP set, aligned B/SE matrix and correlation matrix behind
#'          every model are written to results/tables/condF_inputs/ so the F arithmetic
#'          can be re-checked without MVMR (R/02c_verify_conditional_F.R; upstream steps
#'          such as SNP selection and allele alignment are shared, not re-checked).
#'
#' Inputs:  data/raw/* (via R/utils/sumstats_registry.R); data/instruments/*;
#'          data/reference/eur_w_ld_chr/; data/reference/1kg_v3/EUR.*
#' Outputs: results/tables/feasibility_rg_ldsc.{tsv,xlsx}
#'          results/tables/feasibility_ldsc_intercepts.{tsv,xlsx}
#'          results/tables/feasibility_conditional_F.{tsv,xlsx}
#'          results/tables/condF_inputs/{model}_B.tsv, {model}_pcor.tsv, {model}_condF_unrounded.tsv
#'          data/harmonized/large/ldsc/ (munged sumstats; gitignored)
#'          results/logs/{date}_02b_ses_feasibility.log (+ sessionInfo)
#' Depends on: data.table, GenomicSEM, MVMR, ieugwasr, plinkbinr

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
ensure_pkgs_loaded(c("data.table", "GenomicSEM", "MVMR", "ieugwasr", "plinkbinr"))

P_THRESH  <- 5e-8
CLUMP_KB  <- 10000
CLUMP_R2  <- 0.001
PAL_MAF   <- 0.42
BFILE_EUR <- here::here("data", "reference", "1kg_v3", "EUR")
PLINK_BIN <- plinkbinr::get_plink_exe()
LD_DIR    <- paste0(here::here("data", "reference", "eur_w_ld_chr"), "/")
HM3       <- file.path(LD_DIR, "w_hm3.snplist")
OUT_LDSC  <- here::here("data", "harmonized", "large", "ldsc")
fs::dir_create(OUT_LDSC)

source(here::here("R", "utils", "sumstats_registry.R"))

log_con <- log_run("02b_ses_feasibility")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)  # sessionInfo() is written by close_run_log() at the end

ses_primary <- c("ea_okbay2022", "income_kweon2025", "tdi_panukb", "occupation_akimova2024")
ses_all     <- c(ses_primary, "income_hill2019", "ea_lee2018", "tdi_panukb_eur", "occupation_akimova2024_isei")
short <- c(ea_okbay2022 = "EA", income_kweon2025 = "Income", tdi_panukb = "TDI",
           occupation_akimova2024 = "Occupation", income_hill2019 = "Income_Hill", ea_lee2018 = "EA_Lee",
           tdi_panukb_eur = "TDI_EUR", occupation_akimova2024_isei = "Occupation_ISEI")

bim <- fread(paste0(BFILE_EUR, ".bim"), header = FALSE,
             col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"),
             colClasses = c("character", "character", "numeric", "integer", "character", "character"))
bim <- bim[grepl("^rs[0-9]+$", rsid)]
hm3 <- fread(HM3)

# ---- T2.11: LD score regression --------------------------------------------------

munged <- file.path(OUT_LDSC, paste0(short[ses_all], ".sumstats.gz"))
# Cache key (added 2026-10-03; previously existence-only): the munged files are reused only
# if the reader code, the munge settings, the trait list and every input file's size and
# mtime are unchanged. Raw-file md5 would cost minutes per run; size + mtime catch replaced
# downloads, and the reader-code hash catches parsing changes. Also keyed (2026-10-04) on the
# HapMap3 list (md5), the EUR .bim used for TDI rsID mapping (size + mtime) and package versions.
fp_inputs <- vapply(ses_all, function(k) {
  fi <- file.info(registry[[k]]$file)
  paste(k, basename(registry[[k]]$file), fi$size, format(fi$mtime, "%Y-%m-%dT%H:%M:%S"), sep = "|")
}, character(1))
bim_fi <- file.info(paste0(BFILE_EUR, ".bim"))
fingerprint <- c(
  paste0("reader_md5=", unname(tools::md5sum(here::here("R", "utils", "sumstats_registry.R")))),
  "munge=info.filter:0.9|maf.filter:0.01",
  paste0("hm3_md5=", unname(tools::md5sum(HM3))),
  paste0("bim=", bim_fi$size, "|", format(bim_fi$mtime, "%Y-%m-%dT%H:%M:%S")),
  paste0("GenomicSEM=", packageVersion("GenomicSEM"), "|data.table=", packageVersion("data.table")),
  paste0("input=", fp_inputs))
fp_file <- file.path(OUT_LDSC, "munge_fingerprint.txt")
cache_ok <- all(file.exists(munged)) && file.exists(fp_file) && identical(readLines(fp_file), fingerprint)
logl("LDSC munge cache: {if (cache_ok) 'valid, reusing' else 'missing or stale, rebuilding'}")
if (!cache_ok) {
  file.remove(c(munged[file.exists(munged)], fp_file[file.exists(fp_file)]))
  inputs <- character(0)
  for (k in ses_all) {
    d <- read_registry(k, p_max = 1, rsids = hm3$SNP, bim = bim[rsid %in% hm3$SNP])
    # signs already oriented (higher = higher SES), so rg signs are interpretable directly
    out <- d[, .(SNP = rsid, A1 = ea, A2 = oa, BETA = beta, P = pval, N = n,
                 MAF = ifelse(is.na(eaf), NA_real_, pmin(eaf, 1 - eaf)))]
    if (all(is.na(out$MAF))) out[, MAF := NULL]
    f <- file.path(OUT_LDSC, paste0(short[k], "_input.tsv"))
    fwrite(out, f, sep = "\t", na = "NA")
    inputs <- c(inputs, f)
    logl("ldsc input {k}: {nrow(out)} HapMap3 SNPs")
  }
  owd <- setwd(OUT_LDSC)
  GenomicSEM::munge(files = inputs, hm3 = HM3, trait.names = unname(short[ses_all]),
                    N = rep(NA, length(ses_all)), info.filter = 0.9, maf.filter = 0.01,
                    log.name = "munge_ses")
  setwd(owd)
  file.remove(inputs)
  writeLines(fingerprint, fp_file)
}

owd <- setwd(OUT_LDSC)
ld_res <- GenomicSEM::ldsc(traits = munged, sample.prev = rep(NA, length(ses_all)),
                           population.prev = rep(NA, length(ses_all)), ld = LD_DIR, wld = LD_DIR,
                           trait.names = unname(short[ses_all]), stand = TRUE, ldsc.log = "ldsc_ses")
setwd(owd)

k <- length(ses_all); nm <- unname(short[ses_all])
rg <- ld_res$S_Stand; dimnames(rg) <- list(nm, nm)
# V_Stand is the sampling covariance of vech(S_Stand) (lower triangle, column-wise)
idx <- which(lower.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE)
se_rg <- matrix(NA_real_, k, k, dimnames = list(nm, nm))
se_v <- sqrt(diag(ld_res$V_Stand))
for (r in seq_len(nrow(idx))) { i <- idx[r, 1]; j <- idx[r, 2]; se_rg[i, j] <- se_rg[j, i] <- se_v[r] }
h2 <- diag(ld_res$S); names(h2) <- nm
I <- ld_res$I; dimnames(I) <- list(nm, nm)

pairs <- t(combn(nm, 2))
rg_tab <- data.table(trait_1 = pairs[, 1], trait_2 = pairs[, 2],
                     rg = round(rg[pairs], 3), se = round(se_rg[pairs], 3),
                     ci_low = round(rg[pairs] - 1.96 * se_rg[pairs], 3),
                     ci_high = round(rg[pairs] + 1.96 * se_rg[pairs], 3),
                     cross_trait_intercept = round(I[pairs], 3))
h2_tab <- data.table(trait = nm, h2_snp = round(h2, 4), ldsc_intercept = round(diag(I), 3))
write_result_table(rg_tab, "feasibility_rg_ldsc")
write_result_table(cbind(data.table(trait = nm), as.data.table(round(I, 4))), "feasibility_ldsc_intercepts")
logl("h2 / intercepts:\n{paste(capture.output(print(h2_tab)), collapse = '\n')}")
logl("rg:\n{paste(capture.output(print(rg_tab)), collapse = '\n')}")

# ---- T2.12: conditional F-statistics ----------------------------------------------

inst <- rbindlist(lapply(ses_all, function(key) {
  x <- fread(file.path(DIR_INSTRUMENTS, glue("{key}_instruments.tsv")), colClasses = list(character = "chr"))
  x[, dataset := key]
}))

is_pal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C")

full_cache <- list()
get_full <- function(key, rsids) {
  d <- read_registry(key, p_max = 1, rsids = rsids, bim = bim[rsid %in% rsids])
  d
}

run_model <- function(keys) {
  # union of the exposures' instruments, re-clumped jointly on the minimum p
  u <- inst[dataset %in% keys, .(pval = min(pval)), by = rsid]
  cl <- ieugwasr::ld_clump(data.frame(rsid = u$rsid, pval = u$pval, id = "joint"),
                           clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, clump_p = 1,
                           bfile = BFILE_EUR, plink_bin = PLINK_BIN)
  snps <- cl$rsid
  ref <- get_full(keys[1], snps)[, .(rsid, ref_ea = ea, ref_oa = oa, ref_eaf = eaf)]
  B <- copy(ref)
  for (key in keys) {
    d <- get_full(key, snps)
    m <- merge(ref, d[, .(rsid, ea, oa, beta, se)], by = "rsid")
    m[, aligned := fifelse(ea == ref_ea & oa == ref_oa, beta,
                    fifelse(ea == ref_oa & oa == ref_ea, -beta, NA_real_))]
    B <- merge(B, m[!is.na(aligned), .(rsid, b = aligned, s = se)], by = "rsid")
    setnames(B, c("b", "s"), c(paste0("b_", key), paste0("se_", key)))
  }
  # drop ambiguous palindromes (reference EAF; 1000G-free check is enough here)
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
  # off-diagonal correlations scaled; scale 0 is the independent-sample assumption
  res <- data.table(model = paste(short[keys], collapse = " + "), exposure = short[keys], n_snps = nrow(B))
  for (sc in PCOR_SCALES) {
    pc <- pcor * sc; diag(pc) <- 1
    if (min(eigen(pc, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      res[, (paste0("condF_pcor_x", sc)) := NA_real_]; next
    }
    f <- suppressWarnings(if (sc == 0) MVMR::strength_mvmr(r_input, gencov = 0)
                          else MVMR::strength_mvmr(r_input, gencov = MVMR::phenocov_mvmr(pc, sx)))
    res[, (paste0("condF_pcor_x", sc)) := as.numeric(f[1, ])]
  }
  # unrounded values for the tolerance check in 02c; the reported table is rounded below
  fwrite(res, file.path(DIR_CONDF, paste0(model_id, "_condF_unrounded.tsv")), sep = "\t")
  fc <- grep("^condF_", names(res), value = TRUE)
  res[, (fc) := lapply(.SD, round, 2), .SDcols = fc]
  setnames(res, c("condF_pcor_x0", "condF_pcor_x1"), c("condF_independent", "condF_overlap_adjusted"))
  res
}

DIR_CONDF <- file.path(DIR_TABLES, "condF_inputs"); fs::dir_create(DIR_CONDF)
PCOR_SCALES <- c(0, 0.5, 1, 1.5)
models <- list(ses_primary,
               c("ea_okbay2022", "income_kweon2025"), c("ea_okbay2022", "tdi_panukb"),
               c("ea_okbay2022", "occupation_akimova2024"),
               c("ea_okbay2022", "income_hill2019", "tdi_panukb", "occupation_akimova2024"),
               # EUR-only TDI sensitivity
               c("ea_okbay2022", "tdi_panukb_eur"),
               c("ea_okbay2022", "income_kweon2025", "tdi_panukb_eur", "occupation_akimova2024"),
               # Le-comparable single-adjustment models (added 2026-10-06; plan v0.6 §10.2 B2): project-rule SNP sets
               c("ea_okbay2022", "income_hill2019"), c("ea_okbay2022", "occupation_akimova2024_isei"))
cf <- rbindlist(lapply(models, function(m) { logl("conditional F: {paste(short[m], collapse = ' + ')}"); run_model(m) }))
cf[, below_10_overlap_adjusted := condF_overlap_adjusted < 10]
write_result_table(cf, "feasibility_conditional_F")
logl("conditional F:\n{paste(capture.output(print(cf)), collapse = '\n')}")
logl("Done.")
close_run_log(log_con)
