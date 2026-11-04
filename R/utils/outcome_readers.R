#' outcome_readers.R
#'
#' Outcome-file readers for harmonization (moved unchanged from R/03_harmonize.R on 2026-10-08 so that R/03d can reuse
#' them for proxy SNPs). make_outcome_readers(snp_ref, RAW, logl) returns readers bound to a reference table
#' snp_ref (rsid, chr, pos [GRCh37], a1, a2): rule-1 matching key, rule-3 duplicates, MVP SE from the 95% CI (rule 7),
#' van Zuydam rule 11. Each reader returns one row per snp_ref rsid with status ok | missing | allele_mismatch | duplicate.

make_outcome_readers <- function(snp_ref, RAW, logl) {
  allele_pair_ok <- function(o1, o2, a1, a2) (o1 == a1 & o2 == a2) | (o1 == a2 & o2 == a1)

  # ---- 2. outcome extraction -------------------------------------------------------------------------------------
  # Each reader returns: rsid, ea_y, oa_y, eaf_y, beta_y, se_y, p_y, (se_y_sens), (r2), status
  # status: ok | missing | allele_mismatch | duplicate
  finish <- function(m, key_rows) {
    # m: candidate rows (rsid, alleles...) already restricted to the key; key_rows: rsids with a row at the key
    m <- m[allele_pair_ok(ea_y, oa_y, a1, a2)]
    dup <- m[, .N, by = rsid][N > 1, rsid]
    ok <- m[!rsid %in% dup]
    ok[, intersect(c("chr", "pos", "a1", "a2"), names(ok)) := NULL]
    res <- merge(snp_ref[, .(rsid)], ok, by = "rsid", all.x = TRUE)
    res[, status := fifelse(!is.na(beta_y), "ok", fifelse(rsid %in% dup, "duplicate",
                      fifelse(rsid %in% key_rows, "allele_mismatch", "missing")))]
    res
  }
  by_pos <- function(d) {  # d: chr, pos, ea_y, oa_y + value columns (GRCh37)
    d[, `:=`(chr = sub("^chr", "", as.character(chr)), pos = as.integer(pos), ea_y = toupper(ea_y), oa_y = toupper(oa_y))]
    m <- merge(snp_ref, d, by = c("chr", "pos"))
    finish(m, unique(m$rsid))
  }
  by_rsid <- function(d) {  # d: rsid + alleles + values (GRCh38 files)
    d[, `:=`(ea_y = toupper(ea_y), oa_y = toupper(oa_y))]
    m <- merge(snp_ref, d, by = "rsid")
    finish(m, unique(m$rsid))
  }
  ci_se <- function(lo, hi) (log(as.numeric(hi)) - log(as.numeric(lo))) / (2 * qnorm(0.975))
  gwas_ssf <- function(f, extra = character()) {
    d <- fread(file.path(RAW, f), select = c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                                            "effect_allele_frequency", "beta", "standard_error", "p_value", extra), showProgress = FALSE)
    setnames(d, c("chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta",
                  "standard_error", "p_value"), c("chr", "pos", "ea_y", "oa_y", "eaf_y", "beta_y", "se_y", "p_y"))
    d
  }
  mvp_dbgap <- function(f) {
    d <- fread(file.path(RAW, f), select = c("SNP_ID", "ref", "alt", "ea", "af", "or", "ci", "pval", "r2"), showProgress = FALSE)
    d <- d[SNP_ID %in% snp_ref$rsid]
    stopifnot(all(d$ea == d$alt))
    lims <- tstrsplit(d$ci, ",", fixed = TRUE)
    by_rsid(d[, .(rsid = SNP_ID, ea_y = alt, oa_y = ref, eaf_y = as.numeric(af), beta_y = log(as.numeric(or)),
                  se_y = ci_se(lims[[1]], lims[[2]]), p_y = as.numeric(pval), r2 = as.numeric(r2))])
  }
  mvp_gwascat <- function(f) {
    d <- fread(file.path(RAW, f), select = c("rsid", "effect_allele", "other_allele", "effect_allele_frequency", "odds_ratio",
                                            "ci_lower", "ci_upper", "p_value", "r2"), showProgress = FALSE, na.strings = c("NA", "#NA", ""))
    d <- d[rsid %in% snp_ref$rsid]
    by_rsid(d[, .(rsid, ea_y = effect_allele, oa_y = other_allele, eaf_y = as.numeric(effect_allele_frequency),
                  beta_y = log(as.numeric(odds_ratio)), se_y = ci_se(ci_lower, ci_upper), p_y = as.numeric(p_value),
                  r2 = as.numeric(gsub("\r", "", r2)))])
  }
  outcomes <- list(
    cad_aragam_eur = function() by_pos(gwas_ssf("cad_aragam2022_eur.tsv")),
    cad_aragam_multi = function() by_pos(gwas_ssf("cad_aragam2022_multi.tsv")),
    cad_mvp_phe411_eur = function() mvp_dbgap("cad_mvp_phe411_EUR.txt.gz"),
    pad_mvp_eur = function() mvp_dbgap("pad_mvp_phe440_2_EUR.txt.gz"),
    pad_mvp_afr = function() mvp_dbgap("pad_mvp_phe440_2_AFR.txt.gz"),
    pad_finngen_athscle = function() {
      d <- fread(file.path(RAW, "pad_finngen_r12_athscle.tsv.gz"), select = c("ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"),
                 showProgress = FALSE)
      multi <- grepl(",", d$rsids, fixed = TRUE)
      d1 <- d[!multi & rsids %in% snp_ref$rsid][, rsid := rsids]
      d2 <- d[multi][, .(rsid = unlist(strsplit(rsids, ",", fixed = TRUE))), by = .(ref, alt, pval, beta, sebeta, af_alt)][rsid %in% snp_ref$rsid]
      d <- rbind(d1[, .(rsid, ref, alt, pval, beta, sebeta, af_alt)], d2[, .(rsid, ref, alt, pval, beta, sebeta, af_alt)])
      by_rsid(d[, .(rsid, ea_y = alt, oa_y = ref, eaf_y = af_alt, beta_y = beta, se_y = sebeta, p_y = pval)])
    },
    pad_vanzuydam = function() {
      d <- fread(file.path(RAW, "pad_vanzuydam2021", "final_meta", "pad_primary.txt.gz"), showProgress = FALSE)
      stopifnot(all(c("chr", "bp", "or_allele", "nonor_allele", "eaf", "or", "or_se", "p") %in% names(d)))
      d <- d[, .(chr, pos = bp, ea_y = or_allele, oa_y = nonor_allele, eaf_y = as.numeric(eaf), or = as.numeric(or),
                 or_se = as.numeric(or_se), p_y = as.numeric(p))]
      d[, `:=`(beta_y = log(or), se_y = or_se / or)]
      z <- qnorm(1 - d$p_y / 2)
      d[, se_y_sens := pmax(se_y, fifelse(z > 0, abs(beta_y) / z, NA_real_), na.rm = TRUE)]  # rule 11 sensitivity
      by_pos(d[, .(chr, pos, ea_y, oa_y, eaf_y, beta_y, se_y, p_y, se_y_sens, or_rounded_1 = or == 1)])
    },
    hf_hermes2 = function() {
      d <- gwas_ssf("hf_hermes2_2025_hfall.tsv.gz", extra = "rsid")
      r <- by_pos(d[, .(chr, pos, ea_y, oa_y, eaf_y, beta_y, se_y, p_y, rsid_y = rsid)])
      logl("  HERMES rsID cross-check: {sum(r$status == 'ok' & r$rsid != r$rsid_y, na.rm = TRUE)} matched rows carry a different rsID label")
      r[, rsid_y := NULL]
    },
    hf_mvp_phe428_eur = function() mvp_gwascat("hf_mvp_phe428_EUR.gwascat.tsv.gz"),
    stroke_gigastroke_multi = function() by_pos(gwas_ssf("stroke_gigastroke_any_multi.tsv.gz")),
    stroke_gigastroke_eur = function() by_pos(gwas_ssf("stroke_gigastroke_any_eur.tsv.gz")),
    stroke_mvp_phe433_21_eur = function() mvp_gwascat("stroke_mvp_phe433_21_EUR.gwascat.tsv.gz")
  )
  list(outcomes = outcomes, by_rsid = by_rsid, allele_pair_ok = allele_pair_ok)
}
