#' R/utils/sumstats_registry.R
#'
#' Shared dataset registry and readers for exposure / mediator summary statistics.
#' Sourced by 02_select_instruments.R and 02b_ses_feasibility.R so that column
#' maps, allele conventions and sign flips are defined once.
#'
#' Requires: DIR_DATA_RAW (from 00_setup.R), data.table. `read_panukb()` also
#' needs a 1000G EUR .bim table (columns chr, rsid, cm, pos, a1, a2) to map
#' Pan-UKB chr:pos:alleles to rsIDs.

if (!exists("P_THRESH")) P_THRESH <- 5e-8

# ---- Dataset registry ---------------------------------------------------------
# Column maps are written out explicitly (project rule: no guessing).
# `n_default` is used only where the file has no per-SNP N. `n_source` (optional) labels an N that
# is estimated rather than documented; it is carried into INSTRUMENT_SUMMARY.tsv.
raw <- function(...) file.path(DIR_DATA_RAW, ...)
registry <- list(
  # Exposures (primary)
  ea_okbay2022 = list(role = "exposure", label = "Educational attainment (Okbay 2022 EA4 excl-23andMe)",
    file = raw("Okbay 23andMe Additive Exclusion.txt.gz"),
    cols = c(rsid = "rsID", chr = "Chr", pos = "BP", ea = "Effect_allele", oa = "Other_allele",
             eaf = "EAF_HRC", beta = "Beta", se = "SE", pval = "P"),
    n_default = 765283, units = "SD of educational attainment per allele (SE*sqrt(2pqN) ~ 1.07; not years of schooling)", sign = 1),
  income_kweon2025 = list(role = "exposure", label = "Income Factor (Kweon 2025)",
    file = raw("income_kweon2025_incfactor.txt"),
    cols = c(rsid = "SNP", chr = "CHR", pos = "BP", ea = "A1", oa = "A2",
             eaf = "EAF", beta = "BETA", se = "SE", pval = "P", n = "N"),
    n_default = 668288, units = "SD of the Income Factor per allele", sign = 1),
  tdi_panukb = list(role = "exposure", label = "Townsend Deprivation Index (Pan-UKB 189 irnt, meta)",
    file = raw("tdi_panukb_189_irnt.tsv.bgz"), special = "panukb",
    n_default = 439198, units = "SD (IRNT) per allele, sign-flipped: higher = less deprived", sign = -1),
  # EUR-only stratum of the same Pan-UKB file (sensitivity, added 2026-10-03). The file has no N column:
  # median 1/(2pq SE^2) over common high-confidence variants is 451,551 for EUR and 470,012 for the
  # meta; scaling by the documented meta N (439,198 / 470,012) gives an EUR effective N of ~421,945.
  tdi_panukb_eur = list(role = "exposure (sensitivity)", label = "Townsend Deprivation Index (Pan-UKB 189 irnt, EUR only)",
    file = raw("tdi_panukb_189_irnt.tsv.bgz"), special = "panukb", stratum = "EUR",
    n_default = 421945, n_source = "estimated (median 1/(2pq SE^2), calibrated to the documented meta N)",
    units = "SD (IRNT) per allele, sign-flipped: higher = less deprived", sign = -1),
  occupation_akimova2024 = list(role = "exposure", label = "Occupational status (Akimova 2024 CAMSIS)",
    file = raw("occupation_akimova2024_camsis.tsv"),
    cols = c(rsid = "variant_id", chr = "chromosome", pos = "base_pair_location", ea = "effect_allele",
             oa = "other_allele", eaf = "effect_allele_frequency", beta = "beta", se = "standard_error",
             pval = "p_value", n = "n"),
    n_default = 273157, units = "SD of CAMSIS occupational status per allele (SE*sqrt(2pqN) ~ 0.98)", sign = 1),
  # Added 2026-10-06 (user-approved): Akimova 2024 ISEI (GCST90446163), for the Le-comparable EA4 + ISEI model only.
  # The GWAS Catalog meta.yaml says GRCh38, but positions are GRCh37 (rs10 at 7:92383888, as in the CAMSIS file).
  occupation_akimova2024_isei = list(role = "exposure (comparable model)", label = "Occupational status (Akimova 2024 ISEI)",
    file = raw("occupation_akimova2024_isei.tsv"),
    cols = c(rsid = "variant_id", chr = "chromosome", pos = "base_pair_location", ea = "effect_allele",
             oa = "other_allele", eaf = "effect_allele_frequency", beta = "beta", se = "standard_error",
             pval = "p_value", n = "n"),
    n_default = 271769, units = "SD of ISEI occupational status per allele", sign = 1),
  # Added 2026-10-07 (latent-SES feasibility check; adopted as plan v0.6 §7.4 secondary, pending freeze): Xia 2025 gSES, Genomic SEM common factor of EA,
  # household income, occupational prestige and social deprivation (GCST90566708; UKB + SSGAC; GRCh37).
  ses_gses_xia2025 = list(role = "exposure (shared SES factor; secondary)", label = "General SES factor gSES (Xia 2025)",
    file = raw("ses_gses_xia2025.tsv.gz"),
    cols = c(rsid = "rs_id", chr = "chromosome", pos = "base_pair_location", ea = "effect_allele",
             oa = "other_allele", eaf = "effect_allele_frequency", beta = "beta", se = "standard_error",
             pval = "p_value", n = "n"),
    n_default = 947466, units = "per released gSES factor unit per allele (Genomic SEM; normalization not stated in paper, supplement or deposit metadata)", sign = 1),
  # Exposures (sensitivity comparators)
  ea_lee2018 = list(role = "exposure (sensitivity)", label = "Educational attainment (Lee 2018 EA3 excl-23andMe)",
    file = raw("Lee GWAS EA Excl 23andMe.txt"),
    cols = c(rsid = "MarkerName", chr = "CHR", pos = "POS", ea = "A1", oa = "A2",
             eaf = "EAF", beta = "Beta", se = "SE", pval = "Pval"),
    n_default = 766345, units = "SD of educational attainment per allele (SE*sqrt(2pqN) ~ 1.05)", sign = 1),
  income_hill2019 = list(role = "exposure (sensitivity)", label = "Household income (Hill 2019, UKB only)",
    file = raw("Hill Data", "Household_Income_UKBiobank.txt.gz"),
    cols = c(rsid = "SNP", chr = "Chr", pos = "BPos", ea = "Effect_Allele", oa = "Non_effect_Allele",
             beta = "Beta", se = "Standard_Error_of_Beta", pval = "P"),
    n_default = 286301, units = "standardized log household income per allele", sign = 1),
  # Mediators
  bmi_giant_yengo2018 = list(role = "mediator", label = "BMI (Yengo 2018)",
    file = raw("bmi_giant_yengo2018_combined.txt.gz"), strip_rsid_suffix = TRUE,
    cols = c(rsid = "SNP", chr = "CHR", pos = "POS", ea = "Tested_Allele", oa = "Other_Allele",
             eaf = "Freq_Tested_Allele", beta = "BETA", se = "SE", pval = "P", n = "N"),
    n_default = 681275, units = "SD BMI per allele", sign = 1),
  smoking_wootton2020 = list(role = "mediator", label = "Lifetime smoking index (Wootton 2020)",
    file = raw("Wootton Data", "2019.10.02 Lifetime Smoking GWAS Data Sheet 1.txt"),
    cols = c(rsid = "SNP", chr = "CHR", pos = "BP", ea = "EFFECT_ALLELE", oa = "OTHER_ALLELE",
             eaf = "EAF", beta = "BETA", se = "SE", pval = "P"),
    n_default = 462690, units = "SD lifetime smoking index per allele", sign = 1),
  ldl_glgc_graham2021 = list(role = "mediator", label = "LDL-C (Graham 2021, trans-ancestry)",
    file = raw("ldl_glgc_graham2021_trans_ancestry.tsv.gz"),
    cols = c(rsid = "rsID", chr = "CHROM", pos = "POS_b37", ea = "ALT", oa = "REF", eaf = "POOLED_ALT_AF",
             beta = "METAL_Effect", se = "METAL_StdErr", pval = "METAL_Pvalue", n = "N"),
    n_default = NA, units = "SD (inverse-normal) per ALT allele", sign = 1),
  hdl_glgc_graham2021 = list(role = "mediator", label = "HDL-C (Graham 2021, trans-ancestry)",
    file = raw("hdl_glgc_graham2021_trans_ancestry.tsv.gz"),
    cols = c(rsid = "rsID", chr = "CHROM", pos = "POS_b37", ea = "ALT", oa = "REF", eaf = "POOLED_ALT_AF",
             beta = "METAL_Effect", se = "METAL_StdErr", pval = "METAL_Pvalue", n = "N"),
    n_default = NA, units = "SD (inverse-normal) per ALT allele", sign = 1),
  tg_glgc_graham2021 = list(role = "mediator", label = "Triglycerides (Graham 2021, trans-ancestry)",
    file = raw("tg_glgc_graham2021_trans_ancestry.tsv.gz"),
    cols = c(rsid = "rsID", chr = "CHROM", pos = "POS_b37", ea = "ALT", oa = "REF", eaf = "POOLED_ALT_AF",
             beta = "METAL_Effect", se = "METAL_StdErr", pval = "METAL_Pvalue", n = "N"),
    n_default = NA, units = "SD (inverse-normal log TG) per ALT allele", sign = 1),
  sbp_keaton2024 = list(role = "mediator", label = "Systolic blood pressure (Keaton 2024)",
    file = raw("sbp_keaton2024.tsv.gz"),
    cols = c(rsid = "rs_id", chr = "chromosome", pos = "base_pair_location", ea = "effect_allele",
             oa = "other_allele", eaf = "effect_allele_frequency", beta = "beta", se = "standard_error",
             pval = "p_value", n = "n"),
    n_default = NA, units = "mmHg per allele", sign = 1),
  t2d_diamante_mahajan2022 = list(role = "mediator", label = "Type 2 diabetes (DIAMANTE 2022, trans-ancestry)",
    file = raw("DIAMANTE-TA.sumstat.txt"),
    cols = c(rsid = "rsID", chr = "chromosome(b37)", pos = "position(b37)", ea = "effect_allele",
             oa = "other_allele", beta = "Fixed-effects_beta", se = "Fixed-effects_SE",
             pval = "Fixed-effects_p-value", n = "Effective_sample_size"),
    n_default = NA, units = "log-OR per allele (fixed-effects); R2 is z-based, not liability-scale", sign = 1)
)

# ---- Readers -----------------------------------------------------------------

read_standard <- function(cfg) {
  cm <- cfg$cols
  d <- fread(cfg$file, select = unname(cm), showProgress = FALSE,
             colClasses = list(character = unname(cm[intersect(c("rsid", "chr", "ea", "oa"), names(cm))])))
  setnames(d, unname(cm), names(cm))
  if (isTRUE(cfg$strip_rsid_suffix)) d[, rsid := sub(":.*$", "", rsid)]
  d[, pval := as.numeric(pval)]
  d
}

read_panukb <- function(cfg, p_max = P_THRESH, bim) {
  st <- if (is.null(cfg$stratum)) "meta" else cfg$stratum
  cols <- c(af = paste0("af_", st), beta = paste0("beta_", st), se = paste0("se_", st),
            nlp = paste0("neglog10_pval_", st))
  d <- fread(cfg$file, select = c("chr", "pos", "ref", "alt", unname(cols), "low_confidence_EUR"),
             colClasses = list(character = c("chr", "ref", "alt")), showProgress = FALSE)
  setnames(d, unname(cols), names(cols))
  # Pan-UKB: beta is for the alt allele; drop variants flagged low-confidence in EUR (the dominant stratum)
  d <- d[is.na(low_confidence_EUR) | low_confidence_EUR != "true"]
  d[, `:=`(pval = 10^(-nlp), ea = alt, oa = ref, eaf = af)]
  d <- d[!is.na(pval) & !is.na(beta) & pval < p_max]
  # map to rsID through the EUR bim on chr + pos + allele set
  m <- merge(d, bim[, .(chr, pos, rsid, a1, a2)], by = c("chr", "pos"))
  m <- m[(ea == a1 & oa == a2) | (ea == a2 & oa == a1)]
  m[, .(rsid, chr, pos, ea, oa, eaf, beta, se, pval)]
}

std <- function(d, cfg) {
  if (!"n" %in% names(d)) d[, n := cfg$n_default]
  if (!"eaf" %in% names(d)) d[, eaf := NA_real_]
  d[, `:=`(ea = toupper(ea), oa = toupper(oa), chr = as.character(chr), pos = as.integer(pos),
           beta = as.numeric(beta) * cfg$sign, se = as.numeric(se), eaf = as.numeric(eaf), n = as.numeric(n))]
  d[, .(rsid, chr, pos, ea, oa, eaf, beta, se, pval, n)]
}

# EUR frequency of the effect allele from the reference panel (for files without EAF)
eur_freq <- function(rsids) {
  tf <- tempfile(); writeLines(unique(rsids), paste0(tf, ".snps"))
  system2(PLINK_BIN, c("--bfile", shQuote(BFILE_EUR), "--extract", shQuote(paste0(tf, ".snps")),
                       "--freq", "--out", shQuote(tf), "--silent"), stdout = FALSE, stderr = FALSE)
  fr <- fread(paste0(tf, ".frq"))
  fr[, .(rsid = SNP, ref_a1 = A1, ref_a2 = A2, ref_maf = MAF)]
}

is_palindromic <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C")

#' Read one registry dataset as a standardized data.table
#' (rsid, chr, pos, ea, oa, eaf, beta, se, pval, n), sign applied.
#' p_max filters on p-value (1 = keep all); rsids optionally restricts rows.
read_registry <- function(key, p_max = 1, rsids = NULL, bim = NULL) {
  cfg <- registry[[key]]
  d <- if (identical(cfg$special, "panukb")) read_panukb(cfg, p_max = p_max, bim = bim) else {
    x <- read_standard(cfg); x[!is.na(pval) & pval < p_max]
  }
  d <- std(d, cfg)
  d <- d[!is.na(beta) & !is.na(se) & se > 0 & grepl("^rs[0-9]+$", rsid)]
  if (!is.null(rsids)) d <- d[rsid %in% rsids]
  unique(d, by = "rsid")
}
