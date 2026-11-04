#' 03d_ld_proxies.R
#'
#' Purpose: inputs for two pre-specified sensitivity analyses (plan v1.0 §6; harmonization rules 5 and 9b). No MR
#'   estimate is computed.
#'   A. Residual LD (rule 9b): plink --r2 within each instrument set (1000G EUR; r2 >= 0.001 within 10 Mb). The weaker
#'      SNP of each pair is dropped, using the rule of the pre-registered Li D pruned set (R/02e): pairs taken in
#'      descending r2, and if both SNPs are still present the weaker is dropped. Strength = |beta/se| (univariable
#'      sets) or the maximum |b/se| across exposures (MVMR sets). The script stops unless this rule reproduces the
#'      R/02e Li D pruned set exactly.
#'   B. Proxies (rule 5): for each outcome, every union instrument SNP without a usable outcome row (missing, allele
#'      mismatch or duplicate) gets the highest-r2 proxy present in the outcome, with r2 >= 0.8 within 500 kb in 1000G
#'      EUR (MVP AFR: r2 >= 0.8 in both EUR and AFR, with the same allele phase in both). Palindromic proxies are not
#'      used (strand cannot be checked). The outcome effect is mapped to the target SNP's alleles through the PLINK
#'      in-phase allele pairs; the proxy row's beta, SE and EAF are carried. Proxy-augmented cells are harmonized with
#'      the same code as R/03 (R/utils/harmonize_utils.R).
#' Outputs: data/instruments/ld_pruned/{set}.tsv (kept SNPs) ; results/tables/residual_ld_audit.tsv
#'   data/harmonized/large/{univariable_proxy,mvmr_proxy}/{set}__{outcome}.tsv.gz
#'   results/tables/proxy_log.tsv ; results/tables/harmonization_retention_proxy.tsv
#' Depends on: data.table, plinkbinr

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
log_con <- log_run("03d_ld_proxies")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
source(here::here("R", "utils", "sumstats_registry.R"))
source(here::here("R", "utils", "outcome_readers.R"))
source(here::here("R", "utils", "harmonize_utils.R"))
PLINK_BIN <- plinkbinr::get_plink_exe()
BF <- list(EUR = here::here("data", "reference", "1kg_v3", "EUR"), AFR = here::here("data", "reference", "1kg_v3", "AFR"))
OUT <- file.path(DIR_DATA_HARM, "large"); RAW <- DIR_DATA_RAW
for (d in c("univariable_proxy", "mvmr_proxy")) fs::dir_create(file.path(OUT, d))
fs::dir_create(file.path(DIR_INSTRUMENTS, "ld_pruned"))
tmp <- tempfile("ld"); fs::dir_create(tmp)

plink_r2 <- function(snps, bfile, extra, tag) {
  f <- file.path(tmp, tag); writeLines(unique(snps), paste0(f, ".snps"))
  rc <- system2(PLINK_BIN, c("--bfile", shQuote(bfile), extra, "--out", shQuote(f), "--silent"), stdout = FALSE, stderr = FALSE)
  if (!file.exists(paste0(f, ".ld"))) return(data.table())
  fread(paste0(f, ".ld"))
}

# ---- inputs: instrument sets (same as R/03) ---------------------------------------------------------------------
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
inst <- lapply(c(uni_keys, med_keys), read_inst); names(inst) <- c(names(uni_keys), names(med_keys))
bfiles <- c(list.files(file.path(DIR_TABLES, "condF_inputs"), "_B\\.tsv$", full.names = TRUE),
            list.files(file.path(DIR_TABLES, "condF_inputs_mediators"), "_B\\.tsv$", full.names = TRUE),
            file.path(DIR_TABLES, "pilot_lid_cad_inputs", c("S1_B.tsv", "S1-LDpruned_B.tsv")))
Bs <- lapply(bfiles, fread); names(Bs) <- sub("_B\\.tsv$", "", basename(bfiles))
names(Bs)[names(Bs) %in% c("S1", "S1-LDpruned")] <- paste0("LiD_", names(Bs)[names(Bs) %in% c("S1", "S1-LDpruned")])
snp_ref <- unique(rbindlist(inst)[, .(rsid, chr, pos, a1 = ea_x, a2 = oa_x)], by = "rsid")

# ---- A. residual LD ---------------------------------------------------------------------------------------------
prune <- function(snps, strength, tag) {
  ld <- plink_r2(snps, BF$EUR, c("--extract", shQuote(file.path(tmp, paste0(tag, ".snps"))), "--r2", "--ld-window-kb", "10000",
                                 "--ld-window", "99999", "--ld-window-r2", "0.001"), tag)
  pairs <- if (nrow(ld)) ld[SNP_A != SNP_B & R2 >= 0.001, .(a = SNP_A, b = SNP_B, r2 = R2)] else data.table(a = character(), b = character(), r2 = numeric())
  n_pairs <- nrow(pairs); max_r2 <- if (n_pairs) max(pairs$r2) else NA_real_
  # Same rule as the pre-registered Li D pruned set (R/02e): pairs in descending r2; if both SNPs are still present,
  # drop the weaker one (smaller strength); ties: drop the lexicographically larger rsID.
  pairs <- pairs[order(-r2)]; dropped <- character()
  for (k in seq_len(nrow(pairs))) {
    a <- pairs$a[k]; b <- pairs$b[k]
    if (a %in% dropped || b %in% dropped) next
    d <- if (strength[a] < strength[b]) a else if (strength[b] < strength[a]) b else max(a, b)
    dropped <- c(dropped, d)
  }
  list(dropped = dropped, n_pairs = n_pairs, max_r2 = max_r2)
}
aud <- list()
for (s in names(uni_keys)) {
  x <- inst[[s]]; st <- setNames(abs(x$beta_x / x$se_x), x$rsid)
  writeLines(x$rsid, file.path(tmp, paste0("u_", s, ".snps")))
  p <- prune(x$rsid, st, paste0("u_", s))
  fwrite(x[!rsid %in% p$dropped, .(rsid)], file.path(DIR_INSTRUMENTS, "ld_pruned", paste0(s, ".tsv")), sep = "\t")
  aud[[length(aud) + 1]] <- data.table(type = "univariable", set = s, n = nrow(x), n_pairs_r2_ge_0.001 = p$n_pairs, max_r2 = p$max_r2,
                                       n_dropped = length(p$dropped), n_kept = nrow(x) - length(p$dropped))
}
for (b in setdiff(names(Bs), "LiD_S1-LDpruned")) {
  B <- Bs[[b]]; bc <- grep("^b_", names(B), value = TRUE)
  z <- sapply(bc, function(k) abs(B[[k]] / B[[sub("^b_", "se_", k)]])); st <- setNames(apply(as.matrix(z), 1, max), B$rsid)
  writeLines(B$rsid, file.path(tmp, paste0("m_", b, ".snps")))
  p <- prune(B$rsid, st, paste0("m_", b))
  fwrite(B[!rsid %in% p$dropped, .(rsid)], file.path(DIR_INSTRUMENTS, "ld_pruned", paste0("mvmr_", b, ".tsv")), sep = "\t")
  aud[[length(aud) + 1]] <- data.table(type = "mvmr", set = b, n = nrow(B), n_pairs_r2_ge_0.001 = p$n_pairs, max_r2 = p$max_r2,
                                       n_dropped = length(p$dropped), n_kept = nrow(B) - length(p$dropped))
}
A <- rbindlist(aud); fwrite(A, file.path(DIR_TABLES, "residual_ld_audit.tsv"), sep = "\t")
lid_pre <- nrow(Bs[["LiD_S1"]]) - nrow(Bs[["LiD_S1-LDpruned"]])
lid_kept <- fread(file.path(DIR_INSTRUMENTS, "ld_pruned", "mvmr_LiD_S1.tsv"))$rsid
if (!setequal(lid_kept, Bs[["LiD_S1-LDpruned"]]$rsid)) stop("R/03d pruning rule does not reproduce the pre-registered Li D pruned set")
logl("residual LD: {sum(A$n_pairs_r2_ge_0.001 > 0)} of {nrow(A)} sets have pairs; Li D S1 rule here drops {A[set == 'LiD_S1', n_dropped]} vs pre-registered R/02e set {lid_pre} (kept sets identical)")

# ---- B. proxies -------------------------------------------------------------------------------------------------
ext <- lapply(setNames(nm = sub("\\.tsv\\.gz$", "", list.files(file.path(OUT, "outcome_extracts")))),
              function(o) fread(file.path(OUT, "outcome_extracts", paste0(o, ".tsv.gz")), colClasses = list(character = c("ea_y", "oa_y"))))
need <- unique(unlist(lapply(ext, function(e) e[status != "ok", rsid])))
logl("SNPs needing a proxy in at least one outcome: {length(need)}")
writeLines(need, file.path(tmp, "need.snps"))
ldargs <- c("--r2", "in-phase", "with-freqs", "--ld-snp-list", shQuote(file.path(tmp, "need.snps")), "--ld-window-kb", "500",
            "--ld-window", "99999", "--ld-window-r2", "0.8")
cand <- plink_r2(need, BF$EUR, ldargs, "prox_eur")[SNP_A != SNP_B]
cand_afr <- plink_r2(need, BF$AFR, ldargs, "prox_afr")[SNP_A != SNP_B]
bim <- fread(paste0(BF$EUR, ".bim"), col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"), colClasses = list(character = c(1, 2, 5, 6)))
cand <- merge(cand, bim[, .(SNP_B = rsid, b1 = a1, b2 = a2)], by = "SNP_B")
cand <- cand[!is_palindromic(b1, b2)]
ph <- tstrsplit(cand$PHASE, "/", fixed = TRUE)  # e.g. "TT/CA": target T ~ proxy T; target C ~ proxy A
cand[, `:=`(t1 = substr(ph[[1]], 1, 1), p1 = substr(ph[[1]], 2, 2), t2 = substr(ph[[2]], 1, 1), p2 = substr(ph[[2]], 2, 2))]
stopifnot(all(cand$p1 != cand$p2 & cand$t1 != cand$t2))
afr <- cand_afr[, .(SNP_A, SNP_B, R2_AFR = R2, PHASE_AFR = PHASE)]
cand <- merge(cand, afr, by = c("SNP_A", "SNP_B"), all.x = TRUE)
proxy_ref <- unique(cand[, .(rsid = SNP_B, chr = as.character(CHR_B), pos = as.integer(BP_B), a1 = b1, a2 = b2)], by = "rsid")
rd <- make_outcome_readers(proxy_ref, RAW, logl)
plog <- list(); pext <- list()
for (o in names(ext)) {
  pe <- rd$outcomes[[o]]()[status == "ok"]
  cc <- cand[SNP_A %in% ext[[o]][status != "ok", rsid] & SNP_B %in% pe$rsid]
  if (o == "pad_mvp_afr") cc <- cc[!is.na(R2_AFR) & R2_AFR >= 0.8 & PHASE_AFR == PHASE]
  setorder(cc, SNP_A, -R2); best <- cc[, .SD[1], by = SNP_A]
  m <- merge(best, pe, by.x = "SNP_B", by.y = "rsid")
  # map the proxy's outcome effect allele to the target allele in phase with it
  m[, `:=`(ea_t = fifelse(ea_y == p1, t1, fifelse(ea_y == p2, t2, NA_character_)))]
  m[, oa_t := fifelse(ea_t == t1, t2, t1)]
  stopifnot(all(!is.na(m$ea_t)))
  tgt <- m[, .(rsid = SNP_A, ea_y = ea_t, oa_y = oa_t, eaf_y, beta_y, se_y, p_y, status = "ok", proxy = SNP_B, proxy_r2 = R2,
               proxy_r2_afr = R2_AFR, r2 = if ("r2" %in% names(m)) m$r2 else NA_real_)]
  # the target's alleles must match the instrument's alleles
  chk <- merge(tgt, snp_ref, by = "rsid"); stopifnot(all(rd$allele_pair_ok(chk$ea_y, chk$oa_y, chk$a1, chk$a2)))
  base <- ext[[o]][status == "ok"]; base[, `:=`(proxy = NA_character_, proxy_r2 = NA_real_, proxy_r2_afr = NA_real_)]
  if (!"r2" %in% names(base)) base[, r2 := NA_real_]
  keep_cols <- intersect(names(base), names(tgt))
  pext[[o]] <- rbind(base[, ..keep_cols], tgt[, ..keep_cols])
  plog[[o]] <- tgt[, .(outcome = o, rsid, proxy, proxy_r2, proxy_r2_afr)]
  logl("{o}: {nrow(tgt)} of {sum(ext[[o]]$status != 'ok')} unusable SNPs proxied")
}
fwrite(rbindlist(plog), file.path(DIR_TABLES, "proxy_log.tsv"), sep = "\t")
ret <- list()
tol_for <- function(o) if (o == "pad_vanzuydam") 0.085 else 0.08
for (s in names(uni_keys)) for (o in names(pext)) {
  h <- harmonize(inst[[s]], pext[[o]], tol = tol_for(o), drop_all_pal = o == "pad_mvp_afr")
  fwrite(h, file.path(OUT, "univariable_proxy", paste0(s, "__", o, ".tsv.gz")), sep = "\t")
  r <- retention(h, s, o, "univariable_proxy"); r[, n_final_proxied := sum(h$mr_keep & !is.na(h$proxy))]; ret[[length(ret) + 1]] <- r
}
for (b in names(Bs)) for (o in names(pext)) {
  B <- Bs[[b]]
  x <- B[, c(list(rsid = rsid, ea_x = toupper(ref_ea), oa_x = toupper(ref_oa), eaf_x = as.numeric(ref_eaf)), .SD), .SDcols = patterns("^(b|se)_")]
  x <- merge(x, snp_ref[, .(rsid, chr, pos)], by = "rsid")
  h <- harmonize(x, pext[[o]], tol = tol_for(o), drop_all_pal = o == "pad_mvp_afr")
  fwrite(h, file.path(OUT, "mvmr_proxy", paste0(b, "__", o, ".tsv.gz")), sep = "\t")
  r <- retention(h, b, o, "mvmr_proxy"); r[, n_final_proxied := sum(h$mr_keep & !is.na(h$proxy))]; ret[[length(ret) + 1]] <- r
}
R <- rbindlist(ret); fwrite(R, file.path(DIR_TABLES, "harmonization_retention_proxy.tsv"), sep = "\t")
logl("proxy cells: {nrow(R)}; univariable SNPs added by proxies: {sum(R[type == 'univariable_proxy', n_final_proxied])}")
logl("Done. No MR estimate computed.")
close_run_log(log_con)
