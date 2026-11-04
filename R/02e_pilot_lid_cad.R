#' 02e_pilot_lid_cad.R
#'
#' Purpose: bounded feasibility pilot (spec v2: docs/prereg/pilot_lid_cad_spec.md in the OSF supporting documents). Conditional instrument
#'   strength of MVMR models COMPARABLE to Li D et al. 2025 (Cardiovasc Diabetol 24:241) EA + household income -> CAD,
#'   rebuilt from their stated rules and the same exposure data. Exposure-side strength only: NO causal-effect
#'   estimate is computed, and no outcome beta enters any calculation. The outcome file is used only for SNP
#'   availability and the keep/drop flags of allele harmonization. Sequences are fixed by the spec and must not be
#'   tuned to Li D's 483-SNP count (reported as a check only).
#'
#' Sequences: S1 per-exposure clump -> union -> joint re-clump (min p); S2 union without re-clump (diagnostic, with
#'   LD audit); S3 pooled p<5e-8 SNPs -> single clump (min p); S1-ref = S1 without the CAD availability step.
#'   S1-LDpruned (added 2026-10-06): fixed LD sensitivity, S1 minus one SNP per pair with r2 >= 0.001.
#' Covariance of EA/income estimation errors: rho = c(0, 0.5, 1, 1.5) x 0.230 (LDSC-intercept-implied, from
#'   results/tables/feasibility_ldsc_intercepts.tsv).
#'
#' Inputs:  data/raw/Okbay 23andMe Additive Exclusion.txt.gz; data/raw/Hill Data/Household_Income_UKBiobank.txt.gz;
#'          data/raw/cad_aragam2022_eur.tsv (GCST90132314, md5-verified); data/reference/1kg_v3/EUR.*
#' Outputs: results/tables/pilot_lid_cad_condF.{tsv,xlsx}, pilot_lid_cad_retention.{tsv,xlsx},
#'          pilot_lid_cad_ld_audit.{tsv,xlsx}; results/tables/pilot_lid_cad_inputs/{set}_B.tsv, all_condF_unrounded.tsv; results/logs/{date}_02e_pilot_lid_cad.log (+ sessionInfo)
#' Depends on: data.table, ieugwasr, plinkbinr, MVMR, TwoSampleMR

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
ensure_pkgs_loaded(c("data.table", "ieugwasr", "plinkbinr", "MVMR", "TwoSampleMR"))

P_THRESH <- 5e-8; MAF_MIN <- 0.01; CLUMP_KB <- 10000; CLUMP_R2 <- 0.001; PAL_MAF <- 0.42
BFILE_EUR <- here::here("data", "reference", "1kg_v3", "EUR")
PLINK_BIN <- plinkbinr::get_plink_exe()
DIR_IN <- file.path(DIR_TABLES, "pilot_lid_cad_inputs"); fs::dir_create(DIR_IN)
TMP <- file.path(tempdir(), "pilot"); fs::dir_create(TMP)

log_con <- log_run("02e_pilot_lid_cad")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)

Iq <- fread(file.path(DIR_TABLES, "feasibility_ldsc_intercepts.tsv"))
Im <- as.matrix(Iq[, -1]); dimnames(Im) <- list(Iq$trait, Iq$trait)
RHO_LDSC <- Im["EA", "Income_Hill"] / sqrt(Im["EA", "EA"] * Im["Income_Hill", "Income_Hill"])
RHO <- c(0, 0.5, 1, 1.5) * RHO_LDSC
logl("rho grid: {paste(round(RHO, 4), collapse = ', ')} (LDSC-implied {round(RHO_LDSC, 4)})")

retention <- list(); addr <- function(set, step, n) retention[[length(retention) + 1]] <<- data.table(set = set, step = step, n_snps = n)

# ---- exposures ------------------------------------------------------------------------------------------
bim <- fread(paste0(BFILE_EUR, ".bim"), header = FALSE, col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"),
             colClasses = c("character", "character", "numeric", "integer", "character", "character"))
ea <- fread(here::here("data", "raw", "Okbay 23andMe Additive Exclusion.txt.gz"),
            select = c("rsID", "Chr", "BP", "Effect_allele", "Other_allele", "EAF_HRC", "Beta", "SE", "P"))
setnames(ea, c("rsid", "chr", "pos", "ea", "oa", "eaf", "beta", "se", "p"))
ea <- unique(ea[!is.na(beta) & !is.na(se)], by = "rsid")
hill <- fread(here::here("data", "raw", "Hill Data", "Household_Income_UKBiobank.txt.gz"),
              select = c("SNP", "Chr", "BPos", "Effect_Allele", "Non_effect_Allele", "Beta", "Standard_Error_of_Beta", "P"))
setnames(hill, c("rsid", "chr", "pos", "ea", "oa", "beta", "se", "p"))
hill <- unique(hill[!is.na(beta) & !is.na(se)], by = "rsid")
ea[, `:=`(ea = toupper(ea), oa = toupper(oa))]
hill[, `:=`(ea = toupper(ea), oa = toupper(oa))]

# reference alleles = EA4's where the SNP is in the EA4 file; Hill's otherwise (with 1000G EUR frequency)
gw <- unique(c(ea[p < P_THRESH, rsid], hill[p < P_THRESH, rsid]))
ref <- ea[rsid %in% gw, .(rsid, ref_ea = ea, ref_oa = oa, ref_eaf = eaf, eaf_src = "EA4 EAF_HRC")]
hill_only <- setdiff(hill[rsid %in% gw, rsid], ref$rsid)
if (length(hill_only)) {
  writeLines(hill_only, file.path(TMP, "ho.txt"))
  # paths contain spaces: quote every argument (an unquoted call failed silently in the first run)
  st <- system2(PLINK_BIN, shQuote(c("--bfile", BFILE_EUR, "--extract", file.path(TMP, "ho.txt"), "--freq", "--out", file.path(TMP, "ho"))),
                stdout = FALSE, stderr = FALSE)
  n_in_panel <- sum(hill_only %in% bim$rsid)
  if (n_in_panel > 0 && (st != 0 || !file.exists(file.path(TMP, "ho.frq")))) stop("plink --freq failed (exit ", st, ")")
  frq <- if (file.exists(file.path(TMP, "ho.frq"))) fread(file.path(TMP, "ho.frq")) else data.table(SNP = character(), A1 = character(), A2 = character(), MAF = numeric())
  hr <- merge(hill[rsid %in% hill_only, .(rsid, ref_ea = ea, ref_oa = oa)], frq[, .(rsid = SNP, A1, A2, MAF)], by = "rsid", all.x = TRUE)
  hr[, ref_eaf := fifelse(A1 == ref_ea & A2 == ref_oa, MAF, fifelse(A1 == ref_oa & A2 == ref_ea, 1 - MAF, NA_real_))]
  ref <- rbind(ref, hr[, .(rsid, ref_ea, ref_oa, ref_eaf, eaf_src = "1000G EUR (plink --freq)")])
}
n_nofreq <- ref[is.na(ref_eaf), .N]; ref <- ref[!is.na(ref_eaf)]
logl("GW SNPs (either exposure): {length(gw)}; reference alleles from EA4: {ref[eaf_src == 'EA4 EAF_HRC', .N]}, from Hill+1000G: {ref[eaf_src != 'EA4 EAF_HRC', .N]}; excluded for missing frequency: {n_nofreq}")

align <- function(d, keys) {   # align an exposure table to the reference alleles
  m <- merge(ref[rsid %in% keys], d[, .(rsid, ea, oa, beta, se, p)], by = "rsid")
  m[, b := fifelse(ea == ref_ea & oa == ref_oa, beta, fifelse(ea == ref_oa & oa == ref_ea, -beta, NA_real_))]
  m[!is.na(b), .(rsid, b, se, p)]
}
cand <- function(d, label) {
  x <- merge(d[p < P_THRESH, .(rsid, p)], ref[, .(rsid, ref_eaf)], by = "rsid")
  x <- x[pmin(ref_eaf, 1 - ref_eaf) > MAF_MIN & rsid %in% bim$rsid]
  logl("{label}: {d[p < P_THRESH, .N]} at p<5e-8; {nrow(x)} after MAF>{MAF_MIN} and presence in the 1000G EUR panel")
  x[, .(rsid, p)]
}
clump <- function(x) ieugwasr::ld_clump(data.frame(rsid = x$rsid, pval = x$p, id = "x"), clump_kb = CLUMP_KB,
                                        clump_r2 = CLUMP_R2, clump_p = 1, bfile = BFILE_EUR, plink_bin = PLINK_BIN)$rsid
c_ea <- cand(ea, "EA4"); c_hill <- cand(hill, "Hill income")
k_ea <- clump(c_ea); k_hill <- clump(c_hill)
logl("clumped separately: EA4 {length(k_ea)}, Hill {length(k_hill)} (Li D Table S3 implies ~406 and ~25 before harmonization; check only)")
minp <- function(ids) { a <- ea[rsid %in% ids, .(rsid, p)]; h <- hill[rsid %in% ids, .(rsid, p)]
  rbind(a, h)[, .(p = min(p)), by = rsid] }

sets <- list()
u1 <- union(k_ea, k_hill)
sets$S1 <- clump(minp(u1)); sets$S2 <- u1
sets$S3 <- clump(minp(union(c_ea$rsid, c_hill$rsid)))
for (nm in names(sets)) addr(nm, "constructed (after clumping step)", length(sets[[nm]]))

# ---- build B matrices (both exposures present, alleles aligned) -------------------------------------------
buildB <- function(ids) {
  a <- align(ea, ids); h <- align(hill, ids)
  B <- merge(a[, .(rsid, b_EA = b, se_EA = se)], h[, .(rsid, b_Income = b, se_Income = se)], by = "rsid")
  merge(ref[, .(rsid, ref_ea, ref_oa, ref_eaf)], B, by = "rsid")
}
Bs <- lapply(sets, buildB)
for (nm in names(Bs)) addr(nm, "both exposures present, alleles aligned", nrow(Bs[[nm]]))

# ---- CAD availability (Aragam EUR) via TwoSampleMR harmonization flags only ------------------------------------
cad <- fread(here::here("data", "raw", "cad_aragam2022_eur.tsv"),
             select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value"))
setnames(cad, c("chr", "pos", "ea", "oa", "eaf", "beta", "se", "p"))
cad[, `:=`(chr = as.character(chr), ea = toupper(ea), oa = toupper(oa))]
pos_map <- unique(ea[, .(rsid, chr = as.character(chr), pos)], by = "rsid")
pos_map <- rbind(pos_map, unique(hill[!rsid %in% pos_map$rsid, .(rsid, chr = as.character(chr), pos)], by = "rsid"))
harmonise_keep <- function(B, nm) {
  m <- merge(B[, .(rsid, ref_ea, ref_oa, ref_eaf)], pos_map, by = "rsid")
  o <- merge(m, cad, by = c("chr", "pos"), allow.cartesian = TRUE)
  o <- o[(ea == ref_ea & oa == ref_oa) | (ea == ref_oa & oa == ref_ea) |
         (chartr("ACGT", "TGCA", ea) == ref_ea & chartr("ACGT", "TGCA", oa) == ref_oa) |
         (chartr("ACGT", "TGCA", ea) == ref_oa & chartr("ACGT", "TGCA", oa) == ref_ea)]
  o <- unique(o, by = "rsid")
  addr(nm, "present in Aragam EUR (chr:pos + compatible alleles; no proxies)", nrow(o))
  expo <- data.frame(SNP = o$rsid, effect_allele.exposure = o$ref_ea, other_allele.exposure = o$ref_oa,
                     eaf.exposure = o$ref_eaf, beta.exposure = 1, se.exposure = 1, exposure = "EA", id.exposure = "EA")
  outc <- data.frame(SNP = o$rsid, effect_allele.outcome = o$ea, other_allele.outcome = o$oa, eaf.outcome = o$eaf,
                     beta.outcome = o$beta, se.outcome = o$se, outcome = "CAD", id.outcome = "CAD")
  h <- suppressMessages(TwoSampleMR::harmonise_data(expo, outc, action = 2))
  keep <- h$SNP[h$mr_keep]
  addr(nm, "kept by harmonise_data(action = 2) (ambiguous palindromes / incompatible alleles removed)", length(keep))
  B[rsid %in% keep]
}
final <- list()
for (nm in c("S1", "S2", "S3")) final[[nm]] <- harmonise_keep(Bs[[nm]], nm)
is_pal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C")
final[["S1-ref"]] <- Bs$S1[!(is_pal(ref_ea, ref_oa) & pmin(ref_eaf, 1 - ref_eaf) > PAL_MAF)]
addr("S1-ref", "no CAD step; ambiguous palindromes (MAF>0.42) removed", nrow(final[["S1-ref"]]))

# ---- LD audit (all sets; decides whether S2 can enter the main comparison) --------------------------------------------
ld_audit <- rbindlist(lapply(names(final), function(nm) {
  ids <- final[[nm]]$rsid; f <- file.path(TMP, paste0(nm, ".txt")); writeLines(ids, f)
  out <- file.path(TMP, paste0(nm, "_ld"))
  st <- system2(PLINK_BIN, shQuote(c("--bfile", BFILE_EUR, "--extract", f, "--r2", "--ld-window-kb", CLUMP_KB, "--ld-window", 99999,
                       "--ld-window-r2", CLUMP_R2, "--out", out)), stdout = FALSE, stderr = FALSE)
  # fail loudly: a missing .ld file must never be read as "no LD pairs"
  if (st != 0 || !file.exists(paste0(out, ".ld"))) stop("plink --r2 failed for ", nm, " (exit ", st, ")")
  ld <- fread(paste0(out, ".ld"))
  ld <- ld[R2 >= CLUMP_R2]
  data.table(set = nm, n_snps = length(ids), pairs_r2_ge_0.001 = nrow(ld), max_r2 = if (nrow(ld)) max(ld$R2) else 0,
             snps_involved = if (nrow(ld)) uniqueN(c(ld$SNP_A, ld$SNP_B)) else 0L)
}))
write_result_table(ld_audit, "pilot_lid_cad_ld_audit")
s2_joint_indep <- ld_audit[set == "S2", pairs_r2_ge_0.001] == 0

# ---- fixed LD sensitivity (added 2026-10-06, rule set before running) ---------------------------------------
# From S1, remove SNPs until no pair has r2 >= 0.001 in 1000G EUR (plink --r2, 10 Mb): take pairs in descending r2;
# if both SNPs are still present, drop the one with the larger minimum p across the two exposures (the weaker
# instrument; ties: drop the lexicographically larger rsID). Then recompute conditional F on the pruned set.
ld_pairs <- function(ids, tag) {
  f <- file.path(TMP, paste0(tag, ".txt")); writeLines(ids, f); out <- file.path(TMP, paste0(tag, "_ld"))
  st <- system2(PLINK_BIN, shQuote(c("--bfile", BFILE_EUR, "--extract", f, "--r2", "--ld-window-kb", CLUMP_KB, "--ld-window", 99999,
                                     "--ld-window-r2", CLUMP_R2, "--out", out)), stdout = FALSE, stderr = FALSE)
  if (st != 0 || !file.exists(paste0(out, ".ld"))) stop("plink --r2 failed for ", tag, " (exit ", st, ")")
  ld <- fread(paste0(out, ".ld")); ld[R2 >= CLUMP_R2]
}
pr <- ld_pairs(final$S1$rsid, "S1_for_pruning")[order(-R2)]
mp <- minp(final$S1$rsid); setkey(mp, rsid)
keep <- final$S1$rsid; dropped <- character(0)
for (k in seq_len(nrow(pr))) {
  a <- pr$SNP_A[k]; b <- pr$SNP_B[k]
  if (!(a %in% keep && b %in% keep)) next
  pa <- mp[a, p]; pb <- mp[b, p]
  d <- if (pa > pb) a else if (pb > pa) b else max(a, b)
  keep <- setdiff(keep, d); dropped <- c(dropped, d)
}
final[["S1-LDpruned"]] <- final$S1[rsid %in% keep]
addr("S1-LDpruned", "S1 after removing one SNP per pair with r2 >= 0.001 (weaker instrument dropped)", length(keep))
chk <- ld_pairs(keep, "S1_pruned_check")
ld_audit <- rbind(ld_audit, data.table(set = "S1-LDpruned", n_snps = length(keep), pairs_r2_ge_0.001 = nrow(chk),
                                       max_r2 = if (nrow(chk)) max(chk$R2) else 0,
                                       snps_involved = if (nrow(chk)) uniqueN(c(chk$SNP_A, chk$SNP_B)) else 0L))
write_result_table(ld_audit, "pilot_lid_cad_ld_audit")
fwrite(data.table(rsid = dropped), file.path(DIR_IN, "S1-LDpruned_dropped.tsv"), sep = "\t")
logl("LD sensitivity: {nrow(pr)} pairs in S1; dropped {length(dropped)} SNPs; pruned set {length(keep)}; remaining pairs {nrow(chk)}")
logl("LD audit:\n{paste(capture.output(print(ld_audit)), collapse = '\n')}\nS2 meets joint independence: {s2_joint_indep}")

# ---- conditional F --------------------------------------------------------------------------------------
cf <- rbindlist(lapply(names(final), function(nm) {
  B <- final[[nm]]
  fwrite(B, file.path(DIR_IN, paste0(nm, "_B.tsv")), sep = "\t")
  bx <- as.matrix(B[, .(b_EA, b_Income)]); sx <- as.matrix(B[, .(se_EA, se_Income)])
  r_in <- MVMR::format_mvmr(BXGs = bx, BYG = rep(0, nrow(B)), seBXGs = sx, seBYG = rep(1, nrow(B)), RSID = B$rsid)
  L <- nrow(B)
  rbindlist(lapply(seq_along(RHO), function(i) {
    pc <- matrix(c(1, RHO[i], RHO[i], 1), 2)
    f <- suppressWarnings(if (RHO[i] == 0) MVMR::strength_mvmr(r_in, gencov = 0)
                          else MVMR::strength_mvmr(r_in, gencov = MVMR::phenocov_mvmr(pc, sx)))
    data.table(set = nm, role = if (nm == "S2" && !s2_joint_indep) "diagnostic" else if (nm == "S1-LDpruned") "LD sensitivity" else "main", n_snps = L,
               rho_scale = c(0, 0.5, 1, 1.5)[i], rho = RHO[i], exposure = c("EA", "Income"),
               condF = as.numeric(f[1, ]), condF_df = as.numeric(f[1, ]) * L / (L - 1))
  }))
}))
fwrite(cf, file.path(DIR_IN, "all_condF_unrounded.tsv"), sep = "\t")
cf_out <- copy(cf)[, `:=`(condF = round(condF, 2), condF_df = round(condF_df, 2), rho = round(rho, 4))]
write_result_table(cf_out, "pilot_lid_cad_condF")
ret <- rbindlist(retention); write_result_table(ret, "pilot_lid_cad_retention")
logl("retention:\n{paste(capture.output(print(ret)), collapse = '\n')}")
logl("conditional F:\n{paste(capture.output(print(cf_out)), collapse = '\n')}")
logl("Done.")
close_run_log(log_con)
