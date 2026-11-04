#' 02j_gses_feasibility.R
#'
#' Purpose: bounded feasibility check of an existing latent SES factor (Xia 2025 gSES, GCST90566708) as a candidate
#'   secondary analysis (not adopted). Exposure-side only; no outcome data are used.
#'   (1) genetic correlations of gSES with our exposures (LDSC via GenomicSEM; our munged files are read from the
#'       02b cache, unchanged; gSES is munged into its own folder);
#'   (2) Q_SNP heterogeneity (`q_value`) at the selected gSES instruments (data/instruments/ses_gses_xia2025_instruments.tsv,
#'       selected by R/02 with SELECT_ONLY=ses_gses_xia2025);
#'   (3) instrument overlap with EA4 and Kweon (shared rsIDs).
#' Outputs: results/tables/gses_feasibility_{rg,instruments}.{tsv,xlsx}; results/logs/{date}_02j_gses_feasibility.log
#' Depends on: data.table, GenomicSEM

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
ensure_pkgs_loaded(c("data.table", "GenomicSEM"))
log_con <- log_run("02j_gses_feasibility")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
source(here::here("R", "utils", "sumstats_registry.R"))

OUT_LDSC <- here::here("data", "harmonized", "large", "ldsc")
OUT_G <- file.path(OUT_LDSC, "gses"); fs::dir_create(OUT_G)
LD_DIR <- paste0(here::here("data", "reference", "eur_w_ld_chr"), "/")   # as in R/02b
HM3    <- file.path(LD_DIR, "w_hm3.snplist")
logl("HM3 {HM3}; LD {LD_DIR}")

# munge gSES (own folder)
g <- fread(registry$ses_gses_xia2025$file, select = c("rs_id", "effect_allele", "other_allele", "beta", "p_value", "n", "effect_allele_frequency", "q_value"))
hm3 <- fread(HM3)
if (!file.exists(file.path(OUT_G, "gSES.sumstats.gz"))) {
g_in <- g[rs_id %in% hm3$SNP, .(SNP = rs_id, A1 = effect_allele, A2 = other_allele, BETA = beta, P = p_value, N = n,
                                MAF = pmin(effect_allele_frequency, 1 - effect_allele_frequency))]
f_in <- file.path(OUT_G, "gSES_input.tsv"); fwrite(g_in, f_in, sep = "\t")
owd <- setwd(OUT_G)
GenomicSEM::munge(files = f_in, hm3 = HM3, trait.names = "gSES", N = NA, info.filter = 0.9, maf.filter = 0.01, log.name = "munge_gses")
setwd(owd); file.remove(f_in)
}

traits <- c(EA = file.path(OUT_LDSC, "EA.sumstats.gz"), Income = file.path(OUT_LDSC, "Income.sumstats.gz"),
            TDI = file.path(OUT_LDSC, "TDI.sumstats.gz"), Occupation = file.path(OUT_LDSC, "Occupation.sumstats.gz"),
            Income_Hill = file.path(OUT_LDSC, "Income_Hill.sumstats.gz"), gSES = file.path(OUT_G, "gSES.sumstats.gz"))
stopifnot(all(file.exists(traits)))
owd <- setwd(OUT_G)
ld <- GenomicSEM::ldsc(traits = unname(traits), sample.prev = rep(NA, length(traits)), population.prev = rep(NA, length(traits)),
                       ld = LD_DIR, wld = LD_DIR, trait.names = names(traits), stand = TRUE, ldsc.log = "ldsc_gses")
setwd(owd)
k <- length(traits); nm <- names(traits)
rg <- ld$S_Stand; dimnames(rg) <- list(nm, nm)
idx <- which(lower.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE); se_rg <- matrix(NA_real_, k, k, dimnames = list(nm, nm))
se_v <- sqrt(diag(ld$V_Stand)); for (r in seq_len(nrow(idx))) { i <- idx[r, 1]; j <- idx[r, 2]; se_rg[i, j] <- se_rg[j, i] <- se_v[r] }
I <- ld$I; dimnames(I) <- list(nm, nm)
rg_tab <- data.table(trait = setdiff(nm, "gSES"), rg_with_gSES = round(rg["gSES", setdiff(nm, "gSES")], 3),
                     se = round(se_rg["gSES", setdiff(nm, "gSES")], 3), cross_trait_intercept = round(I["gSES", setdiff(nm, "gSES")], 3))
rg_tab[, h2_gSES := round(diag(ld$S)[nm == "gSES"], 4)]
write_result_table(rg_tab, "gses_feasibility_rg")
logl("rg with gSES:\n{paste(capture.output(print(rg_tab)), collapse = '\n')}")

# instruments: Q_SNP and overlap
ins <- fread(file.path(DIR_INSTRUMENTS, "ses_gses_xia2025_instruments.tsv"))
q <- g[rs_id %in% ins$rsid, .(rsid = rs_id, q_value)]
ins <- merge(ins, q, by = "rsid", all.x = TRUE)
ea <- fread(file.path(DIR_INSTRUMENTS, "ea_okbay2022_instruments.tsv"))$rsid
kw <- fread(file.path(DIR_INSTRUMENTS, "income_kweon2025_instruments.tsv"))$rsid
summ <- data.table(n_instruments = nrow(ins), n_q_lt_5e8 = sum(ins$q_value < 5e-8, na.rm = TRUE),
                   n_q_lt_0.05_bonf = sum(ins$q_value < 0.05 / nrow(ins), na.rm = TRUE),
                   shared_with_EA4 = sum(ins$rsid %in% ea), shared_with_Kweon = sum(ins$rsid %in% kw))
write_result_table(summ, "gses_feasibility_instruments")
logl("instruments:\n{paste(capture.output(print(summ)), collapse = '\n')}")
logl("Done.")
close_run_log(log_con)
