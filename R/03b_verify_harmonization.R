#' 03b_verify_harmonization.R
#'
#' Purpose: verify every harmonized cell written by R/03 (univariable, MVMR, EA4 -> mediator) and R/03d (proxy cells), including the
#'   exceptional paths not covered by the in-script check: van Zuydam (palindrome tolerance 0.085) and MVP AFR (all
#'   palindromes dropped, TwoSampleMR action 3). For each cell, the allele-matched rows are re-harmonized with
#'   TwoSampleMR's internal harmonise_22() at the cell's tolerance and action, and TwoSampleMR's own mr_keep rules
#'   (action 2: drop remove | ambiguous; action 3: drop palindromic | remove | ambiguous). The retained SNP set and
#'   the aligned outcome beta must agree exactly; rows without an allele-matched outcome row must not be retained.
#'   Stops with an error on any mismatch. No MR estimate is computed.
#' Outputs: results/tables/harmonization_verify_tsmr.tsv; results/logs/{date}_03b_verify_harmonization.log
#' Depends on: data.table, TwoSampleMR

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(glue) })
log_con <- log_run("03b_verify_harmonization")
logl <- function(...) { m <- glue::glue(..., .envir = parent.frame()); writeLines(m, log_con); message(m); flush(log_con) }
register_log_on_error(log_con)
OUT <- file.path(DIR_DATA_HARM, "large")
dirs <- c("univariable", "mvmr", "mediation", "univariable_proxy", "mvmr_proxy")
files <- unlist(lapply(dirs, function(d) list.files(file.path(OUT, d), "\\.tsv\\.gz$", full.names = TRUE)))
stopifnot(length(files) == 108 + 516 + 7 + 108 + 516)
cell_rules <- function(o) list(tol = if (o == "pad_vanzuydam") 0.085 else 0.08, action = if (o == "pad_mvp_afr") 3 else 2)

res <- rbindlist(lapply(files, function(f) {
  h <- fread(f, colClasses = list(character = c("ea_x", "oa_x", "ea_y", "oa_y")))
  o <- sub("^.*__", "", sub("\\.tsv\\.gz$", "", basename(f))); r <- cell_rules(o)
  bad_keep <- sum(h$mr_keep & h$status != "ok")
  ok <- h[status == "ok"]
  mine <- setNames(h$beta_y_aligned[h$mr_keep], h$rsid[h$mr_keep])
  if (nrow(ok) == 0) return(data.table(dir = basename(dirname(f)), cell = basename(f), outcome = o, tol = r$tol, action = r$action,
                                       n_ok = 0L, n_mine = length(mine), n_tsmr = 0L, same_snps = length(mine) == 0, max_abs_diff = 0, bad_keep = bad_keep))
  d <- TwoSampleMR:::harmonise_22(ok$rsid, ok$ea_x, ok$oa_x, ok$ea_y, ok$oa_y, rep(0, nrow(ok)), ok$beta_y,
                                  ok$eaf_x, ok$eaf_y, r$tol, r$action)
  keep <- if (r$action == 3) !(d$palindromic | d$remove | d$ambiguous) else !(d$remove | d$ambiguous)
  tsmr <- setNames(d$beta.outcome[keep], d$SNP[keep])
  common <- intersect(names(mine), names(tsmr))
  data.table(dir = basename(dirname(f)), cell = basename(f), outcome = o, tol = r$tol, action = r$action, n_ok = nrow(ok),
             n_mine = length(mine), n_tsmr = length(tsmr), same_snps = setequal(names(mine), names(tsmr)),
             max_abs_diff = if (length(common)) max(abs(mine[common] - tsmr[common])) else 0, bad_keep = bad_keep)
}))
fwrite(res, file.path(DIR_TABLES, "harmonization_verify_tsmr.tsv"), sep = "\t")
logl("cells checked: {nrow(res)} ({paste(names(table(res$dir)), table(res$dir), collapse = ', ')})")
logl("of which van Zuydam {sum(res$outcome == 'pad_vanzuydam')} (tol 0.085), MVP AFR {sum(res$outcome == 'pad_mvp_afr')} (action 3)")
bad <- res[!same_snps | max_abs_diff > 1e-12 | bad_keep > 0]
if (nrow(bad) > 0) { print(bad); stop("harmonization verification failed in ", nrow(bad), " cell(s)") }
logl("All cells agree with TwoSampleMR harmonise_22 (retained SNPs identical; max |beta diff| {signif(max(res$max_abs_diff), 3)}).")
close_run_log(log_con)
