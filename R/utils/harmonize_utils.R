#' harmonize_utils.R
#'
#' Allele alignment (rules 2 and 4) and the retention summary, moved unchanged from R/03_harmonize.R on 2026-10-08 so
#' that R/03d (proxy sensitivity) uses the same verified code. Requires is_palindromic() from sumstats_registry.R.
HARM_MHC <- c(chr = "6", from = 25e6, to = 34e6)  # GRCh37, flag only (sensitivity)

#' x: rsid, ea_x, oa_x, eaf_x (+ any columns to carry); y: outcome extract
harmonize <- function(x, y, tol = 0.08, drop_all_pal = FALSE) {
  h <- merge(x, y, by = "rsid", all.x = TRUE)
  h[is.na(status), status := "missing"]
  h[, `:=`(palindromic = is_palindromic(ea_x, oa_x), drop_reason = fifelse(status == "ok", NA_character_, status))]
  ok <- h$status == "ok"
  same <- ok & h$ea_y == h$ea_x & h$oa_y == h$oa_x
  swap <- ok & h$ea_y == h$oa_x & h$oa_y == h$ea_x
  stopifnot(all(same | swap | !ok))  # finish() guarantees the pair; no strand-complement rescue (rule 2)
  h[, `:=`(beta_y_aligned = fifelse(swap, -beta_y, beta_y), eaf_y_aligned = fifelse(swap, 1 - eaf_y, eaf_y))]
  fA <- fifelse(is.na(h$eaf_x), 0.5, h$eaf_x); fB <- fifelse(is.na(h$eaf_y_aligned), 0.5, h$eaf_y_aligned)
  pal_flip <- ok & h$palindromic & ((fA < 0.5 & fB > 0.5) | (fA > 0.5 & fB < 0.5))
  h[pal_flip, `:=`(beta_y_aligned = -beta_y_aligned, eaf_y_aligned = 1 - eaf_y_aligned)]
  h[, pal_flipped := pal_flip]
  amb <- ok & h$palindromic & ((fA > 0.5 - tol & fA < 0.5 + tol) | (fB > 0.5 - tol & fB < 0.5 + tol))
  h[, ambiguous := amb]
  h[ok & amb, drop_reason := "palindrome_ambiguous"]
  if (drop_all_pal) h[ok & palindromic, drop_reason := "palindrome_dropped_afr"]
  h[, mr_keep := is.na(drop_reason)]
  h[, `:=`(mhc = chr == HARM_MHC[["chr"]] & pos >= as.numeric(HARM_MHC[["from"]]) & pos <= as.numeric(HARM_MHC[["to"]]))]
  h
}
retention <- function(h, set, outcome, type) {
  data.table(type = type, set = set, outcome = outcome, n_selected = nrow(h), n_at_key = sum(h$status != "missing"),
             n_allele_consistent = sum(h$status %in% c("ok", "duplicate")), n_duplicate_dropped = sum(h$status == "duplicate"),
             n_after_duplicates = sum(h$status == "ok"), n_palindromic = sum(h$status == "ok" & h$palindromic),
             n_pal_flipped = sum(h$pal_flipped & h$mr_keep), n_pal_dropped = sum(h$drop_reason %in% c("palindrome_ambiguous", "palindrome_dropped_afr")),
             n_final = sum(h$mr_keep), n_final_mhc = sum(h$mr_keep & h$mhc),
             n_final_mvp_r2_lt_0.8 = if ("r2" %in% names(h)) sum(h$mr_keep & !is.na(h$r2) & h$r2 < 0.8) else NA_integer_)
}
