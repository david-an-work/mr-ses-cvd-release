#' 01_download_sumstats.R
#'
#' Auto-download the openly-hosted GWAS summary statistics for the MR-SES-CVD
#' project, and rewrite MANIFEST.tsv with sha256 / row-count / column /
#' p-value-sanity information for every file (manual + auto).
#'
#' Design:
#'   - Registry (`auto_downloads`) is the source of truth for URLs.
#'     Each row: dataset key, URL, target filename, expected role.
#'   - Manual file registry (`manual_files`) maps existing files in data/raw/
#'     to dataset keys.
#'   - `download_one()` skips already-present files (size + optional sha256
#'     match); otherwise curl'ed with resume + retry.
#'   - `verify_one()` computes sha256, row count, column names, and p-value
#'     column quantiles.
#'   - Results merged into MANIFEST.tsv; log to results/logs/.
#'
#' Inputs:  data/raw/MANIFEST.tsv (existing), data/raw/{manual files}
#' Outputs: data/raw/{dataset}.{ext} for auto-downloaded files;
#'          updated data/raw/MANIFEST.tsv;
#'          results/logs/{YYYY-MM-DD}_01_download_sumstats.log
#' Depends on: here, digest, glue, dplyr, tibble, readr, data.table

suppressPackageStartupMessages({
  library(here); library(digest); library(glue)
  library(dplyr); library(tibble); library(readr); library(data.table)
  library(R.utils) # data.table::fread() calls this transparently for .gz inputs
})

PROJ_ROOT    <- here::here()
DIR_DATA_RAW <- here::here("data", "raw")
DIR_LOGS     <- here::here("results", "logs")
PATH_MANIFEST <- file.path(DIR_DATA_RAW, "MANIFEST.tsv")
if (!dir.exists(DIR_LOGS)) dir.create(DIR_LOGS, recursive = TRUE)

set.seed(20260506)

# ---- Log ---------------------------------------------------------------------

stamp <- format(Sys.time(), "%Y-%m-%d")
log_path <- file.path(DIR_LOGS, glue::glue("{stamp}_01_download_sumstats.log"))
log_con  <- file(log_path, open = "wt")
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), "  ",
                glue::glue(..., .envir = parent.frame()))
  writeLines(msg, log_con); message(msg); flush(log_con)
}
# sessionInfo() is appended by finish_log() at the end of the script and on error.
# (A top-level on.exit() never fires under Rscript; fixed 2026-10-03.)
finish_log <- function() {
  if (tryCatch(isOpen(log_con), error = function(e) FALSE)) {
    writeLines(c("", "# --- sessionInfo() ---"), log_con)
    capture.output(sessionInfo(), file = log_con, append = TRUE)
    close(log_con)
  }
}
options(error = function() {
  try(writeLines(c("", "# --- run FAILED ---", geterrmessage()), log_con), silent = TRUE)
  try(finish_log(), silent = TRUE)
  if (!interactive()) quit(save = "no", status = 1)
})

log_line("Log start. seed=20260506")
log_line("PROJ_ROOT = {PROJ_ROOT}")

# ---- Auto-download registry --------------------------------------------------

# Primary + minimal-sensitivity batch. Supplementary strata (Aragam EUR,
# GIGASTROKE ischemic + EUR, HERMES subtypes, Akimova SIOPS/ISEI) can be added
# to the registry when Phase 5-7 need them.
auto_downloads <- tribble(
  ~dataset,                       ~url,                                                                                                                                                          ~target,                                          ~role,
  # --- Outcomes ---
  "cad_aragam2022_multi",         "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132315/GCST90132315_buildGRCh37.tsv",                             "cad_aragam2022_multi.tsv",                       "outcome primary",
  "hf_hermes2_2025",              "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90726001-GCST90727000/GCST90726617/GCST90726617.tsv.gz",                                     "hf_hermes2_2025_hfall.tsv.gz",                   "outcome primary",
  "stroke_gigastroke_any_multi",  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104534/GCST90104534_buildGRCh37.tsv.gz",                         "stroke_gigastroke_any_multi.tsv.gz",             "outcome primary",
  "pad_vanzuydam2021",            "https://ndownloader.figshare.com/files/24085532",                                                                                                             "pad_vanzuydam2021.tar.gz",                       "outcome primary",
  # --- Exposures ---
  "occupation_akimova2024",       "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446160/GCST90446160.tsv",                                        "occupation_akimova2024_camsis.tsv",              "exposure primary",
  "tdi_panukb",                   "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-189-both_sexes-irnt.tsv.bgz",                                                     "tdi_panukb_189_irnt.tsv.bgz",                    "exposure primary",
  # Kweon 2025 Income Factor (added 2026-10-03): DeSci Nodes dPID 149 (doi:10.62891/aac85602), file DATA/INC_FACTOR_GWAS.txt via its IPFS CID
  "income_kweon2025",             "https://ipfs.desci.com/ipfs/bafybeig2rtgpuxleecxixzwpi7qvnr7aqevgydrd7gobuafwytc7wo5ewq",                                                                   "income_kweon2025_incfactor.txt",                 "exposure primary",
  # --- Mediators ---
  "bmi_giant_yengo2018",          "https://zenodo.org/api/records/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz/content",                                              "bmi_giant_yengo2018_combined.txt.gz",            "mediator",
  "ldl_glgc_graham2021",          "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1.gz",         "ldl_glgc_graham2021_trans_ancestry.tsv.gz",      "mediator",
  "hdl_glgc_graham2021",          "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz",         "hdl_glgc_graham2021_trans_ancestry.tsv.gz",      "mediator",
  "tg_glgc_graham2021",           "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz",       "tg_glgc_graham2021_trans_ancestry.tsv.gz",       "mediator",
  "sbp_keaton2024",               "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90310001-GCST90311000/GCST90310294/GCST90310294.tsv.gz",                                    "sbp_keaton2024.tsv.gz",                          "mediator primary",
  # --- PAD sensitivity (FinnGen R12) ---
  "pad_finngen_r12_athscle",      "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_ATHSCLE.gz",                                                    "pad_finngen_r12_athscle.tsv.gz",                 "outcome sensitivity",
  "pad_finngen_r12_othper",       "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_OTHPER.gz",                                                     "pad_finngen_r12_othper.tsv.gz",                  "outcome sensitivity",
  "pad_finngen_r12_artembthrlow", "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_ARTEMBTHRLOW.gz",                                              "pad_finngen_r12_artembthrlow.tsv.gz",            "outcome sensitivity",
  # --- Added 2026-10-06/07 after user approval of plan v0.6 scope ---
  "stroke_gigastroke_any_eur",    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104539/GCST90104539_buildGRCh37.tsv.gz", "stroke_gigastroke_any_eur.tsv.gz", "outcome sensitivity (EUR LDSC)",
  "occupation_akimova2024_isei",  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446163/GCST90446163.tsv", "occupation_akimova2024_isei.tsv", "exposure (Le-comparable)",
  # MVP gwPheWAS EUR (GIA) via the GWAS Catalog deposit; identical to the dbGaP files (Phe_411 checked 2026-10-07)
  "hf_mvp_phe428_eur",            "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475981/GCST90475981.tsv.gz", "hf_mvp_phe428_EUR.gwascat.tsv.gz", "outcome comparator (UKB-free)",
  "stroke_mvp_phe433_21_eur",     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90477001-GCST90478000/GCST90477991/GCST90477991.tsv.gz", "stroke_mvp_phe433_21_EUR.gwascat.tsv.gz", "outcome comparator (UKB-free)",
  # Latent SES candidate (feasibility check 2026-10-07; not adopted)
  "ses_gses_xia2025",             "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90566001-GCST90567000/GCST90566708/GCST90566708.tsv.gz", "ses_gses_xia2025.tsv.gz", "exposure (latent SES candidate)"
)

log_line("Registry has {nrow(auto_downloads)} auto-download rows.")

# ---- Byte-range registry (members of uncompressed tar archives) --------------

# MVP gwPheWAS (dbGaP phs002453, open access; added 2026-10-03). Phe_440_2 files
# are members of a 149.7 GB uncompressed tar that is not expanded on the FTP.
# Offsets were found by walking the tar headers (2026-10-03); each range is a
# complete .txt.gz. Offsets are only valid for this exact tar; re-derive them if
# dbGaP re-releases it.
MVP_TAR_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/analyses/GIA/",
  "phs002453.MVP_R4.1000G_AGR.GIA.PheCodes_CirculatorySystem_batch2.analysis-PI.MULTI.tar"
)
range_downloads <- tribble(
  ~dataset,                 ~start,         ~size,        ~target,                           ~member,
  "pad_mvp_phe440_2_eur",   5791620096,     592649204,    "pad_mvp_phe440_2_EUR.txt.gz",     "MVP_R4.1000G_AGR.Phe_440_2.EUR.GIA.dbGaP.txt.gz",
  "pad_mvp_phe440_2_afr",   15712993280,    1145560336,   "pad_mvp_phe440_2_AFR.txt.gz",     "MVP_R4.1000G_AGR.Phe_440_2.AFR.GIA.dbGaP.txt.gz",
  "pad_mvp_phe440_2_meta",  70992927744,    875714288,    "pad_mvp_phe440_2_META.txt.gz",    "MVP_R4.1000G_AGR.Phe_440_2.META.GIA.dbGaP.txt.gz",
  "pad_mvp_phe440_2_amr",   127403379200,   866420962,    "pad_mvp_phe440_2_AMR.txt.gz",     "MVP_R4.1000G_AGR.Phe_440_2.AMR.GIA.dbGaP.txt.gz"
) %>% mutate(url = MVP_TAR_URL)
# Phe_411 (ischaemic heart disease; UKB-free CAD comparator), added 2026-10-06 after user approval. Member of the
# circulatory batch1 tar; offset found by walking the tar headers (GNU ././@LongLink names), 2026-10-06.
MVP_TAR_URL_B1 <- sub("batch2", "batch1", MVP_TAR_URL)
range_downloads <- bind_rows(range_downloads, tribble(
  ~dataset,                 ~start,         ~size,        ~target,                           ~member,
  "cad_mvp_phe411_eur",     12125504512,    580351728,    "cad_mvp_phe411_EUR.txt.gz",       "MVP_R4.1000G_AGR.Phe_411.EUR.GIA.dbGaP.txt.gz"
) %>% mutate(url = MVP_TAR_URL_B1))

# ---- Manual file registry ----------------------------------------------------

# The user already downloaded these via source-portal click-through workflows.
# Map on-disk filenames to dataset keys.
manual_files <- tribble(
  ~existing_path,                                                              ~dataset,
  "Okbay 23andMe Additive Exclusion.txt.gz",                                   "ea_okbay2022",
  "Lee GWAS EA Excl 23andMe.txt",                                              "ea_lee2018",
  "Hill Data/Household_Income_UKBiobank.txt.gz",                               "income_hill2019",
  "Wootton Data/2019.10.02 Lifetime Smoking GWAS Data Sheet 1.txt",            "smoking_wootton2020",
  "DIAMANTE-TA.sumstat.txt",                                                   "t2d_diamante_mahajan2022"
)

# ---- Helpers -----------------------------------------------------------------

curl_bin <- Sys.which("curl")
if (curl_bin == "") stop("curl not found on PATH")

#' HEAD-check a URL; return content-length or NA.
head_size <- function(url) {
  res <- suppressWarnings(system2(
    curl_bin,
    c("-sIL", "--max-time", "15", url),
    stdout = TRUE, stderr = FALSE
  ))
  cl_line <- grep("^content-length:", tolower(res), value = TRUE)
  if (length(cl_line) == 0L) return(NA_integer_)
  as.numeric(sub("^content-length:\\s*", "", tail(cl_line, 1), ignore.case = TRUE))
}

#' Download one URL to a target path with resume + retry.
#' Returns TRUE on success, FALSE on any failure.
download_one <- function(url, target_path, expected_size = NA, retries = 3) {
  if (file.exists(target_path)) {
    local_size <- file.info(target_path)$size
    if (!is.na(expected_size) && expected_size > 0L &&
        !is.na(local_size) && abs(local_size - expected_size) < 1024L) {
      log_line("[skip] {basename(target_path)} already present ({local_size} bytes)")
      return(TRUE)
    }
    # No size to compare against (NA or figshare-style 0). Trust that any
    # non-empty file on disk was written by a previous successful run.
    if ((is.na(expected_size) || expected_size == 0L) &&
        !is.na(local_size) && local_size > 0L) {
      log_line("[skip] {basename(target_path)} already present ({local_size} bytes; no size check)")
      return(TRUE)
    }
  }
  dir.create(dirname(target_path), showWarnings = FALSE, recursive = TRUE)
  for (attempt in seq_len(retries)) {
    log_line("[dl {attempt}/{retries}] {basename(target_path)} <- {url}")
    rc <- suppressWarnings(system2(
      curl_bin,
      c("-sSL", "-C", "-", "--retry", "5", "--retry-delay", "10",
        "--connect-timeout", "30", "--max-time", "3600",
        "-o", shQuote(target_path), url),
      stdout = "", stderr = ""
    ))
    if (rc == 0L && file.exists(target_path) && file.info(target_path)$size > 0L) {
      local_size <- file.info(target_path)$size
      # Skip size check when expected_size is missing OR zero (some hosts,
      # notably figshare, redirect through a signed URL and the initial HEAD
      # returns content-length: 0 for the redirect body).
      if (!is.na(expected_size) && expected_size > 0L &&
          abs(local_size - expected_size) > 1024L) {
        log_line("[warn] {basename(target_path)} size {local_size} vs expected {expected_size}; retrying")
        file.remove(target_path)
        next
      }
      log_line("[ok] {basename(target_path)} {local_size} bytes")
      return(TRUE)
    }
  }
  log_line("[fail] {basename(target_path)} after {retries} attempts")
  FALSE
}

#' Download bytes [start, start + size) of `url` to `target_path`, resuming by
#' appending the remaining range after dropped connections. TRUE when the file
#' reaches exactly `size` bytes.
download_range <- function(url, start, size, target_path, max_tries = 60L) {
  have <- if (file.exists(target_path)) file.info(target_path)$size else 0
  if (have == size) {
    log_line("[skip] {basename(target_path)} already present ({have} bytes)")
    return(TRUE)
  }
  if (have > size) { log_line("[warn] {basename(target_path)} larger than expected; restarting"); file.remove(target_path); have <- 0 }
  for (attempt in seq_len(max_tries)) {
    rng <- sprintf("%.0f-%.0f", start + have, start + size - 1)
    log_line("[range {attempt}] {basename(target_path)} bytes {rng}")
    tmp <- paste0(target_path, ".part")
    suppressWarnings(system2(curl_bin,
      c("-sS", "--connect-timeout", "30", "--max-time", "1800", "-r", rng, "-o", shQuote(tmp), url),
      stdout = "", stderr = ""))
    if (file.exists(tmp)) {
      file.append(target_path, tmp); file.remove(tmp)
    }
    have <- if (file.exists(target_path)) file.info(target_path)$size else 0
    if (have == size) { log_line("[ok] {basename(target_path)} {have} bytes"); return(TRUE) }
  }
  log_line("[fail] {basename(target_path)} at {have}/{size} bytes after {max_tries} attempts")
  FALSE
}

#' Compute file sha256 (base R + digest).
file_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

#' Count non-header rows in a text file (auto-detects .gz / .bgz).
#' Shells out to `wc -l` (plain) or `gunzip -c | wc -l` (compressed) since
#' streaming a 20M-row file in R just to count lines is ~30x slower.
count_rows <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  is_gz <- grepl("\\.gz$|\\.bgz$", path, ignore.case = TRUE)
  cmd <- if (is_gz) {
    glue::glue("gunzip -c {shQuote(path)} | wc -l")
  } else {
    glue::glue("wc -l < {shQuote(path)}")
  }
  raw <- tryCatch(
    system2("bash", c("-c", shQuote(cmd)), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  )
  if (length(raw) == 0L) return(NA_integer_)
  # wc output looks like "     12345"; take the first integer token
  first_num <- suppressWarnings(as.integer(
    strsplit(trimws(raw[1]), "\\s+")[[1]][1]
  ))
  if (is.na(first_num)) return(NA_integer_)
  first_num - 1L  # subtract header row
}

#' Peek at first N rows to get column names (auto-detects .gz).
peek_columns <- function(path, n = 5L) {
  if (!file.exists(path)) return(NA_character_)
  tryCatch({
    x <- data.table::fread(path, nrows = n, showProgress = FALSE)
    paste(colnames(x), collapse = "|")
  }, error = function(e) glue::glue("ERROR: {conditionMessage(e)}"))
}

#' Sanity-check the p-value distribution. Returns quantiles of the p-value
#' column as a compact string.
pvalue_sanity <- function(path, pvalue_candidates = c(
  "p", "pval", "p_value", "P", "PVAL", "P_value", "P.value", "pvalue",
  "Pval", "PVAL", "P-value", "p.value", "P_VALUE",
  "MR-MEGA_p-value_association" # DIAMANTE trans-ancestry uses this
)) {
  if (!file.exists(path)) return(NA_character_)
  tryCatch({
    x <- data.table::fread(path, nrows = 200000L, showProgress = FALSE)
    pcol <- intersect(pvalue_candidates, colnames(x))
    if (length(pcol) == 0L) {
      # try fuzzy match
      lower <- tolower(colnames(x))
      pcol <- colnames(x)[grepl("^p(_?val(ue)?)?$", lower)]
    }
    # Pan-UKB stores -log10(p) (e.g. neglog10_pval_meta) instead of p.
    neglog10 <- FALSE
    if (length(pcol) == 0L) {
      pcol <- grep("^neglog10_pval(_meta)?$", colnames(x), value = TRUE)
      neglog10 <- length(pcol) > 0L
    }
    if (length(pcol) == 0L) return("no-pvalue-col")
    pvals <- suppressWarnings(as.numeric(x[[pcol[1]]]))
    if (neglog10) pvals <- 10^(-pvals)
    pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
    if (length(pvals) < 100L) return(glue::glue("only {length(pvals)} valid pvals in first 200k rows"))
    q <- quantile(pvals, probs = c(0, 0.001, 0.01, 0.5, 0.99), na.rm = TRUE)
    glue::glue("col={pcol[1]} n={length(pvals)} min={signif(q[1],2)} q0.001={signif(q[2],2)} q0.01={signif(q[3],2)} median={signif(q[4],3)} q0.99={signif(q[5],3)}")
  }, error = function(e) glue::glue("ERROR: {conditionMessage(e)}"))
}

# ---- Phase A: auto-downloads -------------------------------------------------

log_line("--- Phase A: auto-downloads ---")
auto_downloads <- auto_downloads %>%
  mutate(target_path = file.path(DIR_DATA_RAW, target))

# Sequential download (parallel risks hammering servers; ~19 GB primary batch).
auto_downloads$success <- FALSE
for (i in seq_len(nrow(auto_downloads))) {
  row <- auto_downloads[i, ]
  expected <- head_size(row$url)
  auto_downloads$success[i] <- download_one(row$url, row$target_path, expected)
}

log_line("Phase A complete. Successes: {sum(auto_downloads$success)}/{nrow(auto_downloads)}")

range_downloads <- range_downloads %>% mutate(target_path = file.path(DIR_DATA_RAW, target))
range_downloads$success <- FALSE
for (i in seq_len(nrow(range_downloads))) {
  row <- range_downloads[i, ]
  range_downloads$success[i] <- download_range(row$url, row$start, row$size, row$target_path)
}
log_line("Phase A (byte ranges) complete. Successes: {sum(range_downloads$success)}/{nrow(range_downloads)}")

# Special handling: van Zuydam PAD is a .tar.gz containing the sumstats file.
# Extract if not already extracted.
vz_tar <- file.path(DIR_DATA_RAW, "pad_vanzuydam2021.tar.gz")
vz_dir <- file.path(DIR_DATA_RAW, "pad_vanzuydam2021")
if (file.exists(vz_tar) && !dir.exists(vz_dir)) {
  log_line("Extracting van Zuydam tar.gz")
  dir.create(vz_dir, showWarnings = FALSE)
  system2("tar", c("-xzf", shQuote(vz_tar), "-C", shQuote(vz_dir)))
  extracted <- list.files(vz_dir, recursive = TRUE)
  log_line("van Zuydam contents: {paste(extracted, collapse=', ')}")
}

# ---- Phase B: build unified file inventory -----------------------------------

log_line("--- Phase B: build unified file inventory ---")

# Combine auto + manual into one inventory keyed by dataset.
inventory <- bind_rows(
  auto_downloads %>%
    transmute(dataset, path = target_path, source = "auto"),
  range_downloads %>%
    transmute(dataset, path = target_path, source = "auto (byte range)"),
  manual_files %>%
    mutate(path = file.path(DIR_DATA_RAW, existing_path), source = "manual") %>%
    select(dataset, path, source)
) %>%
  # van Zuydam ships as a tarball; verify the extracted primary sumstats, not the
  # archive. Before 2026-10-03 the checks ran on the .tar.gz and logged a newline-byte
  # count (7,349,496) and fread errors instead of real values (see DEVIATIONS.md).
  mutate(path = ifelse(
    dataset == "pad_vanzuydam2021",
    file.path(vz_dir, "final_meta", "pad_primary.txt.gz"),
    path
  )) %>%
  mutate(
    exists     = file.exists(path),
    size_bytes = ifelse(exists, file.info(path)$size, NA_real_)
  )

log_line("Inventory: {sum(inventory$exists)}/{nrow(inventory)} files present on disk.")
if (file.exists(vz_tar)) {
  log_line("van Zuydam tarball (download provenance) sha256={file_sha256(vz_tar)}")
}

# ---- Phase C: integrity verification ----------------------------------------

log_line("--- Phase C: integrity verification ---")

inventory <- inventory %>%
  mutate(
    sha256      = NA_character_,
    n_rows      = NA_integer_,
    columns     = NA_character_,
    pval_check  = NA_character_
  )

for (i in seq_len(nrow(inventory))) {
  if (!inventory$exists[i]) next
  row <- inventory[i, ]
  log_line("[verify] {row$dataset}  ({basename(row$path)})")
  inventory$sha256[i]     <- file_sha256(row$path)
  inventory$n_rows[i]     <- count_rows(row$path)
  inventory$columns[i]    <- peek_columns(row$path)
  inventory$pval_check[i] <- pvalue_sanity(row$path)
  log_line("  sha256={substr(inventory$sha256[i],1,16)}... n_rows={inventory$n_rows[i]}")
  log_line("  cols={substr(inventory$columns[i],1,120)}")
  log_line("  pval={inventory$pval_check[i]}")
}

# ---- Phase D: merge into MANIFEST.tsv ---------------------------------------

log_line("--- Phase D: merge integrity results into MANIFEST.tsv ---")

manifest <- data.table::fread(
  PATH_MANIFEST, sep = "\t", header = TRUE,
  colClasses = "character", fill = TRUE
)
today <- format(Sys.Date(), "%Y-%m-%d")

for (i in seq_len(nrow(inventory))) {
  if (!inventory$exists[i]) next
  ds <- inventory$dataset[i]
  match_row <- which(manifest$dataset == ds)
  if (length(match_row) == 0L) {
    log_line("[warn] dataset '{ds}' not in MANIFEST -- skipping merge")
    next
  }
  manifest[match_row, `:=`(
    download_date = today,
    sha256        = inventory$sha256[i],
    n_snps        = as.character(inventory$n_rows[i]),
    notes         = paste0(
      manifest[match_row, notes],
      " | verified ", today,
      " | cols=", substr(inventory$columns[i], 1, 200),
      " | pval=", inventory$pval_check[i]
    )
  )]
}

data.table::fwrite(manifest, PATH_MANIFEST, sep = "\t", quote = FALSE, na = "")
log_line("MANIFEST.tsv rewritten with sha256 / n_rows / cols / pval-check for {sum(inventory$exists)} files")

# Also write a compact integrity summary to results/logs/
inventory_out <- inventory %>%
  select(dataset, source, exists, size_bytes, n_rows, sha256, columns, pval_check)
data.table::fwrite(
  inventory_out,
  file.path(DIR_LOGS, glue::glue("{stamp}_01_integrity_summary.tsv")),
  sep = "\t", quote = FALSE, na = ""
)
log_line("Integrity summary written to results/logs/{stamp}_01_integrity_summary.tsv")

log_line("--- 01_download_sumstats.R DONE ---")
finish_log()
