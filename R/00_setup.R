#' 00_setup.R
#'
#' Project-wide setup: paths, seed, package loading, and helper functions.
#' Source this at the top of every analysis script:
#'   source(here::here("R", "00_setup.R"))
#'
#' Inputs:  none
#' Outputs: side effects only (loaded packages, defined globals, defined helpers)
#' Depends on: here, fs (for path management); see required_packages below

# ---- 1. Paths ----------------------------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install with: install.packages('here')")
}

PROJ_ROOT       <- here::here()
DIR_DATA_RAW    <- here::here("data", "raw")
DIR_DATA_HARM   <- here::here("data", "harmonized")
DIR_INSTRUMENTS <- here::here("data", "instruments")
DIR_R           <- here::here("R")
DIR_R_UTILS     <- here::here("R", "utils")
DIR_RESULTS     <- here::here("results")
DIR_TABLES      <- here::here("results", "tables")
DIR_FIGURES     <- here::here("results", "figures")
DIR_LOGS        <- here::here("results", "logs")
DIR_MANUSCRIPT  <- here::here("manuscript")
DIR_NOTEBOOKS   <- here::here("notebooks")

PATH_MANIFEST   <- file.path(DIR_DATA_RAW, "MANIFEST.tsv")
PATH_SOURCES    <- file.path(DIR_DATA_RAW, "SOURCES.md")
PATH_DEVIATIONS <- file.path(DIR_LOGS, "DEVIATIONS.md")

# ---- 2. Seed -----------------------------------------------------------------

PROJECT_SEED <- 20260506L
set.seed(PROJECT_SEED)

# ---- 3. Required packages ----------------------------------------------------

required_packages <- c(
  # core
  "here", "fs", "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
  "data.table", "glue",
  # MR stack
  "TwoSampleMR", "MendelianRandomization", "MRPRESSO", "MVMR", "ieugwasr",
  # plotting
  "ggplot2", "patchwork", "scales",
  # io / reporting
  "openxlsx", "knitr", "gt", "digest"
)

ensure_pkgs_loaded <- function(pkgs = required_packages, quiet = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    warning(
      "Missing packages — install before running analysis: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(setdiff(pkgs, missing))
}

# ---- 4. Helpers --------------------------------------------------------------

#' Read a GWAS summary statistics file with sensible defaults.
#' Auto-detects gzipped files. Returns a data.table.
read_sumstats <- function(path, ...) {
  stopifnot(file.exists(path))
  data.table::fread(path, ...)
}

#' Open a timestamped log for a script and return the connection.
#' Close it with close_run_log(con) as the script's last line (on.exit() at top level never fires under Rscript).
log_run <- function(script_name) {
  fs::dir_create(DIR_LOGS)
  stamp <- format(Sys.time(), "%Y-%m-%d")
  log_path <- file.path(DIR_LOGS, glue::glue("{stamp}_{script_name}.log"))
  con <- file(log_path, open = "wt")
  writeLines(
    c(
      glue::glue("# Run log: {script_name}"),
      glue::glue("# Started: {Sys.time()}"),
      glue::glue("# Seed:    {PROJECT_SEED}"),
      glue::glue("# R:       {R.version.string}"),
      ""
    ),
    con
  )
  con
}

#' Capture sessionInfo() to a log connection. Call at the end of each script.
log_session_info <- function(con) {
  writeLines(c("", "# --- sessionInfo() ---"), con)
  capture.output(sessionInfo(), file = con, append = TRUE)
}

#' Finish a script log: append sessionInfo() and close the connection.
#' Call explicitly as the last line of a script. `on.exit()` must NOT be used
#' for this at script top level: under Rscript it never fires (fixed 2026-10-03).
close_run_log <- function(con) {
  ok <- tryCatch(isOpen(con), error = function(e) FALSE)
  if (ok) { log_session_info(con); close(con) }
  invisible(NULL)
}

#' Make a failing script still write sessionInfo() to its log, then exit non-zero.
register_log_on_error <- function(con) {
  options(error = function() {
    try(writeLines(c("", "# --- run FAILED ---", geterrmessage()), con), silent = TRUE)
    try(close_run_log(con), silent = TRUE)
    if (!interactive()) quit(save = "no", status = 1)
  })
  invisible(NULL)
}

#' Write a tabular result both as TSV (canonical) and XLSX (review-friendly).
write_result_table <- function(df, basename, dir = DIR_TABLES) {
  fs::dir_create(dir)
  tsv_path  <- file.path(dir, paste0(basename, ".tsv"))
  xlsx_path <- file.path(dir, paste0(basename, ".xlsx"))
  readr::write_tsv(df, tsv_path)
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(df, xlsx_path, overwrite = TRUE)
  }
  invisible(c(tsv = tsv_path, xlsx = xlsx_path))
}

#' Save a figure as PDF (canonical) + PNG (preview at 300 dpi).
save_figure <- function(plot, basename, width = 7, height = 5, dir = DIR_FIGURES) {
  fs::dir_create(dir)
  pdf_path <- file.path(dir, paste0(basename, ".pdf"))
  png_path <- file.path(dir, paste0(basename, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 300)
  invisible(c(pdf = pdf_path, png = png_path))
}

#' SHA-256 of a file — used for MANIFEST integrity checks.
file_sha256 <- function(path) {
  stopifnot(file.exists(path))
  digest::digest(file = path, algo = "sha256")
}

# ---- 5. Load packages quietly and announce -----------------------------------

ensure_pkgs_loaded()
suppressPackageStartupMessages({
  invisible(lapply(
    intersect(required_packages, rownames(installed.packages())),
    library, character.only = TRUE
  ))
})

message(glue::glue("[setup] PROJ_ROOT = {PROJ_ROOT}"))
message(glue::glue("[setup] seed     = {PROJECT_SEED}"))
