# 01_pull_eurostat.R — Eurostat pulls (robust to TIME_PERIOD / values naming)
# Usage: source("code/01_pull_eurostat.R") from the repo root

# Optional: auto-detect project root if sourced via an absolute path
if (requireNamespace("rprojroot", quietly = TRUE)) {
  root <- tryCatch(
    rprojroot::find_root(rprojroot::has_file(".git") | rprojroot::has_file("README.md")),
    error = function(e) getwd()
  )
  setwd(root)
}

# Load environment and deps
source("code/00_setup.R")
library(eurostat)
library(data.table)
library(readr)
library(jsonlite)
library(stringr)

# Helper: save with simple metadata
save_with_meta <- function(df, name){
  f_csv <- file.path("data/raw/eurostat", paste0(name, ".csv"))
  f_rds <- file.path("data/raw/eurostat", paste0(name, ".rds"))
  readr::write_csv(df, f_csv)
  saveRDS(df, f_rds)
  meta <- list(dataset = name, extracted_at = as.character(Sys.time()), nrows = nrow(df))
  jsonlite::write_json(meta, paste0(f_csv, ".meta.json"), auto_unbox = TRUE, pretty = TRUE)
}

# Helper: normalize eurostat column names across versions
normalize_es <- function(dt) {
  nm <- names(dt)
  # GEO
  if ("GEO" %in% nm && !("geo" %in% nm)) data.table::setnames(dt, "GEO", "geo")
  # TIME
  if ("TIME_PERIOD" %in% nm && !("time" %in% nm)) data.table::setnames(dt, "TIME_PERIOD", "time")
  if ("TIME" %in% nm && !("time" %in% nm)) data.table::setnames(dt, "TIME", "time")
  # VALUES
  if ("values" %in% nm && !("value" %in% nm)) data.table::setnames(dt, "values", "value")
  if ("OBS_VALUE" %in% nm && !("value" %in% nm)) data.table::setnames(dt, "OBS_VALUE", "value")
  if ("obs_value" %in% nm && !("value" %in% nm)) data.table::setnames(dt, "obs_value", "value")
  if ("obsValue" %in% nm && !("value" %in% nm)) data.table::setnames(dt, "obsValue", "value")
  if ("VALUE" %in% nm && !("value" %in% nm)) data.table::setnames(dt, "VALUE", "value")
  
  # Create integer year from time if needed
  if (!("year" %in% names(dt))) {
    if ("time" %in% names(dt)) {
      if (is.numeric(dt$time)) {
        dt[, year := as.integer(time)]
      } else {
        dt[, year := as.integer(substr(as.character(time), 1, 4))]
      }
    } else {
      stop("No 'time' or 'TIME_PERIOD' (or 'year') column found in eurostat result.")
    }
  }
  # Drop 'time' if present (keep canonical 'year')
  if ("time" %in% names(dt)) dt[, time := NULL]
  dt
}

# -----------------------------
# 1) demo_r_mlifexp (LE by sex)
# -----------------------------
message("Downloading demo_r_mlifexp ...")
dem <- as.data.table(get_eurostat("demo_r_mlifexp", time_format = "num"))
dem <- normalize_es(dem)
# Keep NUTS-2 only (4-char code: e.g., NL41)
dem <- dem[nchar(geo) == 4 & grepl("^[A-Z]{2}..$", geo)]
# Years
dem <- dem[data.table::between(year, 1995, 2023)]
# Harmonize sex
if ("sex" %in% names(dem)) dem[, sex := factor(sex, levels = c("T","M","F"))]
save_with_meta(dem, "demo_r_mlifexp_NUTS2_1995_2023")
message(sprintf("Saved demo_r_mlifexp: %d rows", nrow(dem)))

# ------------------------------------
# 2) tgs00101 (LE at birth, NUTS-2)
# ------------------------------------
message("Downloading tgs00101 ...")
tgs <- as.data.table(get_eurostat("tgs00101", time_format = "num"))
tgs <- normalize_es(tgs)
tgs <- tgs[nchar(geo) == 4 & grepl("^[A-Z]{2}..$", geo)]
tgs <- tgs[data.table::between(year, 2000, 2023)]
if ("sex" %in% names(tgs)) tgs[, sex := factor(sex, levels = c("T","M","F"))]
save_with_meta(tgs, "tgs00101_NUTS2_2000_2023")
message(sprintf("Saved tgs00101: %d rows", nrow(tgs)))

message("Done: Eurostat extracts saved to data/raw/eurostat/.")
