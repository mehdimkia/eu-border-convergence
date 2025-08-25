# 01_pull_eurostat.R — Eurostat pulls + NUTS-2016 harmonization
# Usage: source("code/01_pull_eurostat.R") from the repo root

# -------- Project root (deterministic) --------
if (requireNamespace("rprojroot", quietly = TRUE)) {
  root <- tryCatch(
    rprojroot::find_root(rprojroot::has_file("code/00_setup.R")),
    error = function(e) getwd()
  )
  setwd(root)
}

# -------- Env & packages --------
source("code/00_setup.R")  # expects dirs + EU27 etc.; safe if already sourced
suppressPackageStartupMessages({
  library(eurostat)
  library(data.table)
  library(readr)
  library(jsonlite)
  library(stringr)
})

# Ensure raw/derived dirs exist (idempotent)
dir.create("data/raw/eurostat", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw/crosswalks", recursive = TRUE, showWarnings = FALSE)
dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)

# -------- Helpers --------
save_with_meta <- function(df, name){
  f_csv <- file.path("data/raw/eurostat", paste0(name, ".csv"))
  f_rds <- file.path("data/raw/eurostat", paste0(name, ".rds"))
  readr::write_csv(df, f_csv)
  saveRDS(df, f_rds)
  meta <- list(dataset = name, extracted_at = as.character(Sys.time()), nrows = nrow(df))
  jsonlite::write_json(meta, paste0(f_csv, ".meta.json"), auto_unbox = TRUE, pretty = TRUE)
}

normalize_es <- function(dt) {
  dt <- as.data.table(dt)
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

# =============================
# Raw pulls to data/raw/eurostat
# =============================

# 1) Life expectancy (by sex), NUTS-2, 1995–2023
message("Downloading demo_r_mlifexp ...")
dem <- as.data.table(eurostat::get_eurostat("demo_r_mlifexp", time_format = "num"))
dem <- normalize_es(dem)
# Keep NUTS-2 only (4-char code like NL41)
dem <- dem[nchar(geo) == 4 & grepl("^[A-Z]{2}..$", geo)]
# Years
dem <- dem[data.table::between(year, 1995, 2023)]
# Harmonize sex levels
if ("sex" %in% names(dem)) dem[, sex := factor(sex, levels = c("T","M","F"))]
save_with_meta(dem, "demo_r_mlifexp_NUTS2_1995_2023")
message(sprintf("Saved demo_r_mlifexp: %d rows", nrow(dem)))

# 2) LE at birth, NUTS-2, 2000–2023 (auxiliary, sometimes cleaner for LE@birth)
message("Downloading tgs00101 ...")
tgs <- as.data.table(eurostat::get_eurostat("tgs00101", time_format = "num"))
tgs <- normalize_es(tgs)
tgs <- tgs[nchar(geo) == 4 & grepl("^[A-Z]{2}..$", geo)]
tgs <- tgs[data.table::between(year, 2000, 2023)]
if ("sex" %in% names(tgs)) tgs[, sex := factor(sex, levels = c("T","M","F"))]
save_with_meta(tgs, "tgs00101_NUTS2_2000_2023")
message(sprintf("Saved tgs00101: %d rows", nrow(tgs)))

# 3) Population (total), NUTS-2, 1995–2023  ← NEW
message("Downloading demo_r_pjangrp3 ...")
pop <- as.data.table(eurostat::get_eurostat("demo_r_pjangrp3", time_format = "num"))
pop <- normalize_es(pop)
# Keep NUTS-2 only, sex=T, age=TOTAL, year range
pop <- pop[
  nchar(geo) == 4 & grepl("^[A-Z]{2}..$", geo) &
    sex == "T" & age == "TOTAL" & data.table::between(year, 1995, 2023)
][, .(pop = as.numeric(value)), by = .(geo, year)]

save_with_meta(pop, "demo_r_pjangrp3_TOTAL_NUTS2_1995_2023")
message(sprintf("Saved population: %d rows", nrow(pop)))

message("Done: Eurostat extracts saved to data/raw/eurostat/.")

# =============================
# Harmonize to NUTS-2016 (crosswalk)
# =============================

# Re-read saved LE extract to decouple from workspace mutations
dem <- readRDS('data/raw/eurostat/demo_r_mlifexp_NUTS2_1995_2023.rds')

# If any older NUTS vintage appears, apply crosswalk at:
# data/raw/crosswalks/nuts_crosswalk.csv with cols:
# old_code,new_code,share,year_from,year_to
if (file.exists('data/raw/crosswalks/nuts_crosswalk.csv')) {
  cw <- data.table(readr::read_csv('data/raw/crosswalks/nuts_crosswalk.csv', show_col_types = FALSE))
  dem <- merge(dem, cw, by.x = 'geo', by.y = 'old_code', all.x = TRUE)
  dem[, geo_2016 := ifelse(is.na(new_code), geo, new_code)]
  # If splits exist, distribute value by share; else keep original
  dem[, value := ifelse(is.na(share), value, value * share)]
  dem[, geo := geo_2016][, c('new_code','geo_2016','share') := NULL]
}
saveRDS(dem, 'data/derived/le_panel_raw.rds')
data.table::fwrite(dem, 'data/derived/le_panel_raw.csv')
message('Saved harmonised LE panel to data/derived/le_panel_raw.*')

# --- Harmonize POP the same way and save  ← NEW ---
pop <- readRDS('data/raw/eurostat/demo_r_pjangrp3_TOTAL_NUTS2_1995_2023.rds')

if (file.exists('data/raw/crosswalks/nuts_crosswalk.csv')) {
  cw <- data.table(readr::read_csv('data/raw/crosswalks/nuts_crosswalk.csv', show_col_types = FALSE))
  pop <- merge(pop, cw, by.x = 'geo', by.y = 'old_code', all.x = TRUE)
  pop[, geo_2016 := ifelse(is.na(new_code), geo, new_code)]
  pop[, pop := ifelse(is.na(share), pop, pop * share)]
  pop[, geo := geo_2016][, c('new_code','geo_2016','share') := NULL]
}
saveRDS(pop, 'data/derived/pop_panel_raw.rds')
data.table::fwrite(pop, 'data/derived/pop_panel_raw.csv')
message('Saved harmonised POP panel to data/derived/pop_panel_raw.*')

message("All done: raw extracts written with metadata; LE & POP harmonized to NUTS-2016.")
