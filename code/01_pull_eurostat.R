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
# ---------- Harmonize to NUTS-2016 (correct handling) ----------
# Re-read raw extracts
dem <- readRDS("data/raw/eurostat/demo_r_mlifexp_NUTS2_1995_2023.rds")
pop <- readRDS("data/raw/eurostat/demo_r_pjangrp3_TOTAL_NUTS2_1995_2023.rds")

if (file.exists("data/raw/crosswalks/nuts_crosswalk.csv")) {
  cw <- data.table::fread("data/raw/crosswalks/nuts_crosswalk.csv", showProgress = FALSE)
  setnames(cw, tolower(names(cw)))
  if (!"year_from" %in% names(cw)) cw[, year_from := -Inf]
  if (!"year_to"   %in% names(cw)) cw[, year_to   :=  Inf]
  
  # Identity rows for unmatched codes
  id_codes <- setdiff(unique(dem$geo), unique(cw$old_code))
  if (length(id_codes)) {
    cw <- rbind(cw, data.table(old_code = id_codes, new_code = id_codes, share = 1, year_from = -Inf, year_to = Inf), fill = TRUE)
  }
  
  # Map POP (counts): allocate by share, respect year window, then sum
  pop_map <- merge(pop, cw, by.x = "geo", by.y = "old_code", allow.cartesian = TRUE)
  pop_map <- pop_map[year >= year_from & year <= year_to]
  pop_map[, geo := new_code][, pop := pop * share]
  pop_h <- pop_map[, .(pop = sum(pop, na.rm = TRUE)), by = .(geo, year)]
  
  # Map LE (rates): do NOT scale; later pop-weighted average
  dem_map <- merge(dem, cw, by.x = "geo", by.y = "old_code", allow.cartesian = TRUE)
  dem_map <- dem_map[year >= year_from & year <= year_to]
  dem_map[, geo := new_code]
  setnames(dem_map, "value", "le")
  
  # Attach the mapped population weights at the same split granularity
  # (match by original geo-year-new_code tuple)
  pop_w <- pop_map[, .(old_geo = geo, year, new_code, pop_split = pop, geo = new_code)]
  dem_w <- merge(
    dem_map[, .(old_geo = geo, year, sex, le, new_code)],
    pop_w, by = c("old_geo", "year", "new_code"), all.x = TRUE, allow.cartesian = TRUE
  )
  
  # Weighted collapse to unique (geo=new_code, year, sex)
  dem_h <- dem_w[, .(
    value = sum(le * pop_split, na.rm = TRUE) / sum(pop_split, na.rm = TRUE)
  ), by = .(geo = new_code, year, sex)]
  
} else {
  # No crosswalk: just rename and keep unique
  dem_h <- as.data.table(dem)[, .(value = value), by = .(geo, year, sex)]
  pop_h <- as.data.table(pop)[, .(pop = pop), by = .(geo, year)]
}

# Save derived panels
saveRDS(dem_h, "data/derived/le_panel_raw.rds")
data.table::fwrite(dem_h, "data/derived/le_panel_raw.csv")
saveRDS(pop_h, "data/derived/pop_panel_raw.rds")
data.table::fwrite(pop_h, "data/derived/pop_panel_raw.csv")
message("Saved harmonised LE (pop-weighted) and POP panels to data/derived/.")