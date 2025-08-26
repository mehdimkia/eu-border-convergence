# 02_harmonise_nuts2016.R  — Harmonise LE & POP to NUTS-2016
source("code/00_setup.R")
suppressPackageStartupMessages({ library(data.table); library(readr) })

le_path  <- "data/raw/eurostat/demo_r_mlifexp_NUTS2_1995_2023.rds"
pop_path <- "data/raw/eurostat/demo_r_pjangrp3_TOTAL_NUTS2_1995_2023.rds"
cw_path  <- "data/raw/crosswalks/nuts_crosswalk.csv"

if (!file.exists(le_path))  stop("Missing LE file: ", le_path)
if (!file.exists(pop_path)) stop("Missing POP file: ", pop_path)

dem <- readRDS(le_path)   # expects columns: geo, year, sex, value (LE)
pop <- readRDS(pop_path)  # expects columns: geo, year, pop (counts)

# Ensure canonical column names
setDT(dem); setDT(pop)
setnames(dem, tolower(names(dem)))
setnames(pop, tolower(names(pop)))

# Optional sanity: ensure LE is at birth & units in years (if not enforced upstream)
if ("age" %in% names(dem)) dem <- dem[age %in% c("Y_LT1","Y0")]
if ("unit" %in% names(dem)) dem <- dem[unit %in% c("YR","Y")]

if (file.exists(cw_path)) {
  cw <- fread(cw_path, showProgress = FALSE)
  setnames(cw, tolower(names(cw)))
  if (!"year_from" %in% names(cw)) cw[, year_from := -Inf]
  if (!"year_to"   %in% names(cw)) cw[, year_to   :=  Inf]
  stopifnot(all(c("old_code","new_code","share","year_from","year_to") %in% names(cw)))
  
  # Add identity for any codes not present in the crosswalk
  id_codes <- setdiff(unique(dem$geo), unique(cw$old_code))
  if (length(id_codes)) {
    cw <- rbind(
      cw,
      data.table(old_code = id_codes, new_code = id_codes, share = 1, year_from = -Inf, year_to = Inf),
      fill = TRUE
    )
  }
  
  # --- Map POP (counts): allocate by share, respect year window, then sum
  pop_map <- merge(pop, cw, by.x = "geo", by.y = "old_code", allow.cartesian = TRUE)
  pop_map <- pop_map[year >= year_from & year <= year_to]
  pop_map[, `:=`(geo = new_code, pop_split = pop * share)]
  pop_h <- pop_map[, .(pop = sum(pop_split, na.rm = TRUE)), by = .(geo, year)]
  
  # --- Map LE (rates): do NOT scale; weight with mapped split-pop
  dem_map <- merge(dem, cw, by.x = "geo", by.y = "old_code", allow.cartesian = TRUE)
  dem_map <- dem_map[year >= year_from & year <= year_to]
  dem_map[, geo := new_code]
  setnames(dem_map, "value", "le")
  
  # Bring in the same split weights used for POP mapping
  w <- pop_map[, .(old_geo = geo, year, new_code, pop_split)]
  dem_w <- merge(
    dem_map[, .(old_geo = geo, year, sex, le, new_code)],
    w, by = c("old_geo","year","new_code"), all.x = TRUE, allow.cartesian = TRUE
  )
  
  # Collapse to unique (geo=new_code, year, sex) using population weights
  dem_h <- dem_w[, .(
    value = {
      sw <- sum(pop_split, na.rm = TRUE)
      if (is.finite(sw) && sw > 0) sum(le * pop_split, na.rm = TRUE) / sw else mean(le, na.rm = TRUE)
    }
  ), by = .(geo = new_code, year, sex)]
  
} else {
  # No crosswalk: just ensure uniqueness
  dem_h <- dem[, .(value = value), by = .(geo, year, sex)]
  pop_h <- pop[, .(pop = pop), by = .(geo, year)]
}

# Write derived panels
saveRDS(dem_h, "data/derived/le_panel_raw.rds")
fwrite(dem_h,   "data/derived/le_panel_raw.csv")
saveRDS(pop_h,  "data/derived/pop_panel_raw.rds")
fwrite(pop_h,   "data/derived/pop_panel_raw.csv")

message("Saved harmonised LE (pop-weighted) and POP panels to data/derived/.")
