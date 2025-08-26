# zz_qc_harmonisation.R — QA for NUTS-2016 harmonisation
suppressPackageStartupMessages({
  library(data.table); library(readr); library(stringr)
})

dir.create("outputs/QA", recursive = TRUE, showWarnings = FALSE)

# ---- Paths
le_raw_path   <- "data/raw/eurostat/demo_r_mlifexp_NUTS2_1995_2023.rds"
pop_raw_path  <- "data/raw/eurostat/demo_r_pjangrp3_TOTAL_NUTS2_1995_2023.rds"
le_h_path     <- "data/derived/le_panel_raw.rds"
pop_h_path    <- "data/derived/pop_panel_raw.rds"
cw_path       <- "data/raw/crosswalks/nuts_crosswalk.csv"

# ---- Load
le_h  <- readRDS(le_h_path)   |> as.data.table()
pop_h <- readRDS(pop_h_path)  |> as.data.table()
le_r  <- readRDS(le_raw_path) |> as.data.table()
pop_r <- readRDS(pop_raw_path)|> as.data.table()

setnames(le_h,  tolower(names(le_h)))
setnames(pop_h, tolower(names(pop_h)))
setnames(le_r,  tolower(names(le_r)))
setnames(pop_r, tolower(names(pop_r)))

# ========== 1) UNIQUENESS ==========
dup_le  <- le_h[, .N, by = .(geo, year, sex)][N > 1]
dup_pop <- pop_h[, .N, by = .(geo, year)][N > 1]

fwrite(dup_le,  "outputs/QA/dup_le.csv")
fwrite(dup_pop, "outputs/QA/dup_pop.csv")

cat(sprintf("Uniqueness — LE dups: %d rows; POP dups: %d rows\n",
            nrow(dup_le), nrow(dup_pop)))

# ========== 2) NATIONAL POP TOTALS: BEFORE vs AFTER ==========
ctry_of <- function(x) substr(x, 1, 2)  # note: EL for Greece in Eurostat
pop_r[, ctry := ctry_of(geo)]
pop_h[, ctry := ctry_of(geo)]

tot_r <- pop_r[, .(pop_raw = sum(pop, na.rm = TRUE)), by = .(ctry, year)]
tot_h <- pop_h[, .(pop_h   = sum(pop, na.rm = TRUE)), by = .(ctry, year)]

cmp   <- merge(tot_r, tot_h, by = c("ctry","year"), all = TRUE)
cmp[, rel_diff := abs(pop_h - pop_raw) / pmax(pop_raw, 1)]
viol  <- cmp[is.finite(rel_diff) & rel_diff > 0.001]  # >0.1%

fwrite(cmp,  "outputs/QA/pop_totals_compare.csv")
fwrite(viol, "outputs/QA/pop_totals_violations.csv")

cat(sprintf("Totals — %d country-year pairs exceed 0.1%% difference\n", nrow(viol)))

# ========== 3) LE WEIGHTED-AVERAGE BOUNDS CHECK ==========
if (file.exists(cw_path)) {
  cw <- fread(cw_path, showProgress = FALSE)
  setnames(cw, tolower(names(cw)))
  if (!"year_from" %in% names(cw)) cw[, year_from := -Inf]
  if (!"year_to"   %in% names(cw)) cw[, year_to   :=  Inf]
  
  # Map raw POP with shares to get split weights (same as harmonisation)
  pop_map <- merge(pop_r, cw, by.x = "geo", by.y = "old_code", allow.cartesian = TRUE)
  pop_map <- pop_map[year >= year_from & year <= year_to]
  pop_map[, `:=`(pop_split = pop * share, new_geo = new_code)]
  
  # Map raw LE (do NOT scale), join weights at same split granularity
  le_map <- merge(le_r, cw, by.x = "geo", by.y = "old_code", allow.cartesian = TRUE)
  le_map <- le_map[year >= year_from & year <= year_to]
  le_map[, new_geo := new_code]
  setnames(le_map, "value", "le")
  
  w <- pop_map[, .(old_geo = geo, year, new_geo, pop_split)]
  le_w <- merge(
    le_map[, .(old_geo = geo, year, sex, le, new_geo)],
    w, by = c("old_geo","year","new_geo"), all.x = TRUE, allow.cartesian = TRUE
  )
  
  # Only cases where multiple old geos merge into same new (true aggregate)
  mult <- le_w[, .(n_src = uniqueN(old_geo)), by = .(new_geo, year, sex)][n_src > 1]
  chk  <- merge(le_w, mult, by = c("new_geo","year","sex"))
  
  # Source bounds and harmonised value
  src_bounds <- chk[, .(
    src_min = min(le, na.rm = TRUE),
    src_max = max(le, na.rm = TRUE)
  ), by = .(geo = new_geo, year, sex)]
  
  le_h_sub <- le_h[geo %in% src_bounds$geo]
  test <- merge(src_bounds, le_h_sub, by = c("geo","year","sex"), all.x = TRUE)  # value = LE_harmonised
  
  # Violations if harmonised value lies outside [min, max] by > 1e-6
  viol_le <- test[ (value < src_min - 1e-6) | (value > src_max + 1e-6) ]
  
  fwrite(test,    "outputs/QA/le_bounds_check_full.csv")
  fwrite(viol_le, "outputs/QA/le_bounds_violations.csv")
  
  cat(sprintf("LE bounds — %d merged (geo,year,sex) outside source min–max\n", nrow(viol_le)))
} else {
  cat("LE bounds — skipped (no crosswalk found)\n")
}

# ========== 4) YEAR-WINDOWS APPLIED ==========
if (file.exists(cw_path)) {
  # Identify any applied mappings outside their [year_from, year_to]
  # Build the mapping actually used
  used_map <- merge(le_r[, .(geo, year)], fread(cw_path)[, .(old_code, new_code, year_from, year_to)],
                    by.x = "geo", by.y = "old_code", allow.cartesian = TRUE)
  setnames(used_map, tolower(names(used_map)))
  out_of_window <- used_map[ year < year_from | year > year_to ]
  
  fwrite(out_of_window, "outputs/QA/year_window_violations.csv")
  cat(sprintf("Year windows — %d rows where a mapping would be applied outside its window\n",
              nrow(out_of_window)))
} else {
  cat("Year windows — skipped (no crosswalk found)\n")
}

cat("QA complete. See outputs/QA/ for details.\n")
