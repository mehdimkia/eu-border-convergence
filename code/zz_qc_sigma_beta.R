# zz_qc_sigma_beta.R — QA for 04_sigma_beta.R outputs
suppressPackageStartupMessages({
  library(data.table); library(readr); library(fixest)
})

# ---- Params (match your 04_sigma_beta.R toggles) ----
USE_ANCHOR_WEIGHTS <- TRUE
ANCHOR_YEAR <- 2019
K_YEARS <- 5
weight_label <- if (USE_ANCHOR_WEIGHTS) sprintf("anchor%d", ANCHOR_YEAR) else "yearly"

dir.create("outputs/QA", recursive = TRUE, showWarnings = FALSE)

# ---- Load data/flags ----
le   <- readRDS("data/derived/le_panel_raw.rds") |> as.data.table()
pop  <- readRDS("data/derived/pop_panel_raw.rds") |> as.data.table()
flags <- fread("data/derived/border_flags.csv")

data.table::setnames(flags, tolower(names(flags)))  # <- normalize case
stopifnot(all(c("nuts_id","eu_member") %in% names(flags)))
flags <- flags[, .(geo = nuts_id, eu_member)]

# EU filter via flag (handles EL)
le <- merge(le, flags, by = "geo", all.x = TRUE)
le <- le[eu_member == TRUE]
if ("age" %in% names(le)) le <- le[age %in% c("Y_LT1","Y0")]
le[, value := as.numeric(value)]

# ---- Check EL (Greece) presence in the filtered panel ----
has_EL <- any(substr(le$geo, 1, 2) == "EL")
cat(sprintf("EU filter includes Greece (EL): %s\n", has_EL))

# ---- Build weights (same as main script) ----
if (USE_ANCHOR_WEIGHTS) {
  pa <- copy(pop)[!is.na(pop)]
  pa[, dist := abs(year - ANCHOR_YEAR)]
  setorder(pa, geo, dist, year)
  pa <- pa[, .SD[1L], by = geo][, .(geo, w_pop = pop)]
  le <- merge(le, pa, by = "geo", all.x = TRUE)
  # fallback if any NA
  if (anyNA(le$w_pop)) {
    fallback <- pop[, .(w_pop = mean(pop, na.rm = TRUE)), by = geo]
    le <- merge(le, fallback, by = "geo", all.x = TRUE, suffixes = c("", "_fb"))
    le[, w_pop := fifelse(is.na(w_pop), w_pop_fb, w_pop)]
    le[, w_pop_fb := NULL]
  }
} else {
  le <- merge(le, pop[, .(geo, year, w_pop = pop)], by = c("geo","year"), all.x = TRUE)
}
le[, w_pop := as.numeric(w_pop)]

# ---- Weight sanity ----
w_na <- le[!is.finite(w_pop) | w_pop <= 0]
fwrite(w_na, "outputs/QA/weights_missing_or_nonpos.csv")
cat(sprintf("Weights — non-finite or <=0 rows: %d\n", nrow(w_na)))

# ---- σ-series sanity ----
sig_path <- sprintf("outputs/sigma_series_wpop_%s.csv", weight_label)
if (file.exists(sig_path)) {
  sig <- fread(sig_path)
  sig_bad <- sig[!is.finite(sigma) | sigma <= 0]
  fwrite(sig_bad, "outputs/QA/sigma_nonfinite_or_nonpos.csv")
  cat(sprintf("Sigma — years: %d | bad rows (non-finite/<=0): %d | range: [%.5f, %.5f]\n",
              nrow(sig), nrow(sig_bad), min(sig$sigma, na.rm = TRUE), max(sig$sigma, na.rm = TRUE)))
} else {
  cat("Sigma — MISSING file; did 04_sigma_beta.R run?\n")
}

# Piecewise output presence + finiteness
pw_tbl <- sprintf("outputs/sigma_piecewise_2020_wpop_%s.csv", weight_label)
pw_sum <- sprintf("outputs/sigma_piecewise_2020_summary_wpop_%s.csv", weight_label)
if (file.exists(pw_tbl) && file.exists(pw_sum)) {
  ct <- fread(pw_tbl)
  sm <- fread(pw_sum)
  ok_ct <- all(is.finite(ct$estimate)) && all(is.finite(ct$std_error))
  ok_sm <- all(is.finite(sm$estimate)) && all(is.finite(sm$se))
  cat(sprintf("Sigma piecewise — coef table finite: %s | summary finite: %s\n", ok_ct, ok_sm))
} else {
  cat("Sigma piecewise — MISSING outputs\n")
}

# ---- Annual β: recompute dt and compare to saved model nobs ----
setorder(le, geo, sex, year)
le[, value_lead := shift(value, type = "lead"), by = .(geo, sex)]
le[, year_lead  := shift(year,  type = "lead"), by = .(geo, sex)]
le[, step := year_lead - year]
le <- le[step == 1]
le[, g := (value_lead - value) / step]
le <- le[is.finite(g)]

# Winsorize within (year × sex)
le[, c("lo","hi") := as.list(quantile(g, c(.005,.995), na.rm = TRUE)), by = .(year, sex)]
le[g < lo, g := lo]; le[g > hi, g := hi]; le[, c("lo","hi") := NULL]

# EU-mean baseline
eu_mean <- le[, .(eu_mean = weighted.mean(value, w = w_pop, na.rm = TRUE)), by = .(year, sex)]
le <- merge(le, eu_mean, by = c("year","sex"))
le[, baseline := value - eu_mean]

dt <- le[sex == "T" & is.finite(g) & is.finite(baseline) & is.finite(w_pop) & w_pop > 0]
cat(sprintf("Annual β — constructed sample size (T only): %d\n", nrow(dt)))

mod_path <- sprintf("outputs/beta_eu_mean_Tonly_weighted_%s_annual.rds", weight_label)
if (file.exists(mod_path)) {
  m <- readRDS(mod_path)
  cat(sprintf("Annual β — model nobs: %d | coef(beta baseline): %.6f\n",
              nobs(m), coef(m)["baseline"]))
  if (nobs(m) != nrow(dt)) {
    fwrite(dt[1:200], "outputs/QA/beta_annual_dt_head.csv")
    cat("Annual β — WARNING: model nobs != constructed dt rows. See outputs/QA/beta_annual_dt_head.csv\n")
  }
} else {
  cat("Annual β — MISSING model rds\n")
}

# ---- Rolling β (K-year): existence and sanity ----
roll_csv <- sprintf("outputs/beta_rolling_Tonly_wpop_%s_K%d.csv", weight_label, K_YEARS)
if (file.exists(roll_csv)) {
  rd <- fread(roll_csv)
  bad <- rd[!is.finite(beta) | !is.finite(se) | is.na(center)]
  cat(sprintf("Rolling β — centers: %d | bad rows: %d | beta range: [%.4f, %.4f]\n",
              nrow(rd), nrow(bad),
              ifelse(nrow(rd)>0, min(rd$beta, na.rm = TRUE), NA_real_),
              ifelse(nrow(rd)>0, max(rd$beta, na.rm = TRUE), NA_real_)))
  fwrite(bad, "outputs/QA/beta_rolling_bad_rows.csv")
} else {
  cat("Rolling β — MISSING CSV (maybe not enough years for rolling window)\n")
}

cat("QA complete. See outputs/QA/ for any detailed listings.\n")
