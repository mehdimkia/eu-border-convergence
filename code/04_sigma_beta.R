# 04_sigma_beta.R — σ and β (EU-27; LE at birth; ANNUALIZED log growth; weighted primary)
# Usage: source("code/04_sigma_beta.R") from the repo root

# -------- Toggles --------
USE_ANCHOR_WEIGHTS <- TRUE    # if TRUE, use population from ANCHOR_YEAR for all years
ANCHOR_YEAR        <- 2019    # anchor year to use for weights
ROLLING_USE_K_YEARS <- TRUE   # if TRUE, rolling β uses K-year forward growth (annualized)
K_YEARS            <- 5       # forward horizon for rolling β (if enabled)
RUN_UNWEIGHTED_EU   <- TRUE   # write unweighted EU β (annual)
RUN_WITHIN_COUNTRY  <- TRUE   # write within-country β (annual, weighted)

# -------- Project root (deterministic) --------
if (requireNamespace("rprojroot", quietly = TRUE)) {
  root <- tryCatch(
    rprojroot::find_root(rprojroot::has_file("code/00_setup.R")),
    error = function(e) getwd()
  )
  setwd(root)
}

# -------- Env & packages --------
source("code/00_setup.R")
suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
  library(sandwich)
  library(lmtest)
})

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("figs",    showWarnings = FALSE, recursive = TRUE)

# -------- Load panels & flags --------
le    <- readRDS("data/derived/le_panel_raw.rds")
pop   <- readRDS("data/derived/pop_panel_raw.rds")
flags <- data.table::fread("data/derived/border_flags.csv")  # contains eu_member

setDT(le); setDT(pop); setDT(flags)

# Merge flags + population (year-specific pop; may be replaced by anchor weights below)
le <- merge(le, flags, by.x = "geo", by.y = "NUTS_ID", all.x = TRUE)
le <- merge(le, pop,   by = c("geo","year"),            all.x = TRUE)  # 'pop' is year-specific

# Keep EU-27 & life expectancy at birth (Eurostat uses Y_LT1 or sometimes Y0)
stopifnot("eu_member" %in% names(le))
le <- le[eu_member == TRUE]
if ("age" %in% names(le)) le <- le[age %in% c("Y_LT1","Y0")]

# Sanity conversions
le[, value := as.numeric(value)]
le[, pop   := as.numeric(pop)]

# -------- Build weights (anchor-year or yearly) --------
# If using anchor-year weights: for each geo, take population from ANCHOR_YEAR (nearest-year fallback)
weight_label <- if (USE_ANCHOR_WEIGHTS) sprintf("anchor%d", ANCHOR_YEAR) else "yearly"

if (USE_ANCHOR_WEIGHTS) {
  pa <- copy(pop)[!is.na(pop)]
  pa[, dist := abs(year - ANCHOR_YEAR)]
  setorder(pa, geo, dist, year)
  pa <- pa[, .SD[1L], by = geo]           # nearest-to-anchor per geo
  setnames(pa, "pop", "w_pop")
  pa[, `:=`(year = NULL, dist = NULL)]
  le <- merge(le, pa, by = "geo", all.x = TRUE)
  
  # Fallback: mean(pop) across available years for any remaining NA w_pop (rare)
  if (anyNA(le$w_pop)) {
    fallback <- pop[, .(w_pop = mean(pop, na.rm = TRUE)), by = geo]
    le <- merge(le, fallback, by = "geo", all.x = TRUE, suffixes = c("", "_fb"))
    le[, w_pop := fifelse(is.na(w_pop), w_pop_fb, w_pop)]
    le[, w_pop_fb := NULL]
  }
} else {
  le[, w_pop := pop]  # year-specific weights
}

stopifnot(!all(is.na(le$w_pop)))
le[, w_pop := as.numeric(w_pop)]

# -------- Helpers --------
wsd <- function(x, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  w <- w / sum(w)
  mu <- sum(w * x)
  sqrt(sum(w * (x - mu)^2))
}

theme_pub <- if (exists("theme_pub")) theme_pub else function() theme_minimal(base_size = 12)

# -------- Weighted σ-convergence (SD of log-LE; sex = T) --------
sigma_series <- le[sex == "T" & is.finite(value) & is.finite(w_pop) & w_pop > 0, {
  lx <- log(value)
  .(sigma = wsd(lx, w_pop), n = .N)
}, by = .(year)]

data.table::fwrite(sigma_series, sprintf("outputs/sigma_series_wpop_%s.csv", weight_label))

g_sigma <- ggplot(sigma_series, aes(year, sigma)) +
  geom_line() +
  labs(y = "SD(log LE)", x = NULL,
       title = sprintf("σ-convergence (EU-27; LE at birth; weighted, %s)", weight_label)) +
  theme_pub()
ggplot2::ggsave(sprintf("figs/F2_sigma_series_wpop_%s.png", weight_label), g_sigma, width = 7, height = 4, dpi = 300)
# Legacy filename used in manuscript → points to current weight choice
ggplot2::ggsave("figs/F2_sigma_series_wpop.png", g_sigma, width = 7, height = 4, dpi = 300)

# -------- Segmented σ (change in slope after 2020, robust SE) on ln(σ) --------
sigma_series[, `:=`(
  ln_sigma = log(sigma),
  post     = pmax(0, year - 2020),
  year_c   = year - 1995
)]
m_pw <- lm(ln_sigma ~ year_c + post, data = sigma_series)

# Newey–West (lag=2) robust vcov
vcNW  <- NeweyWest(m_pw, lag = 2)
co    <- coef(m_pw)
sevec <- sqrt(diag(vcNW))
stat  <- co / sevec
pval  <- 2 * pt(abs(stat), df = df.residual(m_pw), lower.tail = FALSE)

ct_df <- data.frame(
  term      = names(co),
  estimate  = unname(co),
  std_error = unname(sevec),
  statistic = unname(stat),
  p_value   = unname(pval),
  stringsAsFactors = FALSE
)
write.csv(ct_df, sprintf("outputs/sigma_piecewise_2020_wpop_%s.csv", weight_label), row.names = FALSE)

# Export pre/post/Δ slopes in % per year (on σ scale via ln-slope)
pre   <- unname(co["year_c"])
del   <- unname(co["post"])
post  <- pre + del
se_pre  <- sevec["year_c"]
se_del  <- sevec["post"]
se_post <- sqrt(se_pre^2 + se_del^2 + 2 * vcNW["year_c", "post"])

sum_df <- data.table(
  metric   = c("pre_2020_slope_pctpy","post_2020_slope_pctpy","delta_slope_pctpy"),
  estimate = 100 * c(pre, post, del),
  se       = 100 * c(se_pre, se_post, se_del)
)[, `:=`(ci_lo = estimate - 1.96 * se,
         ci_hi = estimate + 1.96 * se)]
fwrite(sum_df, sprintf("outputs/sigma_piecewise_2020_summary_wpop_%s.csv", weight_label))

# -------- β-convergence (annualized *log* growth; weighted primary) --------
# Prepare panel for growth and baselines — STRICTLY ANNUAL
setorder(le, geo, sex, year)
le[, value_lead := shift(value, type = "lead"), by = .(geo, sex)]
le[, year_lead  := shift(year,  type = "lead"), by = .(geo, sex)]
le[, step := year_lead - year]
le <- le[step == 1 & is.finite(value) & is.finite(value_lead) & value > 0 & value_lead > 0]

# Annualized *log* growth (g = [log(LE_{t+1}) - log(LE_t)] / 1)
le[, g := (log(value_lead) - log(value)) / step]
le <- le[is.finite(g)]

# Winsorize within (year × sex) at 0.5% / 99.5%
le[, c("lo","hi") := as.list(quantile(g, c(.005, .995), na.rm = TRUE)), by = .(year, sex)]
le[g <  lo, g := lo]
le[g >  hi, g := hi]
le[, c("lo","hi") := NULL]

# Population-weighted EU mean by (year × sex); weights per toggle
eu_mean <- le[, .(eu_mean = weighted.mean(value, w = w_pop, na.rm = TRUE)), by = .(year, sex)]
le <- merge(le, eu_mean, by = c("year","sex"), all.x = TRUE)

# Baseline regressor: log-level relative to EU-mean of same year
le[, baseline := log(value) - log(eu_mean)]

# ---- Choose ONE β sample (T only is primary) ----
SAMPLE <- quote(sex == "T")
dt <- le[eval(SAMPLE) & is.finite(g) & is.finite(baseline) & is.finite(w_pop) & w_pop > 0]

# TWFE with region & year FE; cluster by region; weights = w_pop
bmod_w_annual <- feols(g ~ baseline | geo + year, data = dt, cluster = ~ geo, weights = ~ w_pop)

# Guardrail: should match the sample size used
stopifnot(nobs(bmod_w_annual) == nrow(dt))

# Export baseline SD (per-1-SD effect helper)
sd_log_baseline <- wsd(dt$baseline, dt$w_pop)
data.table(sd_log_baseline = sd_log_baseline,
           n = nrow(dt)) |>
  fwrite(sprintf("outputs/beta_baseline_sd_Tonly_wpop_%s_annual.csv", weight_label))

# Save model and a text table (annual)
saveRDS(bmod_w_annual, sprintf("outputs/beta_eu_mean_Tonly_weighted_%s_annual.rds", weight_label))
fixest::etable(bmod_w_annual, file = sprintf("outputs/beta_eu_mean_Tonly_weighted_%s_annual.txt", weight_label))
# -------------------------------------------------------------------
# A) Unweighted EU-mean β (annual; same log spec, NO weights)
# -------------------------------------------------------------------
if (isTRUE(RUN_UNWEIGHTED_EU)) {
  bmod_unw_annual <- feols(g ~ baseline | geo + year, data = dt, cluster = ~ geo)
  saveRDS(bmod_unw_annual, "outputs/beta_eu_mean_Tonly_unweighted_annual.rds")
  fixest::etable(bmod_unw_annual,
                 file = "outputs/beta_eu_mean_Tonly_unweighted_annual.txt")
  cat(sprintf("UNWEIGHTED EU β: n=%d | beta=%.6f\n",
              nobs(bmod_unw_annual), coef(bmod_unw_annual)["baseline"]),
      file = "outputs/beta_debug_Tonly_unweighted_annual.txt")
}

# -------------------------------------------------------------------
# B) Within-country β (annual; weighted; country-demeaned baseline)
# -------------------------------------------------------------------
if (isTRUE(RUN_WITHIN_COUNTRY)) {
  # Country code (robust, with fallback)
  if ("CNTR_CODE" %in% names(le)) {
    le[, ctry := CNTR_CODE]
  } else if ("country" %in% names(le)) {
    le[, ctry := country]
  } else {
    le[, ctry := substr(geo, 1, 2)]  # NUTS country prefix fallback
  }
  
  # Country mean (weighted to match the primary spec)
  cty_mean <- le[, .(cty_mean = weighted.mean(value, w = w_pop, na.rm = TRUE)),
                 by = .(ctry, year, sex)]
  le <- merge(le, cty_mean, by = c("ctry","year","sex"), all.x = TRUE)
  
  # Country-demeaned log baseline
  le[, baseline_cty := log(value) - log(cty_mean)]
  
  dt_cty <- le[eval(SAMPLE) & is.finite(g) & is.finite(baseline_cty) &
                 is.finite(w_pop) & w_pop > 0]
  
  bmod_cty_w <- feols(g ~ baseline_cty | geo + year,
                      data = dt_cty, cluster = ~ geo, weights = ~ w_pop)
  
  saveRDS(bmod_cty_w, sprintf("outputs/beta_cty_mean_Tonly_weighted_%s_annual.rds", weight_label))
  fixest::etable(bmod_cty_w,
                 file = sprintf("outputs/beta_cty_mean_Tonly_weighted_%s_annual.txt", weight_label))
  
  cat(sprintf("WITHIN-COUNTRY (weighted) β: n=%d | beta=%.6f\n",
              nobs(bmod_cty_w), coef(bmod_cty_w)["baseline_cty"]),
      file = sprintf("outputs/beta_debug_cty_Tonly_wpop_%s_annual.txt", weight_label))
}



# Debug stamp
cat(sprintf("Run: %s | sample=Tonly_wpop_%s_ANNUAL | n=%d | beta_baseline(log-gap)=%.6f\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            weight_label, nrow(dt), coef(bmod_w_annual)["baseline"]),
    file = sprintf("outputs/beta_debug_Tonly_wpop_%s_annual.txt", weight_label))

# -------- Rolling β using K-year forward growth (annualized log growth) --------
bmod_k <- NULL
if (ROLLING_USE_K_YEARS) {
  # Build K-year forward growth on a fresh panel (avoid step==1 filter)
  leK <- readRDS("data/derived/le_panel_raw.rds"); setDT(leK)
  leK <- merge(leK, flags[, .(NUTS_ID, eu_member)], by.x = "geo", by.y = "NUTS_ID", all.x = TRUE)
  leK <- leK[eu_member == TRUE]
  if ("age" %in% names(leK)) leK <- leK[age %in% c("Y_LT1","Y0")]
  
  # Attach weights (anchor-year or yearly)
  if (USE_ANCHOR_WEIGHTS) {
    leK <- merge(leK, pa[, .(geo, w_pop)], by = "geo", all.x = TRUE)
  } else {
    leK <- merge(leK, pop[, .(geo, year, w_pop = pop)], by = c("geo","year"), all.x = TRUE)
  }
  
  leK[, value := as.numeric(value)]
  setorder(leK, geo, sex, year)
  
  # K-year lead
  leK[, value_lead_k := shift(value, n = K_YEARS, type = "lead"), by = .(geo, sex)]
  leK[, year_lead_k  := shift(year,  n = K_YEARS, type = "lead"), by = .(geo, sex)]
  leK[, step_k := year_lead_k - year]
  
  dtk <- leK[sex == "T" & step_k == K_YEARS & is.finite(value) & is.finite(value_lead_k) & value > 0 & value_lead_k > 0]
  
  # Annualized *log* growth over K years: gk = (log(LE_{t+K}) - log(LE_t)) / K
  dtk[, gk := (log(value_lead_k) - log(value)) / K_YEARS]
  dtk <- dtk[is.finite(gk) & is.finite(w_pop) & w_pop > 0]
  
  # Winsorize by start-year
  dtk[, c("lo","hi") := as.list(quantile(gk, c(.005, .995), na.rm = TRUE)), by = .(year)]
  dtk[gk <  lo, gk := lo]
  dtk[gk >  hi, gk := hi]
  dtk[, c("lo","hi") := NULL]
  
  # Weighted EU mean at start year, baseline log-gap
  eu_mean_k <- dtk[, .(eu_mean = weighted.mean(value, w = w_pop, na.rm = TRUE)), by = .(year, sex)]
  dtk <- merge(dtk, eu_mean_k, by = c("year","sex"), all.x = TRUE)
  dtk[, baseline := log(value) - log(eu_mean)]
  
  # Full-sample K-year β (optional to report)
  bmod_k <- feols(gk ~ baseline | geo + year, data = dtk, cluster = ~ geo, weights = ~ w_pop)
  saveRDS(bmod_k, sprintf("outputs/beta_eu_mean_Tonly_weighted_%s_K%d.rds", weight_label, K_YEARS))
  fixest::etable(bmod_k, file = sprintf("outputs/beta_eu_mean_Tonly_weighted_%s_K%d.txt", weight_label, K_YEARS))
  
  # Rolling windows on start year (±2 years around center)
  yrs <- sort(unique(dtk$year))
  if (length(yrs) >= 7) {
    centers <- yrs[yrs >= (min(yrs) + 2) & yrs <= (max(yrs) - 2)]
    roll_dt <- rbindlist(lapply(centers, function(cy) {
      sub <- dtk[year %between% c(cy - 2, cy + 2)]
      if (nrow(sub) < 50) return(data.table(center = cy, beta = NA_real_, se = NA_real_, n = nrow(sub)))
      fit <- tryCatch(
        feols(gk ~ baseline | geo + year, data = sub, cluster = ~ geo, weights = ~ w_pop),
        error = function(e) NULL
      )
      if (is.null(fit)) return(data.table(center = cy, beta = NA_real_, se = NA_real_, n = nrow(sub)))
      co <- coef(fit)["baseline"]
      se <- sqrt(diag(vcov(fit, cluster = "geo")))["baseline"]
      data.table(center = cy, beta = co, se = se, n = nrow(sub))
    }))
    roll_dt[, `:=`(ci_lo = beta - 1.96 * se,
                   ci_hi = beta + 1.96 * se)]
    data.table::fwrite(roll_dt, sprintf("outputs/beta_rolling_Tonly_wpop_%s_K%d.csv", weight_label, K_YEARS))
    
    g_beta <- ggplot(roll_dt, aes(center, beta)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_vline(xintercept = 2020, linetype = 3) +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2) +
      geom_line() + geom_point() +
      labs(x = "Start year (5-year estimation window)",
           y = "β (annualized log growth on log-level gap)",
           title = sprintf("β-convergence (EU-27; LE at birth; weighted, %s; %d-year forward growth) — T only",
                           weight_label, K_YEARS)) +
      theme_pub()
    ggplot2::ggsave(sprintf("figs/F1_beta_rolling_Tonly_wpop_%s_K%d.png", weight_label, K_YEARS),
                    g_beta, width = 7, height = 4, dpi = 300)
    # Legacy filenames for manuscript
    ggplot2::ggsave("figs/F1_beta_rolling_Tonly_wpop.png", g_beta, width = 7, height = 4, dpi = 300)
    ggplot2::ggsave("figs/F1_beta_rolling.png",          g_beta, width = 7, height = 4, dpi = 300)
  } else {
    message("Not enough years for rolling β_K; skipped rolling plot.")
  }
} else {
  message("ROLLING_USE_K_YEARS == FALSE: using annual β for rolling (not recommended).")
}

# -------- Guardrails (magnitude + consistency) --------
beta_a  <- tryCatch(coef(bmod_w_annual)["baseline"], error = function(e) NA_real_)
beta_k5 <- tryCatch(if (!is.null(bmod_k)) coef(bmod_k)["baseline"] else NA_real_, error = function(e) NA_real_)

if (is.finite(beta_a) && is.finite(beta_k5) && abs(beta_a - beta_k5) > 0.03) {
  warning(sprintf("β(K=1)=%0.3f and β(K=%d annualized)=%0.3f differ notably — check annualization.",
                  beta_a, K_YEARS, beta_k5))
}

half_life <- function(beta) log(0.5) / log(1 + beta)  # discrete-time half-life
if (is.finite(beta_a) && abs(beta_a) > 0.10) {
  warning(sprintf("Implausibly fast convergence: β=%0.3f → half-life ≈ %.1f years.", beta_a, half_life(beta_a)))
}

# -------- Provenance note --------
writeLines(sprintf(
  "Weights: %s (anchor year %s = %s)\nβ primary: strict annual steps (log growth on log-level gap)\nRolling β: forward log growth, annualized (K=%d)\nTimestamp: %s",
  ifelse(USE_ANCHOR_WEIGHTS, "anchor-year", "year-specific"),
  ANCHOR_YEAR, weight_label,
  K_YEARS,
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
), sprintf("outputs/weights_info_%s.txt", weight_label))

message("Done: weighted σ (anchor/yearly), β annual (strict, log), rolling β (K-year, log) — weights = ", weight_label, ".")
