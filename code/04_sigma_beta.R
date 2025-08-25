# 04_sigma_beta.R — σ and β (EU-27; LE at birth; ANNUALIZED growth; weighted primary)
# Adds anchor-year weight option + strict annual β + 5-year forward-growth rolling β
# Usage: source("code/04_sigma_beta.R") from the repo root

# -------- Toggles --------
USE_ANCHOR_WEIGHTS <- TRUE   # if TRUE, use population from ANCHOR_YEAR for all years
ANCHOR_YEAR <- 2019          # anchor year to use for weights
ROLLING_USE_K_YEARS <- TRUE  # if TRUE, rolling β uses K-year forward growth (annualized)
K_YEARS <- 5                 # forward horizon for rolling β (if enabled)

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
dir.create("figs", showWarnings = FALSE, recursive = TRUE)

# -------- Load panels & flags --------
le    <- readRDS("data/derived/le_panel_raw.rds")
flags <- data.table::fread("data/derived/border_flags.csv")
pop   <- readRDS("data/derived/pop_panel_raw.rds")

setDT(le); setDT(flags); setDT(pop)

# Merge flags + population
le <- merge(le, flags, by.x = "geo", by.y = "NUTS_ID", all.x = TRUE)
le <- merge(le, pop,  by = c("geo","year"), all.x = TRUE)  # 'pop' here is year-specific; we may replace with anchor weights below

# Keep EU-27 & life expectancy at birth
le <- le[substr(geo, 1, 2) %in% EU27]
if ("age" %in% names(le)) le <- le[age == "Y_LT1"]

# Sanity conversions
le[, value := as.numeric(value)]
le[, pop   := as.numeric(pop)]

# -------- Build weight vector --------
# If using anchor-year weights: for each geo, take population from ANCHOR_YEAR
# If missing, take the nearest available year to ANCHOR_YEAR (tie -> earliest)
if (USE_ANCHOR_WEIGHTS) {
  pa <- copy(pop)[!is.na(pop)]
  pa[, dist := abs(year - ANCHOR_YEAR)]
  setorder(pa, geo, dist, year)
  pa <- pa[, .SD[1L], by = geo]
  setnames(pa, "pop", "w_pop")
  pa[, year := NULL]; pa[, dist := NULL]
  le <- merge(le, pa, by = "geo", all.x = TRUE)
  # Guard: if any w_pop still NA, fall back to mean(pop) across years for that geo (rare)
  if (anyNA(le$w_pop)) {
    fallback <- pop[, .(w_pop = mean(pop, na.rm = TRUE)), by = geo]
    le <- merge(le, fallback, by = "geo", all.x = TRUE, suffixes = c("", "_fb"))
    le[, w_pop := fifelse(is.na(w_pop), w_pop_fb, w_pop)]
    le[, w_pop_fb := NULL]
  }
  weight_label <- paste0("anchor", ANCHOR_YEAR)
} else {
  le[, w_pop := pop]
  weight_label <- "yearly"
}

stopifnot(!all(is.na(le$w_pop)))
le[, w_pop := as.numeric(w_pop)]

# -------- Weighted σ-convergence (log-LE; T only) --------
wsd <- function(x, w) {
  w <- if (missing(w) || is.null(w)) rep(1/length(x), length(x)) else w/sum(w, na.rm = TRUE)
  mu <- sum(w * x, na.rm = TRUE)
  sqrt(sum(w * (x - mu)^2, na.rm = TRUE))
}

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
ggplot2::ggsave(sprintf("figs/F2_sigma_series_wpop%s.png",
                        ifelse(USE_ANCHOR_WEIGHTS, "", "_yearly")),
                g_sigma, width = 7, height = 4, dpi = 300)
# Also write legacy filename used in manuscript (points to current weight choice)
ggplot2::ggsave("figs/F2_sigma_series_wpop.png", g_sigma, width = 7, height = 4, dpi = 300)

# -------- Segmented σ (change in slope after 2020, robust SE) --------
sigma_series[, `:=`(post = pmax(0, year - 2020),
                    year_c = year - 1995)]
m_pw <- lm(sigma ~ year_c + post, data = sigma_series)

# Robust vcov and tidy table (manual to avoid dimnames issues)
vcNW <- NeweyWest(m_pw, lag = 2)
terms <- names(coef(m_pw))
est   <- as.numeric(coef(m_pw))
se    <- sqrt(diag(vcNW))
stat  <- est / se
pval  <- 2*pt(abs(stat), df = stats::df.residual(m_pw), lower.tail = FALSE)

ct_df <- data.frame(
  term     = terms,
  estimate = est,
  std_error= se,
  statistic= stat,
  p_value  = pval,
  stringsAsFactors = FALSE
)
write.csv(ct_df, sprintf("outputs/sigma_piecewise_2020_wpop_%s.csv", weight_label), row.names = FALSE)

# Also export pre/post/Δ slopes with 95% CIs
pre  <- unname(coef(m_pw)["year_c"])
post <- unname(coef(m_pw)["year_c"] + coef(m_pw)["post"])
del  <- unname(coef(m_pw)["post"])

se_pre  <- sqrt(vcNW["year_c","year_c"])
se_post <- sqrt(vcNW["year_c","year_c"] + vcNW["post","post"] + 2*vcNW["year_c","post"])
se_del  <- sqrt(vcNW["post","post"])

sum_df <- data.table(
  metric   = c("pre_2020_slope","post_2020_slope","delta_slope"),
  estimate = c(pre, post, del),
  se       = c(se_pre, se_post, se_del)
)[, `:=`(
  ci_lo = estimate - 1.96*se,
  ci_hi = estimate + 1.96*se
)]
fwrite(sum_df, sprintf("outputs/sigma_piecewise_2020_summary_wpop_%s.csv", weight_label))

# -------- β-convergence (annualized growth; weighted primary) --------
# Prepare panel for growth and baselines — STRICTLY ANNUAL
setorder(le, geo, sex, year)
le[, value_lead := shift(value, type = "lead"), by = .(geo, sex)]
le[, year_lead  := shift(year,  type = "lead"), by = .(geo, sex)]
le[, step := year_lead - year]
le <- le[step == 1]  # enforce strictly annual steps

# Annualized growth (annual diffs)
le[, g := (value_lead - value) / step]
le <- le[is.finite(g)]

# Winsorize within (year × sex) at 0.5% / 99.5%
le[, c("lo","hi") := as.list(quantile(g, c(.005,.995), na.rm = TRUE)), by = .(year, sex)]
le[g < lo, g := lo]
le[g > hi, g := hi]
le[, c("lo","hi") := NULL]

# Population-weighted EU mean baseline by (year × sex); weights per toggle
eu_mean <- le[, .(eu_mean = weighted.mean(value, w = w_pop, na.rm = TRUE)),
              by = .(year, sex)]
le <- merge(le, eu_mean, by = c("year","sex"))
le[, baseline := value - eu_mean]

# ---- Choose ONE β sample (T only is primary) ----
SAMPLE <- quote(sex == "T")

dt <- le[eval(SAMPLE) &
         is.finite(g) & is.finite(baseline) &
         is.finite(w_pop) & w_pop > 0]

# TWFE model with region & year FE; cluster by region; weights = w_pop (ANNUAL)
bmod_w_annual <- feols(g ~ baseline | geo + year, data = dt, cluster = ~geo, weights = ~w_pop)

# Guardrail
stopifnot(nobs(bmod_w_annual) == nrow(dt))

# Save model and a text table (annual)
saveRDS(bmod_w_annual, sprintf("outputs/beta_eu_mean_Tonly_weighted_%s_annual.rds", weight_label))
fixest::etable(bmod_w_annual, file = sprintf("outputs/beta_eu_mean_Tonly_weighted_%s_annual.txt", weight_label))

# Debug stamp
cat(sprintf("Run: %s | sample=Tonly_wpop_%s_ANNUAL | n=%d | beta_baseline=%.6f\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            weight_label, nrow(dt), coef(bmod_w_annual)["baseline"]),
    file = sprintf("outputs/beta_debug_Tonly_wpop_%s_annual.txt", weight_label))

# -------- Rolling β using K-year forward growth (annualized) --------
if (ROLLING_USE_K_YEARS) {
  # Build K-year forward growth on the original (pre-step-filter) panel
  leK <- readRDS("data/derived/le_panel_raw.rds")
  setDT(leK); leK <- merge(leK, flags[, .(NUTS_ID, CNTR_CODE)], by.x="geo", by.y="NUTS_ID", all.x=TRUE)
  leK <- leK[substr(geo, 1, 2) %in% EU27]
  if ("age" %in% names(leK)) leK <- leK[age == "Y_LT1"]
  leK <- merge(leK, if (USE_ANCHOR_WEIGHTS) pa[, .(geo, w_pop)] else pop[, .(geo, year, pop)][, .(geo, w_pop = pop)], by="geo", all.x=TRUE)
  leK[, value := as.numeric(value)]
  setorder(leK, geo, sex, year)
  leK[, value_lead_k := shift(value, n = K_YEARS, type = "lead"), by = .(geo, sex)]
  leK[, year_lead_k  := shift(year,  n = K_YEARS, type = "lead"), by = .(geo, sex)]
  leK[, step_k := year_lead_k - year]
  dtk <- leK[sex=="T" & step_k == K_YEARS]
  dtk[, gk := (value_lead_k - value) / K_YEARS]
  dtk <- dtk[is.finite(gk) & is.finite(w_pop) & w_pop > 0]
  # Winsorize by start-year
  dtk[, c("lo","hi") := as.list(quantile(gk, c(.005,.995), na.rm = TRUE)), by = .(year)]
  dtk[gk < lo, gk := lo]; dtk[gk > hi, gk := hi]; dtk[, c("lo","hi") := NULL]
  # Weighted EU mean at start year
  eu_mean_k <- dtk[, .(eu_mean = weighted.mean(value, w = w_pop, na.rm = TRUE)), by = .(year, sex)]
  dtk <- merge(dtk, eu_mean_k, by = c("year","sex"))
  dtk[, baseline := value - eu_mean]
  # Full-sample K-year β (optional to report)
  bmod_k <- feols(gk ~ baseline | geo + year, data = dtk, cluster = ~geo, weights = ~w_pop)
  saveRDS(bmod_k, sprintf("outputs/beta_eu_mean_Tonly_weighted_%s_K%d.rds", weight_label, K_YEARS))
  fixest::etable(bmod_k, file = sprintf("outputs/beta_eu_mean_Tonly_weighted_%s_K%d.txt", weight_label, K_YEARS))
  # Rolling windows on start year
  yrs <- sort(unique(dtk$year))
  if (length(yrs) >= 7) {
    centers <- yrs[yrs >= (min(yrs) + 2) & yrs <= (max(yrs) - 2)]
    roll_dt <- rbindlist(lapply(centers, function(cy){
      sub <- dtk[year %between% c(cy - 2, cy + 2)]
      if (nrow(sub) < 50) return(data.table(center = cy, beta = NA_real_, se = NA_real_, n = nrow(sub)))
      fit <- tryCatch(
        feols(gk ~ baseline | geo + year, data = sub, cluster = ~geo, weights = ~w_pop),
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
           y = "β (annualized; years/year per 1-year gap)",
           title = sprintf("β-convergence (EU-27; LE at birth; weighted, %s; %d-year forward growth) — T only",
                           weight_label, K_YEARS)) +
      theme_pub()
    ggplot2::ggsave(sprintf("figs/F1_beta_rolling_Tonly_wpop_%s_K%d.png", weight_label, K_YEARS),
                    g_beta, width = 7, height = 4, dpi = 300)
    # Also write legacy filenames to point to the K-year figure by default
    ggplot2::ggsave("figs/F1_beta_rolling_Tonly_wpop.png", g_beta, width = 7, height = 4, dpi = 300)
    ggplot2::ggsave("figs/F1_beta_rolling.png",          g_beta, width = 7, height = 4, dpi = 300)
  } else {
    message("Not enough years for rolling β_K; skipped rolling plot.")
  }
} else {
  message("ROLLING_USE_K_YEARS == FALSE: using annual β for rolling (not recommended).")
}
# -------- Provenance note --------
writeLines(sprintf("Weights: %s (anchor year %s = %s)\nβ primary: strict annual steps; Rolling β: %s (K=%d)\nTimestamp: %s",
                   ifelse(USE_ANCHOR_WEIGHTS, "anchor-year", "year-specific"),
                   ANCHOR_YEAR, weight_label,
                   ifelse(ROLLING_USE_K_YEARS, "forward growth", "annual"),
                   K_YEARS,
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
           sprintf("outputs/weights_info_%s.txt", weight_label))

message("Done: weighted σ (anchor), β annual (strict), and rolling β (K-year) exported (weights = ", weight_label, ").")
