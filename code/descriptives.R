# 03_descriptives.R — panel audit & descriptives for LE convergence project
# Run from repo root

# ---- Project root (same as your other scripts) ----
if (requireNamespace("rprojroot", quietly = TRUE)) {
  root <- tryCatch(
    rprojroot::find_root(rprojroot::has_file(".git") | rprojroot::has_file("README.md")),
    error = function(e) getwd()
  )
  setwd(root)
}
dir.create("outputs", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

# ---- Setup ----
source("code/00_setup.R")
library(data.table)
library(ggplot2)
theme_set(theme_minimal(base_size = 11))

# ---- Load & merge ----
le <- readRDS("data/derived/le_panel_raw.rds")
flags <- data.table::fread("data/derived/border_flags.csv")
setDT(le); setDT(flags)
le <- merge(le, flags, by.x = "geo", by.y = "NUTS_ID", all.x = TRUE)

# Optional population weights
pop_path <- "data/derived/pop_panel.rds"
has_pop <- file.exists(pop_path)
if (has_pop) {
  pop <- readRDS(pop_path); setDT(pop)
  pop <- pop[, .(geo, year, pop = as.numeric(pop))]
  le <- merge(le, pop, by = c("geo","year"), all.x = TRUE)
}

# ---- Filters: EU-27, LE at birth ----
le <- le[substr(geo, 1, 2) %in% EU27]
if ("age" %in% names(le)) le <- le[age == "Y_LT1"]

# ---- Data dictionary (what variables we use) ----
dict <- data.table::data.table(
  variable = c("geo","year","sex","value","border","pop"),
  type     = c("factor","int","factor","numeric","0/1","numeric"),
  meaning  = c("NUTS-2 region code","calendar year","M/F/T",
               "life expectancy at birth (years)",
               "1 if border region (GISCO)","population (if available)")
)
fwrite(dict, "outputs/descriptives_data_dictionary.csv")

# ---- Panel audit: counts & coverage ----
panel_counts <- le[, .(
  n = .N,
  years_min = min(year), years_max = max(year),
  n_years = uniqueN(year)
), by = .(geo, sex)]
fwrite(panel_counts, "outputs/descriptives_panel_counts_by_geo_sex.csv")

by_year <- le[, .(n_regions = uniqueN(geo),
                  n_obs = .N), by = .(year, sex)]
fwrite(by_year, "outputs/descriptives_counts_by_year_sex.csv")

# ---- Lead/step & growth g ----
setorder(le, geo, sex, year)
le[, value_lead := shift(value, type = "lead"), by = .(geo, sex)]
le[, year_lead  := shift(year,  type = "lead"), by = .(geo, sex)]
le[, step := year_lead - year]

step_tab <- le[, .N, by = step][order(step)]
fwrite(step_tab, "outputs/descriptives_step_distribution.csv")

# Keep valid transitions and compute annualized growth
le_valid <- le[step >= 1 & is.finite(value) & is.finite(value_lead)]
le_valid[, g := (value_lead - value) / step]

# Winsorize g within (year, sex) at 0.5% tails (same as main spec)
le_valid[, `:=`(
  ql = quantile(g, 0.005, na.rm = TRUE),
  qh = quantile(g, 0.995, na.rm = TRUE)
), by = .(year, sex)]
le_valid[g < ql, g := ql]
le_valid[g > qh, g := qh]
le_valid[, c("ql","qh") := NULL]

# ---- EU mean & baseline (unweighted + weighted if pop available) ----
eu_mean_unw <- le_valid[, .(eu_mean = mean(value, na.rm = TRUE)), by = .(year, sex)]
le_valid <- merge(le_valid, eu_mean_unw, by = c("year","sex"))
le_valid[, baseline_unw := value - eu_mean]

if (has_pop) {
  eu_mean_w <- le_valid[, .(
    eu_mean_w = weighted.mean(value, w = pop, na.rm = TRUE)
  ), by = .(year, sex)]
  le_valid <- merge(le_valid, eu_mean_w, by = c("year","sex"))
  le_valid[, baseline_w := value - eu_mean_w]
}

# ---- Descriptive summaries ----
# 1) LE by year (mean, sd, min, max), Total sex primary
le_year_T <- le_valid[sex == "T",
                      .(mean_LE = mean(value), sd_LE = sd(value),
                        p10 = quantile(value, .10), p50 = median(value),
                        p90 = quantile(value, .90), min = min(value), max = max(value)),
                      by = year]
fwrite(le_year_T, "outputs/descriptives_LE_by_year_T.csv")

# 2) Growth g distribution (before/after winsor — here only after)
g_summ <- le_valid[, .(
  n = .N,
  mean = mean(g), sd = sd(g),
  p01 = quantile(g, .01), p05 = quantile(g, .05),
  p50 = median(g), p95 = quantile(g, .95), p99 = quantile(g, .99),
  min = min(g), max = max(g)
), by = sex][order(sex)]
fwrite(g_summ, "outputs/descriptives_g_summary_by_sex.csv")

# 3) Baseline distribution around EU mean
base_cols <- c("baseline_unw", if (has_pop) "baseline_w" else NULL)
base_long <- melt(le_valid[, c("year","sex", base_cols), with = FALSE],
                  id.vars = c("year","sex"), variable.name = "baseline_type", value.name = "baseline")
base_summ <- base_long[, .(
  mean = mean(baseline), sd = sd(baseline),
  p10 = quantile(baseline, .10), p50 = median(baseline),
  p90 = quantile(baseline, .90), min = min(baseline), max = max(baseline)
), by = .(sex, baseline_type)]
fwrite(base_summ, "outputs/descriptives_baseline_summary.csv")

# 4) Weight coverage
if (has_pop) {
  w_cov <- le_valid[, .(
    share_with_pop = mean(is.finite(pop)),
    sum_pop = sum(pop, na.rm = TRUE)
  ), by = year]
  fwrite(w_cov, "outputs/descriptives_weight_coverage_by_year.csv")
}

# ---- Quick plots (appendix-quality) ----
# Histogram of g (Total)
p1 <- ggplot(le_valid[sex=="T"], aes(g)) +
  geom_histogram(bins = 60) +
  labs(title = "Distribution of annualized growth g (Total)", x = "g (years per year)", y = "Count")
ggsave("figs/D1_hist_g_T.png", p1, width = 7, height = 4, dpi = 300)

# Histogram of baseline (unweighted)
p2 <- ggplot(le_valid[sex=="T"], aes(baseline_unw)) +
  geom_histogram(bins = 60) +
  labs(title = "Baseline (EU mean gap), unweighted (Total)", x = "LE - EU mean (years)", y = "Count")
ggsave("figs/D2_hist_baseline_unw_T.png", p2, width = 7, height = 4, dpi = 300)

# EU mean LE trend (unweighted, with weighted overlay if available)
p3 <- ggplot(unique(le_valid[sex=="T", .(year, eu_mean)]), aes(year, eu_mean)) +
  geom_line() +
  labs(title = "EU mean life expectancy (unweighted, Total)", x = NULL, y = "Years")
ggsave("figs/D3_eu_mean_trend_unw_T.png", p3, width = 7, height = 4, dpi = 300)

if (has_pop) {
  p3w <- ggplot(unique(le_valid[sex=="T", .(year, eu_mean_w)]), aes(year, eu_mean_w)) +
    geom_line() +
    labs(title = "EU mean life expectancy (population-weighted, Total)", x = NULL, y = "Years")
  ggsave("figs/D3b_eu_mean_trend_wpop_T.png", p3w, width = 7, height = 4, dpi = 300)
}

# σ series (unweighted; weighted if available)
wsd <- function(x, w = NULL) {
  if (is.null(w)) {
    x <- x[is.finite(x)]
    return(sd(x))
  } else {
    w <- w / sum(w, na.rm = TRUE)
    m <- sum(w * x, na.rm = TRUE)
    return(sqrt(sum(w * (x - m)^2, na.rm = TRUE)))
  }
}

sig_unw <- le_valid[sex=="T", .(sigma = wsd(log(value))), by = year]
fwrite(sig_unw, "outputs/descriptives_sigma_unweighted.csv")
p4 <- ggplot(sig_unw, aes(year, sigma)) + geom_line() +
  labs(title = "σ-convergence (SD of log LE) — unweighted, Total", y = "SD(log LE)", x = NULL)
ggsave("figs/D4_sigma_series_unw_T.png", p4, width = 7, height = 4, dpi = 300)

if (has_pop) {
  sig_w <- le_valid[sex=="T", .(sigma = wsd(log(value), w = pop)), by = year]
  fwrite(sig_w, "outputs/descriptives_sigma_weighted.csv")
  p5 <- ggplot(sig_w, aes(year, sigma)) + geom_line() +
    labs(title = "σ-convergence (SD of log LE) — weighted, Total", y = "SD(log LE)", x = NULL)
  ggsave("figs/D5_sigma_series_wpop_T.png", p5, width = 7, height = 4, dpi = 300)
}

# Boxplot of g by year (Total) — shows COVID tails visually
p6 <- ggplot(le_valid[sex=="T"], aes(factor(year), g)) +
  geom_boxplot(outlier.alpha = .2) +
  labs(title = "Annualized growth g by year (Total)", x = "Year", y = "g (years per year)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("figs/D6_boxplot_g_by_year_T.png", p6, width = 9, height = 4.8, dpi = 300)

message("Descriptives exported: check outputs/ and figs/ for CSVs and PNGs.")
