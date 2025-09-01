# BETA CONVERGENCE DIAGNOSTICS
# Run this after your main script to investigate the large beta coefficient

# -------- Load the processed data --------
# Assuming dt is still in memory from your main script
# If not, you'll need to reload and reprocess the data

cat("=== DETAILED DIAGNOSTICS FOR LARGE BETA ===\n\n")

# 1. Examine the relationship visually
library(ggplot2)

# Sample a subset for visualization (too many points otherwise)
set.seed(123)
dt_sample <- dt[sample(.N, min(1000, .N))]

# Scatterplot: growth vs log(initial LE)
p1 <- ggplot(dt_sample, aes(x = log_le, y = g)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(x = "Log(Life Expectancy)", y = "Growth Rate",
       title = "Raw Relationship: Growth vs Log(Initial LE)") +
  theme_minimal()

print(p1)
ggsave("figs/diagnostic_growth_vs_loglevel.png", p1, width = 8, height = 6)

# 2. Check correlation
cor_raw <- cor(dt$log_le, dt$g, use = "complete.obs")
cat(sprintf("Raw correlation between log(LE) and growth: %.4f\n", cor_raw))

# 3. Alternative specifications to compare
cat("\n=== TESTING ALTERNATIVE SPECIFICATIONS ===\n")

# Spec 1: No fixed effects (raw correlation)
mod_no_fe <- lm(g ~ log_le, data = dt, weights = dt$w_pop)
beta_no_fe <- coef(mod_no_fe)["log_le"]
cat(sprintf("Beta with NO fixed effects: %.6f\n", beta_no_fe))

# Spec 2: Only year fixed effects
library(fixest)
mod_year_fe <- feols(g ~ log_le | year, data = dt, cluster = ~ geo, weights = ~ w_pop)
beta_year_fe <- coef(mod_year_fe)["log_le"]
cat(sprintf("Beta with YEAR fixed effects only: %.6f\n", beta_year_fe))

# Spec 3: Only geo fixed effects  
mod_geo_fe <- feols(g ~ log_le | geo, data = dt, cluster = ~ geo, weights = ~ w_pop)
beta_geo_fe <- coef(mod_geo_fe)["log_le"]
cat(sprintf("Beta with GEO fixed effects only: %.6f\n", beta_geo_fe))

# Spec 4: Alternative functional form - levels instead of logs
dt[, level_le := value]  # Use levels instead of logs
mod_levels <- feols(g ~ level_le | geo + year, data = dt, cluster = ~ geo, weights = ~ w_pop)
beta_levels <- coef(mod_levels)["level_le"]
cat(sprintf("Beta using LE LEVELS (not logs): %.6f\n", beta_levels))

# 5. Check if the issue is with very small variation in log(LE)
le_stats <- dt[, .(
  le_min = min(value, na.rm = TRUE),
  le_max = max(value, na.rm = TRUE), 
  le_range = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
  le_cv = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE),
  log_le_min = min(log_le, na.rm = TRUE),
  log_le_max = max(log_le, na.rm = TRUE),
  log_le_range = max(log_le, na.rm = TRUE) - min(log_le, na.rm = TRUE),
  log_le_cv = sd(log_le, na.rm = TRUE) / mean(log_le, na.rm = TRUE)
)]

cat("\n=== LIFE EXPECTANCY VARIATION ===\n")
print(le_stats)

# 6. Check for different patterns by time period
dt[, period := ifelse(year < 2008, "Pre-crisis", 
                      ifelse(year < 2020, "Post-crisis", "COVID"))]

# Estimate beta by period
beta_by_period <- dt[, {
  if (.N > 100) {
    mod <- tryCatch(
      feols(g ~ log_le | geo + year, data = .SD, cluster = ~ geo, weights = ~ w_pop),
      error = function(e) NULL
    )
    if (!is.null(mod)) {
      data.table(beta = coef(mod)["log_le"], n = .N)
    } else {
      data.table(beta = NA_real_, n = .N)
    }
  } else {
    data.table(beta = NA_real_, n = .N)
  }
}, by = period]

cat("\n=== BETA BY TIME PERIOD ===\n")
print(beta_by_period)

# 7. Try a different convergence specification: Gap-to-mean approach
cat("\n=== GAP-TO-MEAN SPECIFICATION ===\n")
dt[, le_gap := value - eu_mean]  # Gap in levels
mod_gap <- feols(g ~ le_gap | geo + year, data = dt, cluster = ~ geo, weights = ~ w_pop)
beta_gap <- coef(mod_gap)["le_gap"]
cat(sprintf("Beta for gap-to-mean (levels): %.6f\n", beta_gap))

# 8. Check what happens with de-meaned specification
dt[, log_le_demeaned := log_le - mean(log_le, na.rm = TRUE), by = year]
mod_demean <- feols(g ~ log_le_demeaned | geo + year, data = dt, cluster = ~ geo, weights = ~ w_pop)
beta_demean <- coef(mod_demean)["log_le_demeaned"]
cat(sprintf("Beta with year-demeaned log(LE): %.6f\n", beta_demean))

# 9. Summary of all specifications
cat("\n=== SUMMARY OF ALL SPECIFICATIONS ===\n")
spec_summary <- data.table(
  Specification = c("No FE", "Year FE only", "Geo FE only", "Geo + Year FE (original)", 
                    "Levels (not logs)", "Gap-to-mean", "Year-demeaned"),
  Beta = c(beta_no_fe, beta_year_fe, beta_geo_fe, -0.396200, 
           beta_levels, beta_gap, beta_demean),
  Realistic = c("No", "No", "No", "No", "TBD", "TBD", "TBD")
)
print(spec_summary)

# 10. Recommended specification based on literature
cat("\n=== RECOMMENDED FIXES ===\n")
cat("Based on the diagnostics above:\n")
cat("1. If gap-to-mean gives reasonable results, use that specification\n")
cat("2. Consider using LE levels instead of logs if variation is too small\n") 
cat("3. The issue may be that log(LE) variation is too compressed\n")
cat("4. Try longer time differences (e.g., 5-year growth) to get more variation\n")

# 11. Quick test with 5-year differences
cat("\n=== TESTING 5-YEAR DIFFERENCES ===\n")
setorder(dt, geo, year)
dt[, value_5y := shift(value, n = 5, type = "lead"), by = geo]
dt[, year_5y := shift(year, n = 5, type = "lead"), by = geo]
dt5 <- dt[!is.na(value_5y) & (year_5y - year) == 5]
dt5[, g_5y := ((value_5y / value)^(1/5) - 1)]  # Annualized 5-year growth
dt5[, log_le_5y := log(value)]

if (nrow(dt5) > 100) {
  mod_5y <- feols(g_5y ~ log_le_5y | geo + year, data = dt5, cluster = ~ geo, weights = ~ w_pop)
  beta_5y <- coef(mod_5y)["log_le_5y"] 
  half_life_5y <- if(beta_5y < 0) -log(2)/log(1 + beta_5y) else NA
  cat(sprintf("Beta with 5-year growth: %.6f (half-life: %.1f years)\n", beta_5y, half_life_5y))
}

cat("\n=== DIAGNOSIS COMPLETE ===\n")
cat("Look at the results above to identify which specification gives realistic beta coefficients.\n")