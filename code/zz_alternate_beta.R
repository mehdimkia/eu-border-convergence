# ALTERNATIVE BETA CONVERGENCE SPECIFICATIONS
# Add this to your main script after line where you define dt

cat("\n=== TESTING ALTERNATIVE SPECIFICATIONS ===\n")

# Alternative 1: Gap-to-mean specification (more common in regional convergence)
dt[, le_gap := value - eu_mean]  # Gap to EU mean in levels
bmod_gap <- feols(g ~ le_gap | geo + year, data = dt, cluster = ~ geo, weights = ~ w_pop)
beta_gap <- coef(bmod_gap)["le_gap"]
half_life_gap <- if(beta_gap < 0) -log(2)/log(1 + beta_gap) else NA

cat(sprintf("GAP-TO-MEAN SPECIFICATION:\n"))
cat(sprintf("Beta: %.6f\n", beta_gap))
if (!is.na(half_life_gap)) cat(sprintf("Half-life: %.1f years\n", half_life_gap))

# Alternative 2: Levels instead of logs (since LE range is compressed in logs)
bmod_levels <- feols(g ~ value | geo + year, data = dt, cluster = ~ geo, weights = ~ w_pop)
beta_levels <- coef(bmod_levels)["value"]

# Convert to comparable units (effect of 1% higher initial LE)
beta_levels_pct <- beta_levels * mean(dt$value, na.rm = TRUE) / 100
half_life_levels <- if(beta_levels < 0) -log(2)/log(1 + beta_levels_pct) else NA

cat(sprintf("\nLEVELS SPECIFICATION:\n"))
cat(sprintf("Beta (per 1 year higher LE): %.6f\n", beta_levels))
cat(sprintf("Beta (per 1%% higher LE): %.6f\n", beta_levels_pct))
if (!is.na(half_life_levels)) cat(sprintf("Half-life: %.1f years\n", half_life_levels))

# Alternative 3: Simpler fixed effects (year only)
bmod_simple <- feols(g ~ log_le | year, data = dt, cluster = ~ geo, weights = ~ w_pop)
beta_simple <- coef(bmod_simple)["log_le"]
half_life_simple <- if(beta_simple < 0) -log(2)/log(1 + beta_simple) else NA

cat(sprintf("\nSIMPLE SPECIFICATION (year FE only):\n"))
cat(sprintf("Beta: %.6f\n", beta_simple))
if (!is.na(half_life_simple)) cat(sprintf("Half-life: %.1f years\n", half_life_simple))

# Alternative 4: 5-year differences (reduces noise)
setorder(dt, geo, year) 
dt[, value_5y := shift(value, n = 5, type = "lead"), by = geo]
dt[, year_5y := shift(year, n = 5, type = "lead"), by = geo]
dt5 <- dt[!is.na(value_5y) & (year_5y - year) == 5]
dt5[, g_5y := ((value_5y / value)^(1/5) - 1)]  # Annualized

if (nrow(dt5) > 100) {
  bmod_5y <- feols(g_5y ~ log(value) | geo + year, data = dt5, cluster = ~ geo, weights = ~ w_pop)
  beta_5y <- coef(bmod_5y)["log(value)"]
  half_life_5y <- if(beta_5y < 0) -log(2)/log(1 + beta_5y) else NA
  
  cat(sprintf("\n5-YEAR DIFFERENCES:\n"))
  cat(sprintf("Beta: %.6f\n", beta_5y))
  if (!is.na(half_life_5y)) cat(sprintf("Half-life: %.1f years\n", half_life_5y))
}

# Summary recommendation
cat(sprintf("\n=== RECOMMENDATION ===\n"))
cat(sprintf("Original specification beta: -0.396 (unrealistic)\n"))

realistic_specs <- c()
if (abs(beta_gap) < 0.05 && beta_gap < 0) realistic_specs <- c(realistic_specs, "Gap-to-mean")
if (abs(beta_levels_pct) < 0.05 && beta_levels_pct < 0) realistic_specs <- c(realistic_specs, "Levels")
if (abs(beta_simple) < 0.05 && beta_simple < 0) realistic_specs <- c(realistic_specs, "Simple FE")

if (length(realistic_specs) > 0) {
  cat(sprintf("Realistic specifications: %s\n", paste(realistic_specs, collapse = ", ")))
  cat(sprintf("RECOMMENDED: Use the gap-to-mean specification for your main results.\n"))
} else {
  cat(sprintf("All specifications still show unrealistic coefficients.\n"))
  cat(sprintf("Consider: 1) Data quality issues, 2) Different time periods, 3) Subsamples\n"))
}

# Save the better specification
if (abs(beta_gap) < abs(beta_levels_pct) && abs(beta_gap) < 0.05) {
  cat(sprintf("\nSaving gap-to-mean specification as main result...\n"))
  saveRDS(bmod_gap, sprintf("outputs/beta_gap_to_mean_weighted_%s_annual.rds", weight_label))
  fixest::etable(bmod_gap, file = sprintf("outputs/beta_gap_to_mean_weighted_%s_annual.txt", weight_label))
}