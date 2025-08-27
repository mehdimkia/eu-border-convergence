# 04_sigma_beta_DEBUG_FIX.R
# Debugging script to identify why beta is too large

# =============================================================================
# SETUP
# =============================================================================
library(data.table)
library(fixest)

# Load data
le <- readRDS("data/derived/le_panel_raw.rds")
flags <- fread("data/derived/border_flags.csv")
setDT(le); setDT(flags)

# Merge and filter
le <- merge(le, flags, by.x = "geo", by.y = "NUTS_ID", all.x = TRUE)
le <- le[eu_member == TRUE]
if ("age" %in% names(le)) le <- le[age %in% c("Y_LT1", "Y0")]

# Calculate growth (strictly annual)
setorder(le, geo, sex, year)
le[, value_lead := shift(value, type = "lead"), by = .(geo, sex)]
le[, year_lead := shift(year, type = "lead"), by = .(geo, sex)]
le[, step := year_lead - year]
le <- le[step == 1]

# =============================================================================
# CALCULATE DIFFERENT GROWTH SPECIFICATIONS
# =============================================================================
cat("=== DEBUGGING BETA CONVERGENCE ===\n\n")

# 1. ABSOLUTE GROWTH (your current method)
le[, g_absolute := value_lead - value]
cat("Absolute growth (years per year):\n")
print(summary(le$g_absolute))

# 2. PERCENTAGE GROWTH (more standard)
le[, g_percent := 100 * (value_lead - value) / value]
cat("\nPercentage growth (% per year):\n")
print(summary(le$g_percent))

# 3. LOG GROWTH (econometrically preferred)
le[, g_log := log(value_lead) - log(value)]
cat("\nLog growth (≈ % for small changes):\n")
print(summary(le$g_log))

# =============================================================================
# CALCULATE DIFFERENT BASELINE SPECIFICATIONS
# =============================================================================

# Get EU means (unweighted for simplicity in debug)
eu_mean <- le[, .(eu_mean = mean(value, na.rm = TRUE)), by = .(year, sex)]
le <- merge(le, eu_mean, by = c("year", "sex"))

# 1. ABSOLUTE BASELINE (your current method)
le[, baseline_absolute := value - eu_mean]
cat("\n\nAbsolute baseline (years from EU mean):\n")
print(summary(le$baseline_absolute))

# 2. PERCENTAGE BASELINE
le[, baseline_percent := 100 * (value - eu_mean) / eu_mean]
cat("\nPercentage baseline (% from EU mean):\n")
print(summary(le$baseline_percent))

# 3. LOG BASELINE
le[, baseline_log := log(value) - log(eu_mean)]
cat("\nLog baseline:\n")
print(summary(le$baseline_log))

# =============================================================================
# RUN REGRESSIONS WITH DIFFERENT SPECIFICATIONS
# =============================================================================
cat("\n\n=== REGRESSION RESULTS (Total sex only) ===\n")

# Filter to total sex
dt <- le[sex == "T"]
results <- list()

# Define specifications
specs <- list(
  list(g = "g_absolute", b = "baseline_absolute", label = "Absolute/Absolute (CURRENT)"),
  list(g = "g_percent", b = "baseline_percent", label = "Percent/Percent"),
  list(g = "g_log", b = "baseline_log", label = "Log/Log"),
  list(g = "g_percent", b = "baseline_absolute", label = "Percent/Absolute"),
  list(g = "g_absolute", b = "baseline_percent", label = "Absolute/Percent")
)

# Run each specification
for (i in seq_along(specs)) {
  spec <- specs[[i]]
  
  # Create regression data
  reg_data <- dt[is.finite(get(spec$g)) & is.finite(get(spec$b))]
  
  # Skip if not enough data
  if (nrow(reg_data) < 100) next
  
  # Run regression
  formula <- as.formula(paste(spec$g, "~", spec$b, "| geo + year"))
  mod <- feols(formula, data = reg_data, cluster = ~geo)
  
  # Extract results
  beta <- coef(mod)[spec$b]
  se <- sqrt(vcov(mod)[spec$b, spec$b])
  
  # Calculate implied half-life (years to halve the gap)
  if (spec$g == "g_log" && spec$b == "baseline_log") {
    # For log-log: beta is elasticity
    half_life <- -log(2) / beta
  } else if (spec$g == "g_percent" && spec$b == "baseline_percent") {
    # For percent-percent: beta is % reduction per 1% gap
    half_life <- -log(2) / (beta / 100)
  } else if (spec$g == "g_absolute" && spec$b == "baseline_absolute") {
    # For absolute-absolute: rough approximation
    # Assume average gap is ~2 years
    avg_gap <- 2
    half_life <- (avg_gap / 2) / abs(beta)
  } else {
    half_life <- NA
  }
  
  cat(sprintf("\n%s:\n", spec$label))
  cat(sprintf("  Beta = %.4f (SE = %.4f)\n", beta, se))
  cat(sprintf("  Interpretation: %s\n", 
              if (spec$g == "g_absolute" && spec$b == "baseline_absolute") {
                sprintf("1-year gap → %.3f years/year convergence", abs(beta))
              } else if (spec$g == "g_percent" && spec$b == "baseline_percent") {
                sprintf("1%% gap → %.3f%% annual convergence", abs(beta))
              } else if (spec$g == "g_log" && spec$b == "baseline_log") {
                sprintf("1%% gap → %.3f%% annual convergence (elasticity)", abs(beta) * 100)
              } else {
                "Mixed units - interpretation unclear"
              }))
  
  if (!is.na(half_life) && half_life > 0) {
    cat(sprintf("  Implied half-life: %.1f years\n", half_life))
  }
  
  results[[spec$label]] <- list(beta = beta, se = se, n = nrow(reg_data))
}

# =============================================================================
# SHOW THE PROBLEM WITH ABSOLUTE/ABSOLUTE
# =============================================================================
cat("\n\n=== WHY ABSOLUTE/ABSOLUTE IS PROBLEMATIC ===\n")

# Example calculation
avg_le <- 80  # typical life expectancy
gap <- 2      # 2-year gap from mean
beta_abs <- -0.411

cat(sprintf("\nExample: Region with LE = %.0f (%.0f years below EU mean of %.0f)\n", 
            avg_le - gap, gap, avg_le))
cat(sprintf("Your beta = %.3f implies:\n", beta_abs))
cat(sprintf("  Annual growth boost = %.3f * %.0f = %.3f years/year\n", 
            abs(beta_abs), gap, abs(beta_abs) * gap))
cat(sprintf("  Time to close gap = %.0f / %.3f = %.1f years\n", 
            gap, abs(beta_abs) * gap, gap / (abs(beta_abs) * gap)))

cat("\nThis is impossibly fast! Real convergence takes decades, not 2-3 years.\n")

# =============================================================================
# RECOMMENDED SPECIFICATION
# =============================================================================
cat("\n=== RECOMMENDED APPROACH ===\n")
cat("Use log-log specification:\n")
cat("  g = log(LE_t+1) - log(LE_t)\n")
cat("  baseline = log(LE_t) - log(EU_mean_t)\n")
cat("\nThis gives beta as an elasticity (% convergence per 1% gap)\n")
cat("Typical values: beta ∈ [-0.01, -0.05], implying 1-5% convergence per year\n")

# =============================================================================
# CORRECT YOUR ORIGINAL CODE
# =============================================================================
cat("\n=== TO FIX YOUR CODE ===\n")
cat("Replace these lines:\n")
cat('  le[, g := (value_lead - value) / step]\n')
cat('  le[, baseline := value - eu_mean]\n')
cat("\nWith:\n")
cat('  le[, g := log(value_lead) - log(value)]  # log growth\n')
cat('  le[, baseline := log(value) - log(eu_mean)]  # log baseline\n')
cat("\nOr for percentage specification:\n")
cat('  le[, g := 100 * (value_lead - value) / value]  # % growth\n')
cat('  le[, baseline := 100 * (value - eu_mean) / eu_mean]  # % baseline\n')

# Save results summary
sink("outputs/beta_specification_debug.txt")
cat("Beta Convergence Specification Analysis\n")
cat("=====================================\n\n")
for (i in seq_along(results)) {
  r <- results[[i]]
  cat(sprintf("%s: beta = %.4f (SE = %.4f), N = %d\n", 
              names(results)[i], r$beta, r$se, r$n))
}
sink()

cat("\n\nDebug complete. Check outputs/beta_specification_debug.txt\n")