# 05_event_study.R — Border × Year contrast (EU-27; LE at birth)
# Usage: source("code/05_event_study.R") from the repo root

# -------- Toggles --------
WEIGHTED <- TRUE                  # primary: weighted by population
USE_ANCHOR_WEIGHTS <- TRUE        # match primary spec used in 04_sigma_beta.R
ANCHOR_YEAR <- 2019
EXCLUDE_EU_NON_EU_TOUCH <- FALSE  # sensitivity: exclude regions touching non-EU neighbors

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
  library(broom)
  library(readr)
})

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("figs",    showWarnings = FALSE, recursive = TRUE)

# -------- Load & merge --------
le    <- readRDS("data/derived/le_panel_raw.rds")   |> as.data.table()
pop   <- readRDS("data/derived/pop_panel_raw.rds")  |> as.data.table()
flags <- data.table::fread("data/derived/border_flags.csv")

# normalise flag names (NUTS_ID vs nuts_id, etc.)
data.table::setnames(flags, tolower(names(flags)))
stopifnot(all(c("nuts_id","border","eu_member") %in% names(flags)))

# Merge flags + population (year-specific pop now; anchor handled below)
le <- merge(le, flags, by.x = "geo", by.y = "nuts_id", all.x = TRUE)
le <- merge(le, pop,   by = c("geo","year"),          all.x = TRUE)

# Restrict to EU-27 via flag, LE at birth, total sex
le <- le[eu_member == TRUE]
if ("age" %in% names(le)) le <- le[age %in% c("Y_LT1","Y0")]
dt <- le[sex == "T"]

# Optional exclusion: regions touching non-EU borders
if (EXCLUDE_EU_NON_EU_TOUCH && "eu_noneu_touch" %in% names(dt)) {
  dt <- dt[eu_noneu_touch == FALSE]
}

# Design: dynamic diff-in-diff style around 2019
dt[, post_year := year - 2019]
dt <- dt[!is.na(border) & !is.na(post_year)]

# Weights (anchor-year or yearly)
suffix_w <- if (WEIGHTED) "_Tonly_wpop" else "_Tonly_unw"
weight_label <- if (WEIGHTED && USE_ANCHOR_WEIGHTS) sprintf("_%d", ANCHOR_YEAR) else ""
if (WEIGHTED) {
  # build w: anchor-year or yearly
  if (USE_ANCHOR_WEIGHTS) {
    pa <- copy(pop)[!is.na(pop)]
    pa[, dist := abs(year - ANCHOR_YEAR)]
    setorder(pa, geo, dist, year)
    pa <- pa[, .SD[1L], by = geo][, .(geo, w = as.numeric(pop))]
    dt <- merge(dt, pa, by = "geo", all.x = TRUE)
  } else {
    dt[, w := fifelse(is.finite(pop) & pop > 0, as.numeric(pop), NA_real_)]
  }
  if (any(is.na(dt$w))) {
    message("NOTE: ", sum(is.na(dt$w)), " observations have NA weights and will be dropped. ",
            "Consider anchor-year weights if you need full pre-period coverage.")
  }
} else {
  dt[, w := 1.0]
}

# TWFE with region & year FE; cluster by region
# IMPORTANT: ref = 0 makes 2019 the reference year (correct labeling)
esmod <- feols(value ~ i(post_year, border, ref = 0) | geo + year,
               data = dt, cluster = ~ geo, weights = ~ w)

saveRDS(esmod, paste0("outputs/eventstudy_total", suffix_w, weight_label, ".rds"))
fixest::etable(esmod, file = paste0("outputs/eventstudy_total", suffix_w, weight_label, ".txt"))

# Extract coefficients and 95% CI for plotting
coefs <- broom::tidy(esmod)
coefs <- coefs[grepl("^post_year::", coefs$term), ]
coefs$year_rel <- as.integer(gsub("post_year::(-?\\d+).*", "\\1", coefs$term))
coefs$ci_lo <- coefs$estimate - 1.96 * coefs$std.error
coefs$ci_hi <- coefs$estimate + 1.96 * coefs$std.error  # <-- fixed (no extra 1.96)
readr::write_csv(coefs, paste0("outputs/eventstudy_coefs", suffix_w, weight_label, ".csv"))

# Pre-trend joint Wald test: all pre-2019 coefficients = 0 (cluster-robust)
pre_terms <- coefs$term[coefs$year_rel < 0]
pre_terms <- intersect(pre_terms, names(coef(esmod)))
if (length(pre_terms) > 0) {
  b <- coef(esmod)[pre_terms]
  V <- vcov(esmod, cluster = "geo")[pre_terms, pre_terms, drop = FALSE]
  # numerical stability
  invV <- tryCatch(chol2inv(chol(V)), error = function(e) MASS::ginv(V))
  stat <- as.numeric(t(b) %*% invV %*% b)
  df   <- length(pre_terms)
  pval <- pchisq(stat, df = df, lower.tail = FALSE)
  out <- data.frame(n_pre = df, statistic = stat, df = df, p_value = pval,
                    terms = paste(pre_terms, collapse = "; "))
  readr::write_csv(out, paste0("outputs/eventstudy_pretrend_wald", suffix_w, weight_label, ".csv"))
}

# Plot
p <- ggplot(coefs, aes(year_rel, estimate)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2) +
  geom_line() + geom_point() +
  labs(x = "Years relative to 2019",
       y = "Border − Non-border (LE, years)",
       title = paste0("Event-style contrast (EU-27; LE at birth; ref = 2019; ",
                      ifelse(WEIGHTED, paste0("weighted", weight_label), "unweighted"), ")")) +
  theme_pub()

ggplot2::ggsave(paste0("figs/F3_event_study", suffix_w, weight_label, ".png"),
                p, width = 7, height = 4, dpi = 300)

message("Done: event-study model, coefficients, pre-trend test, and plot exported: ",
        suffix_w, weight_label)
