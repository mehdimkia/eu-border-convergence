# 05_event_study.R — Border × Year contrast (EU27-only, plot included)
# Usage: source("code/05_event_study.R") from the repo root

# Optional: auto-detect project root
if (requireNamespace("rprojroot", quietly = TRUE)) {
  root <- tryCatch(
    rprojroot::find_root(rprojroot::has_file(".git") | rprojroot::has_file("README.md")),
    error = function(e) getwd()
  )
  setwd(root)
}

source("code/00_setup.R")
library(data.table)
library(fixest)
library(ggplot2)
library(broom)
library(readr)

# Load data
le <- readRDS("data/derived/le_panel_raw.rds")
flags <- data.table::fread("data/derived/border_flags.csv")
setDT(le); setDT(flags)

# Merge and restrict to EU27
le <- merge(le, flags, by.x = "geo", by.y = "NUTS_ID", all.x = TRUE)
le <- le[substr(geo, 1, 2) %in% EU27]

# Quick diagnostic: how many NA flags remain?
na_ct <- le[is.na(border), .N]
if (na_ct > 0) message("Note: ", na_ct, " rows have NA border flag after EU27 filter (will be dropped).")

# Keep total sex
dt <- le[sex == "T"]

# Remove rows with missing RHS components
dt[, post_year := year - 2019]
dt <- dt[!is.na(border) & !is.na(post_year)]

# TWFE (interpretive), cluster by region
esmod <- feols(value ~ i(post_year, border, ref = -1) | geo + year, data = dt, cluster = ~geo)
saveRDS(esmod, "outputs/eventstudy_total.rds")
fixest::etable(esmod, file = "outputs/eventstudy_total.txt")

# Extract coefficients and 95% CI
coefs <- broom::tidy(esmod)
coefs <- coefs[grepl("^post_year::", coefs$term), ]
coefs$year_rel <- as.integer(gsub("post_year::(-?\\d+).*", "\\1", coefs$term))
coefs$ci_lo <- coefs$estimate - 1.96 * coefs$std.error
coefs$ci_hi <- coefs$estimate + 1.96 * coefs$std.error
readr::write_csv(coefs, "outputs/eventstudy_coefs.csv")

# Plot
p <- ggplot(coefs, aes(year_rel, estimate)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2) +
  geom_line() +
  geom_point() +
  labs(x = "Years relative to 2019", y = "Border − Non-border (LE difference)",
       title = "Event-style contrast (EU27; ref=2019; region & year FE)") +
  theme_pub()

ggplot2::ggsave("figs/F3_event_study.png", p, width = 7, height = 4, dpi = 300)

message("Done: event-study model, coefficients, and plot exported.")
