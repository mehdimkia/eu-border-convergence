# 05_event_study.R — Border × Year contrast (EU27; LE at birth)
# Usage: source("code/05_event_study.R")

if (requireNamespace("rprojroot", quietly = TRUE)) {
  root <- tryCatch(
    rprojroot::find_root(rprojroot::has_file(".git") | rprojroot::has_file("README.md")),
    error = function(e) getwd()
  ); setwd(root)
}

source("code/00_setup.R")
library(data.table); library(fixest); library(ggplot2); library(broom); library(readr)

le <- readRDS("data/derived/le_panel_raw.rds")
flags <- data.table::fread("data/derived/border_flags.csv")
setDT(le); setDT(flags)
le <- merge(le, flags, by.x="geo", by.y="NUTS_ID", all.x=TRUE)

# Keep EU27, LE at birth, total sex
le <- le[substr(geo,1,2) %in% EU27]
if ("age" %in% names(le)) le <- le[age == "Y_LT1"]
dt <- le[sex=="T"]

dt[, post_year := year - 2019]
dt <- dt[!is.na(border) & !is.na(post_year)]

esmod <- feols(value ~ i(post_year, border, ref=-1) | geo + year, data = dt, cluster = ~geo)
saveRDS(esmod, "outputs/eventstudy_total.rds")
fixest::etable(esmod, file = "outputs/eventstudy_total.txt")

coefs <- broom::tidy(esmod)
coefs <- coefs[grepl("^post_year::", coefs$term), ]
coefs$year_rel <- as.integer(gsub("post_year::(-?\\d+).*", "\\1", coefs$term))
coefs$ci_lo <- coefs$estimate - 1.96 * coefs$std.error
coefs$ci_hi <- coefs$estimate + 1.96 * coefs$std.error
readr::write_csv(coefs, "outputs/eventstudy_coefs.csv")

p <- ggplot(coefs, aes(year_rel, estimate)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2) +
  geom_line() + geom_point() +
  labs(x = "Years relative to 2019", y = "Border − Non-border (LE, years)",
       title = "Event-style contrast (EU-27; LE at birth; ref=2019; region & year FE)") +
  theme_pub()
ggplot2::ggsave("figs/F3_event_study.png", p, width = 7, height = 4, dpi = 300)

message("Done: event-study for EU-27, LE at birth.")
