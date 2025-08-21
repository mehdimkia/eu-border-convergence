# 04_sigma_beta.R — σ and β (EU-27; LE at birth; ANNUALIZED growth)
# Usage: source("code/04_sigma_beta.R") from the repo root

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
library(data.table)
library(fixest)
library(ggplot2)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("figs", showWarnings = FALSE, recursive = TRUE)

# -------- Load & merge --------
le <- readRDS("data/derived/le_panel_raw.rds")
flags <- data.table::fread("data/derived/border_flags.csv")
setDT(le); setDT(flags)
le <- merge(le, flags, by.x = "geo", by.y = "NUTS_ID", all.x = TRUE)

# Keep EU-27 & life expectancy at birth
le <- le[substr(geo, 1, 2) %in% EU27]
if ("age" %in% names(le)) le <- le[age == "Y_LT1"]

# -------- σ-convergence (log-LE; equal-weight for now) --------
wsd <- function(x, w) sqrt(sum(w * (x - sum(w*x))^2))
sigma_series <- le[sex == "T", {
  w <- rep(1/.N, .N)
  lx <- log(value)
  list(sigma = wsd(lx, w), n = .N)
}, by = .(year)]
data.table::fwrite(sigma_series, "outputs/sigma_series.csv")

g_sigma <- ggplot(sigma_series, aes(year, sigma)) +
  geom_line() +
  labs(y = "SD(log LE)", x = NULL, title = "σ-convergence (EU-27; LE at birth)") +
  theme_pub()
ggplot2::ggsave("figs/F2_sigma_series.png", g_sigma, width = 7, height = 4, dpi = 300)

# -------- β-convergence (annualized growth) --------
setorder(le, geo, sex, year)
le[, value_lead := shift(value, type = "lead"), by = .(geo, sex)]
le[, year_lead  := shift(year,  type = "lead"), by = .(geo, sex)]
le[, step := year_lead - year]
le <- le[step >= 1]
le[, g := (value_lead - value) / step]

# Trim extreme 0.5% tails for robustness
qs <- quantile(le$g, c(0.005, 0.995), na.rm = TRUE)
le[g < qs[1], g := qs[1]]
le[g > qs[2], g := qs[2]]

# EU-mean baseline by year × sex
eu_mean <- le[, .(eu_mean = mean(value, na.rm = TRUE)), by = .(year, sex)]
le <- merge(le, eu_mean, by = c("year","sex"))
le[, baseline := value - eu_mean]

# ---- Choose ONE β sample and reuse it everywhere ----
SAMPLE <- quote(sex == "T")                  # (A) Total only  ← main paper
# SAMPLE <- quote(sex %in% c("M","F","T"))   # (B) Pooled across sexes

dt <- le[eval(SAMPLE) & is.finite(g) & is.finite(baseline)]

# TWFE model with region & year FE; cluster by region
bmod <- feols(g ~ baseline | geo + year, data = dt, cluster = ~geo)

# Guardrails: saved and printed models must match exactly
stopifnot(nobs(bmod) == nrow(dt))

# Tag artifacts by sample to avoid mix-ups
tag <- if (identical(sort(unique(dt$sex)), "T")) "Tonly" else "AllSexes"

# Save model (RDS) and write a fresh text table (explicit overwrite)
saveRDS(bmod, file.path("outputs", sprintf("beta_eu_mean_%s.rds", tag)))
txt <- capture.output(fixest::etable(bmod))
writeLines(txt, file.path("outputs", sprintf("beta_eu_mean_%s.txt", tag)))

# Debug stamp with the exact coefficient & sample we just fit
cat(sprintf("Run: %s | sample=%s | n=%d | beta_baseline=%.6f\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            tag, nrow(dt), coef(bmod)["baseline"]),
    file = file.path("outputs", sprintf("beta_debug_%s.txt", tag)))

# -------- Rolling 5-year β (centered) --------
yrs <- sort(unique(dt$year))
if (length(yrs) >= 7) {
  centers <- yrs[yrs >= (min(yrs)+2) & yrs <= (max(yrs)-2)]
  roll_res <- lapply(centers, function(cy){
    sub <- dt[year >= (cy-2) & year <= (cy+2)]
    if (nrow(sub) < 50) return(data.table(center=cy, beta=NA_real_, se=NA_real_, n=nrow(sub)))
    fit <- tryCatch(feols(g ~ baseline | geo + year, data = sub, cluster = ~geo), error=function(e) NULL)
    if (is.null(fit)) return(data.table(center=cy, beta=NA_real_, se=se, n=nrow(sub)))
    co <- coef(fit)["baseline"]
    se <- sqrt(diag(vcov(fit, cluster = "geo")))["baseline"]
    data.table(center=cy, beta=co, se=se, n=nrow(sub))
  })
  roll_dt <- data.table::rbindlist(roll_res)
  roll_dt[, ci_lo := beta - 1.96*se]
  roll_dt[, ci_hi := beta + 1.96*se]
  data.table::fwrite(roll_dt, file.path("outputs", sprintf("beta_rolling_%s.csv", tag)))

  g_beta <- ggplot(roll_dt, aes(center, beta)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 2020, linetype = 3) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2) +
    geom_line() + geom_point() +
    labs(x = "Center year (5-year window)",
         y = "β (annualized; years/year per 1-year gap)",
         title = sprintf("β-convergence (EU-27; LE at birth; annualized) — %s", tag)) +
    theme_pub()
  ggplot2::ggsave(file.path("figs", sprintf("F1_beta_rolling_%s.png", tag)),
                  g_beta, width = 7, height = 4, dpi = 300)
  # keep legacy filename for continuity
  ggplot2::ggsave("figs/F1_beta_rolling.png", g_beta, width = 7, height = 4, dpi = 300)
} else {
  message("Not enough years for rolling β; skipped rolling plot.")
}

message("Done: σ (log) and β (annualized) exported. Artifacts tagged as: ", tag)
