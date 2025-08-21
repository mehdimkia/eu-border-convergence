source("code/00_setup.R")
library(data.table)
le <- readRDS("data/derived/le_panel_raw.rds") |> as.data.table()
le <- le[substr(geo,1,2) %in% EU27 & sex=="T"]; if ("age" %in% names(le)) le <- le[age=="Y_LT1"]
setorder(le, geo, year)
le[, value_lead := shift(value, type="lead"), by=.(geo)]
le[, year_lead  := shift(year,  type="lead"), by=.(geo)]
le[, step := year_lead - year]
le <- le[step >= 1]
le[, g := (value_lead - value) / step]

cat("\n== Level (LE, years) ==\n"); print(summary(le$value))
cat("\n== Annualized growth g (years/year) ==\n")
print(quantile(le$g, c(.01,.05,.5,.95,.99), na.rm=TRUE))  # expect ~ -0.3..+0.3

bmod <- readRDS("outputs/beta_eu_mean_total.rds")
beta_hat <- coef(bmod)["baseline"]
se_beta  <- sqrt(diag(vcov(bmod, cluster="geo")))["baseline"]
eu_mean <- le[, .(eu_mean = mean(value, na.rm=TRUE)), by=.(year)]
le <- merge(le, eu_mean, by="year"); le[, baseline := value - eu_mean]
sd_baseline <- sd(le$baseline, na.rm=TRUE)
cat("\nβ per 1-year:", round(beta_hat,3), "[SE", round(se_beta,3), "]",
    "\nβ per 1-SD:  ", round(beta_hat*sd_baseline,3), "[SE", round(se_beta*sd_baseline,3), "]\n")
