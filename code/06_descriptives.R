# code/06_descriptives.R  — coverage stats + EU means for §3.1

suppressPackageStartupMessages({
  library(data.table); library(readr)
})

# Load derived panels + flags
le    <- readRDS("data/derived/le_panel_raw.rds")   |> as.data.table()
pop   <- readRDS("data/derived/pop_panel_raw.rds")  |> as.data.table()
flags <- fread("data/derived/border_flags.csv")     |> as.data.table()

# Normalise flag names to lower-case (as produced by 03_border_flags.R)
setnames(flags, tolower(names(flags)))

# Merge flags (use whatever columns exist; 'nuts_id' is the join key)
flag_cols <- intersect(c("nuts_id","eu_member","border","eu_noneu_touch","cntr_code"), names(flags))
le <- merge(le, flags[, ..flag_cols], by.x = "geo", by.y = "nuts_id", all.x = TRUE)

# Filter to EU-27, LE at birth, Total sex
le <- le[eu_member == TRUE]
if ("age" %in% names(le)) le <- le[age %in% c("Y_LT1","Y0")]
le <- le[sex == "T"]

# Universe & coverage
years     <- sort(unique(le$year))
n_years   <- length(years)
regions   <- sort(unique(le$geo))
N_regions <- length(regions)
N_obs     <- nrow(le)
balanced_max  <- N_regions * n_years
coverage_pct  <- 100 * N_obs / balanced_max

# Border counts (time-invariant)
by_geo       <- unique(le[, .(geo, border, eu_noneu_touch)])
N_border     <- by_geo[isTRUE(border), .N]
N_nonborder  <- by_geo[!isTRUE(border), .N]
P_border     <- 100 * N_border / (N_border + N_nonborder)
P_nonborder  <- 100 - P_border
N_ext_touch  <- by_geo[isTRUE(eu_noneu_touch), .N]

# Anchor-2019 population weights for EU means (nearest available year if 2019 missing)
ANCHOR_YEAR <- 2019
pa <- copy(pop)[!is.na(pop)]
pa[, dist := abs(year - ANCHOR_YEAR)]
setorder(pa, geo, dist, year)
pa <- pa[, .SD[1L], by = geo][, .(geo, w_pop = as.numeric(pop))]

le <- merge(le, pa, by = "geo", all.x = TRUE)

# Population-weighted EU mean LE by year
eu_mean <- le[, .(eu_mean_le = weighted.mean(as.numeric(value), w = as.numeric(w_pop), na.rm = TRUE)),
              by = year]

# Pull the few headline numbers
grab <- function(y) eu_mean[year == y, eu_mean_le][1]
LE_1995 <- grab(1995)
LE_2019 <- grab(2019)
LE_2020 <- grab(2020)
LE_2023 <- grab(2023)

# Round & write summary
sumDT <- data.table(
  N_regions   = N_regions,
  N_obs       = N_obs,
  coverage_pct= round(coverage_pct, 1),
  N_border    = N_border,
  P_border    = round(P_border, 1),
  N_nonborder = N_nonborder,
  P_nonborder = round(P_nonborder, 1),
  N_ext_touch = N_ext_touch,
  LE_1995     = round(LE_1995, 2),
  LE_2019     = round(LE_2019, 2),
  LE_2020     = round(LE_2020, 2),
  LE_2023     = round(LE_2023, 2)
)

fwrite(sumDT, "outputs/descriptives_summary.csv")
print(sumDT)
