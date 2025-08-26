# zz_qc_event_study.R — QA for 05_event_study.R outputs
suppressPackageStartupMessages({
  library(data.table); library(readr); library(fixest)
})

# ---- Params: match the toggles you used in 05_event_study.R ----
WEIGHTED <- TRUE
USE_ANCHOR_WEIGHTS <- TRUE
ANCHOR_YEAR <- 2019

# Filenames mirror 05_event_study.R
suffix_w     <- if (WEIGHTED) "_Tonly_wpop" else "_Tonly_unw"
weight_label <- if (WEIGHTED && USE_ANCHOR_WEIGHTS) sprintf("_%d", ANCHOR_YEAR) else ""

dir.create("outputs/QA", recursive = TRUE, showWarnings = FALSE)

# ---- Paths ----
rds_path   <- paste0("outputs/eventstudy_total", suffix_w, weight_label, ".rds")
txt_path   <- paste0("outputs/eventstudy_total", suffix_w, weight_label, ".txt")
coef_csv   <- paste0("outputs/eventstudy_coefs",  suffix_w, weight_label, ".csv")
wald_csv   <- paste0("outputs/eventstudy_pretrend_wald", suffix_w, weight_label, ".csv")
fig_path   <- paste0("figs/F3_event_study", suffix_w, weight_label, ".png")

# ---- 1) Model exists & readable ----
if (!file.exists(rds_path)) stop("Missing model RDS: ", rds_path)
m <- readRDS(rds_path)
cat(sprintf("Model loaded: nobs=%d | fe: %s\n", nobs(m), paste0(m$fixef_vars, collapse = ", ")))

# ---- 2) Coefs CSV exists & matches model ----
if (!file.exists(coef_csv)) stop("Missing coef CSV: ", coef_csv)
co <- data.table::fread(coef_csv)
# pull coefs from model via broom and compare by term
bm <- as.data.table(broom::tidy(m))
bm <- bm[grepl("^post_year::", term)]
co <- co[grepl("^post_year::", term)]

# basic checks
same_terms <- setequal(co$term, bm$term)
cat("Coef CSV terms match model terms: ", same_terms, "\n")

# align and compare estimates within small tolerance
setkey(co, term); setkey(bm, term)
merged <- bm[co, nomatch = 0]
diff_est <- max(abs(merged$estimate - merged$i.estimate), na.rm = TRUE)
cat(sprintf("Max |estimate(model)-estimate(csv)|: %.3e\n", diff_est))

# CI sanity
bad_ci <- co[!(is.finite(ci_lo) & is.finite(ci_hi)) | ci_hi <= ci_lo]
fwrite(bad_ci, "outputs/QA/event_ci_issues.csv")
cat(sprintf("CI sanity — bad rows: %d\n", nrow(bad_ci)))

# Pre/post coverage around 2019
co[, year_rel := as.integer(gsub("post_year::(-?\\d+).*", "\\1", term))]
has_pre  <- any(co$year_rel < 0)
has_post <- any(co$year_rel > 0)
rng_pre  <- if (has_pre) paste0("[", min(co$year_rel[co$year_rel<0]), ",", max(co$year_rel[co$year_rel<0]), "]") else "none"
rng_post <- if (has_post) paste0("[", min(co$year_rel[co$year_rel>0]), ",", max(co$year_rel[co$year_rel>0]), "]") else "none"
cat(sprintf("Coverage — pre: %s | post: %s\n", rng_pre, rng_post))

# ---- 3) Pre-trend Wald CSV exists & p-value sane ----
if (!file.exists(wald_csv)) {
  warning("Missing pre-trend Wald CSV: ", wald_csv)
} else {
  w <- fread(wald_csv)
  ok_p <- is.finite(w$p_value) && w$p_value >= 0 && w$p_value <= 1
  cat(sprintf("Pre-trend Wald — rows=%d | p in [0,1]: %s | p=%.4g\n",
              nrow(w), ok_p, if (nrow(w)) w$p_value[1] else NA_real_))
}

# ---- 4) Figure exists ----
cat(sprintf("Figure exists: %s\n", file.exists(fig_path)))
