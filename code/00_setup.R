# 00_setup.R  ───────────────────────────────────────────────────────────────────

# Packages
pkgs <- c(
  "data.table","dplyr","stringr","readr","lubridate",
  "eurostat","giscoR","sf","fixest","sandwich","lmtest",
  "ggplot2","jsonlite","broom"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "),
    "\nInstall them once via install.packages(c(...)) and re-run.",
    call. = FALSE
  )
}

# Attach the few you’ll likely use unqualified
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(fixest)
})

# Options
options(
  stringsAsFactors   = FALSE,
  scipen             = 999,
  eurostat_cache     = TRUE,
  eurostat_cache_dir = "data/raw/eurostat"
)

# Folders
dirs <- c(
  "data/raw/eurostat","data/raw/gisco","data/raw/crosswalks",
  "data/derived","figs","outputs"
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Wire caches to the packages
giscoR::gisco_set_cache_dir("data/raw/gisco")

# Constants
EU27 <- c("AT","BE","BG","HR","CY","CZ","DK","EE","FI","FR","DE",
          "GR","HU","IE","IT","LV","LT","LU","MT","NL","PL","PT",
          "RO","SK","SI","ES","SE")

# Theme
theme_pub <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      plot.title.position= "plot",
      plot.caption       = ggplot2::element_text(hjust = 0)
    )
}
