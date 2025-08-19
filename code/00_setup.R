# Setup ----
# install.packages(c('eurostat','data.table','dplyr','fixest','sandwich','lmtest','sf','tmap','clubSandwich'))
pkgs <- c('eurostat','data.table','dplyr','fixest','sandwich','lmtest','sf','tmap','clubSandwich')
invisible(lapply(pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)))
options(stringsAsFactors = FALSE, scipen = 999)
dir.create('data/raw', recursive = TRUE, showWarnings = FALSE)
dir.create('data/derived', recursive = TRUE, showWarnings = FALSE)
dir.create('figs', recursive = TRUE, showWarnings = FALSE)
dir.create('outputs', recursive = TRUE, showWarnings = FALSE)
