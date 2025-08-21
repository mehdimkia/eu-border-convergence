# Setup ----
pkgs <- c('data.table','dplyr','stringr','readr','lubridate','eurostat','giscoR','sf','fixest','sandwich','lmtest','ggplot2','jsonlite','broom')
invisible(lapply(pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)))
options(stringsAsFactors = FALSE, scipen = 999, eurostat_cache = TRUE)

dir.create('data/raw/eurostat', recursive = TRUE, showWarnings = FALSE)
dir.create('data/raw/gisco', recursive = TRUE, showWarnings = FALSE)
dir.create('data/raw/crosswalks', recursive = TRUE, showWarnings = FALSE)
dir.create('data/derived', recursive = TRUE, showWarnings = FALSE)
dir.create('figs', recursive = TRUE, showWarnings = FALSE)
dir.create('outputs', recursive = TRUE, showWarnings = FALSE)

EU27 <- c('AT','BE','BG','HR','CY','CZ','DK','EE','FI','FR','DE','GR','HU','IE','IT','LV','LT','LU','MT','NL','PL','PT','RO','SK','SI','ES','SE')

theme_pub <- function(){
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = element_blank(),
                   plot.title.position = "plot",
                   plot.caption = ggplot2::element_text(hjust = 0))
}
