# Harmonise to NUTS-2016 ----
source('code/00_setup.R'); library(data.table); library(readr)

dem <- readRDS('data/raw/eurostat/demo_r_mlifexp_NUTS2_1995_2023.rds')

# If any older NUTS vintage appears, provide crosswalk at data/raw/crosswalks/nuts_crosswalk.csv
# cols: old_code,new_code,share,year_from,year_to
if (file.exists('data/raw/crosswalks/nuts_crosswalk.csv')){
  cw <- data.table(readr::read_csv('data/raw/crosswalks/nuts_crosswalk.csv', show_col_types = FALSE))
  dem <- merge(dem, cw, by.x='geo', by.y='old_code', all.x=TRUE)
  dem[, geo_2016 := ifelse(is.na(new_code), geo, new_code)]
  dem[, value := ifelse(is.na(share), value, value*share)]
  dem[, geo := geo_2016][, c('new_code','geo_2016','share'):=NULL]
}

saveRDS(dem, 'data/derived/le_panel_raw.rds')
data.table::fwrite(dem, 'data/derived/le_panel_raw.csv')

message('Saved harmonised panel to data/derived/.')
