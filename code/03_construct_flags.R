# Border flags via GISCO ----
source('code/00_setup.R'); library(sf); library(dplyr); library(data.table); library(giscoR)

nuts <- giscoR::gisco_get_nuts(year = 2016, resolution = "60")
nuts2 <- nuts %>% dplyr::filter(LEVL_CODE == 2) %>%
  dplyr::select(NUTS_ID, CNTR_CODE, geometry)

tch <- sf::st_touches(nuts2)
border_flag <- logical(nrow(nuts2))

for (i in seq_len(nrow(nuts2))){
  nbrs <- unlist(tch[[i]])
  if (length(nbrs) == 0) { border_flag[i] <- FALSE; next }
  border_flag[i] <- any(nuts2$CNTR_CODE[nbrs] != nuts2$CNTR_CODE[i])
}

nuts2$border <- border_flag
nuts2$eu_member <- nuts2$CNTR_CODE %in% EU27

eu_noneu_touch <- logical(nrow(nuts2))
for (i in seq_len(nrow(nuts2))){
  nbrs <- unlist(tch[[i]])
  if (length(nbrs) == 0) { eu_noneu_touch[i] <- FALSE; next }
  eu_noneu_touch[i] <- any( (nuts2$CNTR_CODE[nbrs] %in% EU27) != (nuts2$CNTR_CODE[i] %in% EU27) )
}
nuts2$eu_noneu_touch <- eu_noneu_touch

flags <- nuts2 %>% sf::st_drop_geometry() %>% as.data.table()
data.table::fwrite(flags, 'data/derived/border_flags.csv')
sf::st_write(nuts2, 'data/derived/nuts2_2016_flags.gpkg', delete_dsn = TRUE)

message('Saved border flags to data/derived/.')
