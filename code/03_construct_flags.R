# Border flags via GISCO ----
source("code/00_setup.R")
suppressPackageStartupMessages({ library(sf); library(dplyr); library(data.table); library(giscoR) })

# Get NUTS-2 (2016) and clean
nuts <- giscoR::gisco_get_nuts(year = 2016, resolution = "20")
nuts2 <- nuts |>
  dplyr::filter(LEVL_CODE == 2) |>
  dplyr::select(NUTS_ID, CNTR_CODE, geometry) |>
  sf::st_make_valid()

# EU membership (GISCO uses EL for Greece)
EU27_gisco <- union(EU27, "EL")
nuts2$eu_member <- nuts2$CNTR_CODE %in% EU27_gisco

# Use GEOS for reliable touches
old_s2 <- sf::sf_use_s2()
sf::sf_use_s2(FALSE)
tch <- sf::st_touches(nuts2)
sf::sf_use_s2(old_s2)

# Flags
idx <- seq_len(nrow(nuts2))
border_flag <- vapply(idx, function(i) {
  nbr <- unlist(tch[[i]]); if (!length(nbr)) return(FALSE)
  any(nuts2$CNTR_CODE[nbr] != nuts2$CNTR_CODE[i])
}, logical(1))

eu_noneu_touch <- vapply(idx, function(i) {
  nbr <- unlist(tch[[i]]); if (!length(nbr)) return(FALSE)
  any(nuts2$eu_member[nbr] != nuts2$eu_member[i])
}, logical(1))

nuts2$border <- border_flag
nuts2$eu_noneu_touch <- eu_noneu_touch

# Save neighbor edges (useful QA)
edges <- data.table::rbindlist(lapply(idx, function(i) {
  nbr <- unlist(tch[[i]]); if (!length(nbr)) return(NULL)
  data.table::data.table(src = nuts2$NUTS_ID[i], dst = nuts2$NUTS_ID[nbr])
}))
data.table::fwrite(edges, "data/derived/nuts2_touches_edges.csv")

# Outputs
flags <- nuts2 |> sf::st_drop_geometry() |> as.data.table()
data.table::fwrite(flags, "data/derived/border_flags.csv")
sf::st_write(nuts2, "data/derived/nuts2_2016_flags.gpkg", delete_dsn = TRUE)

message("Saved border flags and GPKG to data/derived/.")
