# Data notes

- Place Eurostat extracts in `data/raw/eurostat/`. The scripts will create directories as needed.
- GISCO shapes are cached by giscoR; we also save a frozen copy to `data/raw/gisco/` for reproducibility.
- If a NUTS crosswalk is needed (old -> NUTS-2016), place a CSV at `data/raw/crosswalks/nuts_crosswalk.csv` with columns:
  `old_code,new_code,share,year_from,year_to`.
