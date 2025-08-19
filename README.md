# EU Border vs Non-Border Life Expectancy Convergence (1995–2023)

**Goal.** Update β/σ-convergence across EU-27 NUTS-2 regions through 2023 and contrast border vs non-border trajectories with a 2020 break.

## Repo layout
- `code/` – analysis scripts
- `data/raw/` – untouched Eurostat/GISCO extracts
- `data/derived/` – harmonised panels & flags
- `figs/` – exported figures (F1–F3) and appendix plots
- `manuscript/` – manuscript files
- `outputs/` – tables and logs
- `docs/` – notes, protocols

## Environment
- R >= 4.3 (primary), packages: eurostat, data.table, dplyr, fixest, sandwich, lmtest, sf, tmap, clubSandwich
- Optional Python (for GIS): geopandas, shapely

## Pipeline
1. `00_setup.R` – install/load packages; set paths; figure style.
2. `01_pull_eurostat.R` – download `demo_r_mlifexp` and `tgs00101`; freeze extracts.
3. `02_harmonize_nuts.R` – crosswalk to NUTS-2016; build population weights.
4. `03_construct_flags.R` – border flags from GISCO; EU–non‑EU indicator.
5. `04_sigma_beta.R` – compute σ series; estimate β to EU and national means.
6. `05_event_study.R` – border×year event‑study; pre‑trend checks.
7. `06_robustness.R` – LOCO; alternative σ; spatial diagnostics.

## Conventions
- Units: NUTS-2 (NUTS-2016), EU-27.
- Outcomes: LE at birth (total & sex-specific); optional LE@65 for robustness.
- Weights: population in year t (primary).
- Figures: F1 β; F2 σ; F3 event-study. Keep causal language out.

(c) 2025-08-19 — Draft skeleton.
