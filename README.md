# Code for: Joint Developmental Trajectories of Cognitive Function and PM2.5 Exposure Among Middle-Aged and Older Adults in China

Reproducible analysis code for the manuscript. The pipeline links four waves of the China Health and Retirement Longitudinal Study (CHARLS; 2011, 2013, 2015, 2018) with city-level annual PM2.5 estimates from the CHAP dataset, fits group-based multi-trajectory models (GBMTM) jointly on cognitive score and PM2.5 exposure, and profiles the resulting trajectory groups with multinomial logistic regression.

## Data

Participant-level data are not distributed in this repository.

- CHARLS survey data: https://charls.pku.edu.cn (registration and data-use agreement required; participant-level data cannot be redistributed).
- CHAP PM2.5 dataset: https://data.tpdc.ac.cn (Wei et al. 2021, Remote Sensing of Environment 252, 112136).
- Province boundaries for Figure 3 (download to `data/china_provinces_datav.json`): https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json

The Python figure script (`python/fig_map_5class.py`) reads `tables/` and `data/` relative to the repository root; run it from the repo root after the R pipeline has written the region composition table.

Place the required CHARLS files under `data/charls/` following the layout expected by `R/config.R` (see the `PATHS` list), and the merged city-year PM2.5 file at `data/pollution.xlsx`. Environment variables `CHARLS_DIR` and `POLLUTION_FILE` can override these locations.

## Requirements

- R (4.5+) with packages: `gbmt`, `doParallel`, `foreach`, `dplyr`, `mice`, `nnet`, `haven`, `ggplot2`
- Python 3 with `pandas`, `geopandas`, `matplotlib` (regional map figure only)

## Pipeline

Run order (scripts in `R/`):

| Script | Output | Manuscript element |
| --- | --- | --- |
| `run_stage1.R` | `trajectory_long.csv` | cohort construction (age >= 45, >= 3 valid waves on both indicators) |
| `run_stage2.R` | `model_comparison.csv`, `gbmtm_ng{2..6}.rds`, `final_trajectory_assignment.csv`, `group_reordering.csv` | Table 1; group labels G1-G5 ordered by z(mean PM2.5) - z(cognitive change) |
| `run_descriptives.R` | sample and wave trend tables | Section 3.1 |
| `run_group_profile.R` | group profile table | Section 3.3 (class sizes, mean PM2.5, fitted trajectories) |
| `run_regional.R` | region composition table | Section 3.3 (four economic regions) |
| `run_stage3a_baseline.R` / `run_stage3a_add_region.R` | baseline table with Bonferroni-corrected tests | Table 2 |
| `run_build_baseline.R` | `baseline_5class.csv` (group assignment joined onto baseline covariates; IDs zero-padded to 12 digits) | input for Sections 3.1, 3.3-3.6 |
| `run_stage3c.R` | missingness assessment, MICE (m = 10), imputation density figure | Section 3.4, Figure S1 |
| `run_stage3d.R` | multinomial logistic models 1-3, forest plot | Section 3.5, Figure 4, Table S2 |
| `run_stage3e.R` | sensitivity analysis in the high-confidence subsample (posterior probability >= 0.70) | Section 3.6 |
| `python/fig_map_5class.py` | regional choropleth figure | Figure 3 |

A single sweep is reproduced by running the scripts in the order above with `Rscript`; shared settings (seed = 2026, reference group, variable lists) live in `R/config.R`.

## Notes

- Fixed seed (2026) is set in `R/config.R`; EM starting values follow the defaults of the `gbmt` package.
- Model selection criteria: BIC comparison plus minimum classification-quality thresholds (entropy >= 0.70, APPA >= 0.70, smallest class >= 10%), implemented in `run_stage2.R`.
- The group-ordering rule and the region mapping of provinces follow the definitions described in the Methods section of the manuscript.
