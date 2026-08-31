# HRS–CHARLS peak expiratory flow statistical analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22207621.svg)](https://doi.org/10.5281/zenodo.22207621)

This repository contains only the statistical analysis code used to examine whether lower-than-expected baseline peak expiratory flow (PEF) preceded subsequent physician diagnosis of chronic lung disease in the Health and Retirement Study (HRS) and the China Health and Retirement Longitudinal Study (CHARLS).

## Scope and disclosure

- No manuscript files are included.
- No participant-level, raw, derived, or aggregate result data are included.
- No figures, tables, reports, logs, credentials, or machine-specific paths are included.
- HRS and CHARLS source files must be obtained independently from their official data custodians and used under the applicable data-use terms.
- All generated files are written under `outputs/`, which is excluded from version control.

## Repository contents

| Script | Purpose |
|---|---|
| `code/01_prepare_hrs.py` | Constructs the HRS analytic dataset and preliminary model outputs from the harmonised HRS public-use file. |
| `code/02_prepare_charls.py` | Constructs the CHARLS analytic dataset and preliminary model outputs from the harmonised CHARLS public-use file. |
| `code/03_primary_analysis.R` | Rebuilds cohort-relative PEF measures and runs survey-weighted primary and subgroup models. |
| `code/04_build_sensitivity_datasets.py` | Constructs datasets for attrition and delayed-outcome sensitivity analyses. |
| `code/05_sensitivity_analysis.R` | Runs prespecified sensitivity analyses, including censoring weights and multiple imputation. |
| `code/06_absolute_risk_spline_analysis.R` | Estimates standardised absolute risks, risk differences, and restricted cubic spline associations. |
| `code/07_create_main_forest.R` | Generates the main cohort and subgroup forest plot. |
| `code/08_create_sensitivity_forest.R` | Generates the sensitivity-analysis forest plot. |
| `code/09_create_clinical_value_figures.R` | Generates absolute-risk and spline figures. |
| `run_pipeline.py` | Runs the analysis scripts in the intended order. |

## Requirements

- Python 3.10 or later
- R 4.3 or later
- Python packages listed in `requirements.txt`
- R packages installed by `install_r_packages.R`

Install dependencies:

```bash
python -m pip install -r requirements.txt
Rscript install_r_packages.R
```

## Input files

The pipeline expects:

1. The RAND HRS longitudinal public-use Stata file containing the variables referenced in `code/01_prepare_hrs.py` and `code/04_build_sensitivity_datasets.py`.
2. The Harmonized CHARLS public-use Stata file containing the variables referenced in `code/02_prepare_charls.py` and `code/04_build_sensitivity_datasets.py`.

The files are not distributed here. Obtain them from the official [HRS data portal](https://hrs.isr.umich.edu/data-products) and [CHARLS data portal](https://charls.charlsdata.com/), subject to registration and data-use requirements.

## Run the full pipeline

From the repository root:

```bash
python run_pipeline.py --hrs-data /path/to/randhrs.dta --charls-data /path/to/harmonized_charls.dta
```

To run analyses without generating publication figures:

```bash
python run_pipeline.py --hrs-data /path/to/randhrs.dta --charls-data /path/to/harmonized_charls.dta --skip-figures
```

Individual scripts may also be run from the repository root after setting `PEF_PROJECT_ROOT`, `HRS_DATA_FILE`, and `CHARLS_DATA_FILE` in the process environment.

## Outputs

Generated participant-level analytic files, statistical results, and figures are written under `outputs/`. This directory is intentionally ignored by Git and must not be committed or uploaded.

## Reproducibility notes

- The cohorts are analysed separately to preserve their survey designs and follow-up structures.
- Low PEF is defined within cohort and sex from the lower tail of residuals after modelling log PEF using age, age squared, and measured height; the HRS model also includes baseline wave.
- Primary associations use survey-weighted modified Poisson regression with design-robust standard errors.
- The code implements relative effects, standardised absolute risks, risk differences, spline analyses, subgroup analyses, inverse-probability-of-censoring weighting, delayed-outcome analyses, alternative exposure definitions, and multiple imputation.

## Licence and citation

The code is released under the MIT License. Version 1.0.0 is archived at [doi:10.5281/zenodo.22207621](https://doi.org/10.5281/zenodo.22207621). Citation metadata are provided in `CITATION.cff` and `.zenodo.json`.
