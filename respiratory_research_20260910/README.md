# Reproducing the Respiratory Research revision

This code produces the upgraded, locked-baseline analyses dated 10 September 2026. The earlier public archive is a different analysis version. No participant-level data or fitted objects containing participant data are distributed here.

## Inputs

Obtain authorised CHARLS and HRS files from their custodians. `build_upgrade.py` expects the following relative structure under the `--data-root` directory (edit the path lookup if your licensed copies have different names):

- `CHARLS date/Harmonized CHARLS/H_CHARLS_D_Data/H_CHARLS_D_Data.dta`
- `HRS data/01_rand_longitudinal/` containing `randhrs1992_2022v1.dta`
- `HRS data/02_rand_fat/2006/`, `2010/`, `2014/`, `2018/`, each containing the corresponding RAND FAT Stata file.

Keep participant-level working directories private. The build script reads only the variables required for the analyses, but the intermediate CSV and RDS files contain individual data and must not be uploaded with the paper.

## Execution

Python requires pandas and numpy. R requires survey, mice and ggplot2 (with their dependencies). The accompanying manuscript supplement records the executed environment. From this directory, with Python and Rscript available:

```powershell
python build_upgrade.py --data-root "/path/to/authorised-data" --out "local_work/data"
Rscript analyze_upgrade.R "local_work/data" "local_work/results"
Rscript predict_upgrade.R "local_work/results"
Rscript make_figures.R "local_work/results" "local_work/figures"
```

The analysis and imputation stage may take substantial time. Prediction uses fixed random seeds and five repeats of five cluster-grouped folds. Seeds, transformations, variables and model formulae are explicit in the scripts.

`validate_upgrade.R` additionally reproduces the original estimates and requires the archived individual-level `charls_formal_dataset.csv` and `hrs_formal_dataset.csv`. These are not redistributed. If authorised local copies are available, run `Rscript validate_upgrade.R "local_work/results" "PATH_TO_ARCHIVED_DATA"`. This validator also writes two private `_validation_local.csv` files for internal checks; do not distribute them.

## Outputs and interpretation

The manuscript supplementary aggregate_results folder contains table and figure sources, model coefficients, reference transformations, per-imputation denominators and validation-fold counts. Do not distribute any `*_local.csv` or `.rds` from local_work: even reference fit objects can contain model frames.

Association references use all eligible locked baselines. Prediction references are fitted separately within each training fold; final coefficient CSVs describe full-data fits for reproducibility. Logistic risk is exp(eta)/(1+exp(eta)), where eta is the coefficient-weighted design row, including the intercept. Factors use R treatment contrasts: smoking never is reference, sex 1 is reference, HRS race 1 and earliest represented wave are reference; CHARLS prediction education uses below primary, primary, middle, high or above, with below primary as reference. For each sex, prediction PEF z equals (log(PEF) minus the matching reference linear predictor) divided by residual_sd. Age is years, height is metres, BMI is kg/m², grip10 is kg/10, and CES-D retains its cohort-specific scale. These are research models, not validated clinical calculators.

The sensitivity exercise for misclassification is deterministic and assumption based. Inadmissible corrected risks are retained and marked. Repeat ranges for validation are not confidence intervals. The analyses are post hoc upgrades; see analysis_amendment.md.

## Version and execution helper

This directory contains the 10 September 2026 analysis revision. The root-level
pipeline and Zenodo version 1.0.0 belong to the earlier analysis. The old DOI does
not archive this revision. Cite the Git commit containing this directory.

Install Python dependencies with `python -m pip install -r requirements.txt`
and R dependencies with `Rscript install_r_packages.R`.
Run `python run_pipeline.py --data-root /path/to/authorised-data --work-dir /path/to/private-work`.
Use `--dry-run` to inspect the commands without reading data or running analyses.
The helper requires all four specified HRS FAT waves for the smoking-intensity
analysis. Arial and Cairo graphics support are used for figure export.

The statistical scripts are unchanged from the locally verified manuscript
analysis. SHA256SUMS.txt identifies those scripts. No data, fitted objects,
results, figures or manuscript files are included in this repository.
