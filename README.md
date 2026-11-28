# Peak expiratory flow and subsequent self-reported chronic lung disease

Analysis code for the CHARLS and HRS cohorts. Run Python with pandas and numpy, and R with survey, mice and ggplot2.

Place the authorised Harmonized CHARLS file under `CHARLS date/Harmonized CHARLS/H_CHARLS_D_Data/H_CHARLS_D_Data.dta`, the RAND HRS longitudinal file under `HRS data/01_rand_longitudinal`, the four RAND FAT files under `HRS data/02_rand_fat/{2006,2010,2014,2018}` (one Stata file per year), and `trk2022tr_r.dta` under `HRS data/04_tracker`.

Run `python run_pipeline.py --data-root PATH_TO_DATA --work-dir PRIVATE_OUTPUT`. Use a private working directory for participant-level intermediate files. The pipeline prepares the cohorts, fits survey-weighted association and imputation models, generates held-out predictions and paired design-based intervals, evaluates measurement restrictions and longitudinal PEF patterns, and draws the figures. The results directory contains the model summaries, exposure-history contrasts, reference coefficients, standardised risks and paired intervals. Use `--dry-run` to inspect the commands.

The sensitivity analyses report the CHARLS comparison of persistent low PEF with low PEF at the later assessment only, and the HRS weight comparison in a common physical-measure sample. For the full prediction models, `prediction_bootstrap.R` refits both models and the PEF reference equations during one grouped fivefold validation cycle per bootstrap draw. It writes draw-level paired differences in AUC and Brier score to the private results directory.
