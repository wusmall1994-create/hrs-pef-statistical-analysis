options(survey.lonely.psu = "adjust")
suppressPackageStartupMessages({
  library(survey)
  library(mice)
})

root <- normalizePath(Sys.getenv("PEF_PROJECT_ROOT", unset = "."), mustWork = FALSE)
in_dir <- file.path(root, "outputs", "sensitivity")
formal_dir <- file.path(root, "outputs", "formal_analysis")
dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)

weighted_quantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  x[which(cumsum(w) / sum(w) >= p)[1]]
}

add_exposures <- function(d, cohort) {
  d$pef_residual <- NA_real_
  d$low_pef_adj <- NA_integer_
  d$low_pef_q10 <- NA_integer_
  d$pef_lower_1sd <- NA_real_
  d$pef_residual_w <- NA_real_
  d$low_pef_weighted <- NA_integer_
  d$low_pef_raw_weighted <- NA_integer_
  wname <- if (cohort == "HRS") "respondent_weight" else "weight"
  for (s in sort(unique(d$sex))) {
    ix <- which(d$sex == s & is.finite(d$pef) & is.finite(d$age) & is.finite(d$height_m))
    f <- if (cohort == "HRS") log(pef) ~ age + I(age^2) + height_m + factor(wave) else log(pef) ~ age + I(age^2) + height_m
    fit <- lm(f, data = d[ix, ])
    used <- as.integer(rownames(model.frame(fit)))
    r <- residuals(fit)
    cut20 <- unname(quantile(r, .20, na.rm = TRUE))
    cut10 <- unname(quantile(r, .10, na.rm = TRUE))
    d$pef_residual[used] <- r
    d$low_pef_adj[used] <- as.integer(r <= cut20)
    d$low_pef_q10[used] <- as.integer(r <= cut10)
    d$pef_lower_1sd[used] <- -r / sd(r, na.rm = TRUE)

    fw <- lm(f, data = d[ix, ], weights = d[ix, wname])
    usedw <- as.integer(rownames(model.frame(fw)))
    rw <- residuals(fw)
    ww <- d[usedw, wname]
    cutw <- weighted_quantile(rw, ww, .20)
    rawcut <- weighted_quantile(d$pef[usedw], ww, .20)
    d$pef_residual_w[usedw] <- rw
    d$low_pef_weighted[usedw] <- as.integer(rw <= cutw)
    d$low_pef_raw_weighted[usedw] <- as.integer(d$pef[usedw] <= rawcut)
  }
  d
}

prep <- function(d, cohort) {
  d <- add_exposures(d, cohort)
  d$sex_f <- factor(d$sex, levels = c(1, 2), labels = c("Male", "Female"))
  d$smoking_f <- factor(d$smoking, levels = c("never", "former", "current"), labels = c("Never", "Former", "Current"))
  if (cohort == "HRS") {
    d$race_f <- factor(d$race)
    d$wave_f <- factor(d$wave)
  } else {
    d$rural_f <- factor(d$race_context)
    d$education_f <- factor(d$education)
    d$age_c10 <- (d$age - 60) / 10
    d$age_c10_sq <- d$age_c10^2
    d$height_c10 <- (d$height_m - 1.65) / .10
    d$bmi_c5 <- (d$bmi - 24) / 5
    d$cesd_per5 <- d$cesd / 5
  }
  d
}

covariates <- function(cohort) {
  if (cohort == "HRS") {
    c("age", "I(age^2)", "sex_f", "height_m", "wave_f", "race_f", "education_years", "smoking_f", "bmi",
      "hypertension", "diabetes", "heart", "stroke", "cancer", "cesd")
  } else {
    c("age_c10", "age_c10_sq", "sex_f", "height_c10", "rural_f", "education_f", "smoking_f", "bmi_c5",
      "hypertension", "diabetes", "heart", "stroke", "cancer", "cesd_per5")
  }
}

make_design <- function(d, cohort, weight_name = "analysis_weight") {
  if (cohort == "HRS") {
    svydesign(ids = ~survey_half_sample, strata = ~survey_stratum,
              weights = as.formula(paste0("~", weight_name)), nest = TRUE, data = d)
  } else {
    svydesign(ids = ~communityID, weights = as.formula(paste0("~", weight_name)), data = d)
  }
}

fit_effect <- function(d, cohort, exposure, analysis, outcome = "event", weight_name = NULL) {
  if (is.null(weight_name)) weight_name <- if (cohort == "HRS") "respondent_weight" else "weight"
  d$analysis_weight <- d[[weight_name]]
  des <- make_design(d, cohort, "analysis_weight")
  rhs <- paste(c(exposure, covariates(cohort)), collapse = " + ")
  fit <- svyglm(as.formula(paste(outcome, "~", rhs)), design = des,
                family = quasipoisson(link = "log"), na.action = na.omit)
  b <- coef(fit)[exposure]; se <- sqrt(vcov(fit)[exposure, exposure])
  mf <- model.frame(fit)
  data.frame(cohort = cohort, analysis = analysis, exposure = exposure,
             n = nrow(mf), events = sum(model.response(mf)), rr = exp(b),
             ci_low = exp(b - 1.96 * se), ci_high = exp(b + 1.96 * se),
             p = 2 * pnorm(abs(b / se), lower.tail = FALSE), stringsAsFactors = FALSE)
}

results <- data.frame()
attrition <- data.frame()

for (cohort in c("CHARLS", "HRS")) {
  file <- if (cohort == "HRS") "hrs_first_baseline_all_followup_status.csv" else "charls_baseline_all_followup_status.csv"
  d <- prep(read.csv(file.path(in_dir, file), check.names = FALSE), cohort)
  base_weight <- if (cohort == "HRS") "respondent_weight" else "weight"
  d$analysis_weight <- d[[base_weight]]
  des <- make_design(d, cohort)

  for (g in c("Overall", "Non-low PEF", "Low PEF")) {
    subdes <- if (g == "Overall") des else if (g == "Low PEF") subset(des, low_pef_adj == 1) else subset(des, low_pef_adj == 0)
    for (st in c("observed", "death", "other_nonresponse")) {
      expr <- if (st == "observed") "followup_observed" else if (st == "death") "death_before_followup" else "I(as.numeric(followup_status == 'other_nonresponse'))"
      m <- as.numeric(svymean(as.formula(paste0("~", expr)), subdes, na.rm = TRUE))[1]
      status_n <- if (st == "observed") sum(subdes$variables$followup_observed == 1) else if (st == "death")
        sum(subdes$variables$death_before_followup == 1) else sum(subdes$variables$followup_status == "other_nonresponse")
      attrition <- rbind(attrition, data.frame(cohort = cohort, low_pef_group = g, status = st,
                                               group_n = nrow(subdes$variables), status_n = status_n, weighted_pct = 100 * m))
    }
  }

  # Restrict to a fully observed baseline covariate set for stable response modeling.
  needed <- unique(c("followup_observed", "event", "low_pef_adj", covariates(cohort)))
  plain_needed <- gsub("I\\(|\\)|\\^2", "", needed)
  plain_needed <- plain_needed[plain_needed %in% names(d)]
  cc <- complete.cases(d[, setdiff(plain_needed, "event"), drop = FALSE])
  dcc <- d[cc, ]
  dcc$analysis_weight <- dcc[[base_weight]]
  descc <- make_design(dcc, cohort)
  response_formula <- as.formula(paste("followup_observed ~", paste(c("low_pef_adj", covariates(cohort)), collapse = " + ")))
  rfit <- svyglm(response_formula, design = descc, family = quasibinomial(), na.action = na.omit)
  pobs <- as.numeric(predict(rfit, newdata = dcc, type = "response"))
  pnum <- as.numeric(svymean(~followup_observed, descc, na.rm = TRUE))
  dcc$ipcw <- pnum / pmax(pmin(pobs, .995), .05)
  qs <- quantile(dcc$ipcw, c(.01, .99), na.rm = TRUE)
  dcc$ipcw <- pmin(pmax(dcc$ipcw, qs[1]), qs[2])
  obs <- dcc[dcc$followup_observed == 1 & !is.na(dcc$event), ]
  obs$analysis_weight <- obs[[base_weight]]
  results <- rbind(results, fit_effect(obs, cohort, "low_pef_adj", "earliest_baseline_complete_case", weight_name = "analysis_weight"))
  obs$analysis_weight <- obs[[base_weight]] * obs$ipcw
  results <- rbind(results, fit_effect(obs, cohort, "low_pef_adj", "loss_to_followup_IPCW", weight_name = "analysis_weight"))

  # Exposure-definition sensitivity among the observed baseline cohort.
  observed <- d[d$followup_observed == 1 & !is.na(d$event), ]
  for (x in c("low_pef_weighted", "low_pef_q10", "low_pef_raw_weighted", "pef_lower_1sd")) {
    results <- rbind(results, fit_effect(observed, cohort, x, paste0("alternative_exposure_", x)))
  }
  results <- rbind(results, fit_effect(observed[observed$pef >= 50 & observed$pef <= 800, ], cohort,
                                       "low_pef_adj", "PEF_50_to_800"))

  if (cohort == "CHARLS") {
    results <- rbind(results,
      fit_effect(observed[observed$full_effort == 1, ], cohort, "low_pef_adj", "full_effort_only"),
      fit_effect(observed[observed$repeatable_40 == 1, ], cohort, "low_pef_adj", "repeatable_top2_within_40"),
      fit_effect(observed[observed$repeatable_10pct == 1, ], cohort, "low_pef_adj", "repeatable_top2_within_10pct"),
      fit_effect(observed[is.na(observed$lung_medication) | observed$lung_medication == 0, ], cohort,
                 "low_pef_adj", "exclude_baseline_lung_medication")
    )
  }
}

# Delayed-incidence sensitivity reduces the likelihood that low PEF merely reflects
# an immediately diagnosable but unrecognized baseline condition.
hrs_lag <- prep(read.csv(file.path(in_dir, "hrs_four_year_lagged.csv")), "HRS")
charls_lag <- prep(read.csv(file.path(in_dir, "charls_five_year_lagged.csv")), "CHARLS")
for (x in c("low_pef_adj", "pef_lower_1sd")) {
  results <- rbind(results,
    fit_effect(hrs_lag, "HRS", x, "four_year_delayed_incidence"),
    fit_effect(charls_lag, "CHARLS", x, "five_year_delayed_incidence")
  )
}

# Persistent HRS report: positive at the next interview and again two years later.
hrs_confirm <- read.csv(file.path(formal_dir, "hrs_formal_dataset.csv"), check.names = FALSE)
hrs_confirm$sex_f <- factor(hrs_confirm$sex, levels = c(1, 2), labels = c("Male", "Female"))
hrs_confirm$smoking_f <- factor(hrs_confirm$smoking, levels = c("never", "former", "current"), labels = c("Never", "Former", "Current"))
hrs_confirm$race_f <- factor(hrs_confirm$race); hrs_confirm$wave_f <- factor(hrs_confirm$wave)
hrs_confirm$confirm_observed <- tolower(as.character(hrs_confirm$confirm_observed)) %in% c("true", "1")
hrs_confirm <- hrs_confirm[hrs_confirm$confirm_observed & !is.na(hrs_confirm$event_confirmed), ]
for (x in c("low_pef_adj", "pef_lower_1sd")) {
  results <- rbind(results, fit_effect(hrs_confirm, "HRS", x, "persistent_two_consecutive_reports", outcome = "event_confirmed"))
}

# Multiple-imputation sensitivity for missing M3 covariates.
mi_effect <- function(cohort, exposure, m = 10) {
  file <- if (cohort == "HRS") "hrs_formal_dataset.csv" else "charls_formal_dataset.csv"
  d <- read.csv(file.path(formal_dir, file), check.names = FALSE)
  if (cohort == "HRS") {
    d$sex_f <- factor(d$sex, levels = c(1, 2), labels = c("Male", "Female")); d$wave_f <- factor(d$wave)
    d$race_f <- factor(d$race); d$smoking_f <- factor(d$smoking, levels = c("never", "former", "current"))
    for (v in c("hypertension", "diabetes", "heart", "stroke", "cancer")) d[[v]] <- factor(d[[v]], levels = c(0, 1))
    keep <- c("event", exposure, "age", "sex_f", "height_m", "wave_f", "race_f", "education_years", "smoking_f", "bmi",
              "hypertension", "diabetes", "heart", "stroke", "cancer", "cesd", "respondent_weight", "survey_stratum", "survey_half_sample")
    f <- as.formula(paste("event ~", paste(c(exposure, covariates("HRS")), collapse = " + ")))
    protected <- c("event", exposure, "age", "sex_f", "height_m", "wave_f", "respondent_weight", "survey_stratum", "survey_half_sample")
  } else {
    d$sex_f <- factor(d$sex, levels = c(1, 2), labels = c("Male", "Female")); d$rural_f <- factor(d$race_context)
    d$education_f <- factor(d$education); d$smoking_f <- factor(d$smoking, levels = c("never", "former", "current"))
    d$age_c10 <- (d$age - 60) / 10; d$age_c10_sq <- d$age_c10^2; d$height_c10 <- (d$height_m - 1.65) / .10
    d$bmi_c5 <- (d$bmi - 24) / 5; d$cesd_per5 <- d$cesd / 5
    for (v in c("hypertension", "diabetes", "heart", "stroke", "cancer")) d[[v]] <- factor(d[[v]], levels = c(0, 1))
    keep <- c("event", exposure, "age_c10", "age_c10_sq", "sex_f", "height_c10", "rural_f", "education_f", "smoking_f", "bmi_c5",
              "hypertension", "diabetes", "heart", "stroke", "cancer", "cesd_per5", "weight", "communityID")
    f <- as.formula(paste("event ~", paste(c(exposure, covariates("CHARLS")), collapse = " + ")))
    protected <- c("event", exposure, "age_c10", "age_c10_sq", "sex_f", "height_c10", "weight", "communityID")
  }
  dat <- d[, keep]
  method <- make.method(dat)
  method[intersect(protected, names(method))] <- ""
  pred <- make.predictorMatrix(dat)
  pred[, intersect(c("respondent_weight", "survey_stratum", "survey_half_sample", "weight", "communityID"), colnames(pred))] <- 0
  pred[intersect(protected, rownames(pred)), ] <- 0
  imp <- mice(dat, m = m, maxit = 5, method = method, predictorMatrix = pred, printFlag = FALSE,
              seed = ifelse(cohort == "HRS", 20260829, 20260830))
  q <- u <- numeric(m)
  for (i in seq_len(m)) {
    z <- complete(imp, i)
    des <- if (cohort == "HRS") {
      svydesign(ids = ~survey_half_sample, strata = ~survey_stratum, weights = ~respondent_weight, nest = TRUE, data = z)
    } else svydesign(ids = ~communityID, weights = ~weight, data = z)
    fit <- svyglm(f, design = des, family = quasipoisson(link = "log"))
    q[i] <- coef(fit)[exposure]; u[i] <- vcov(fit)[exposure, exposure]
  }
  qbar <- mean(q); total <- mean(u) + (1 + 1/m) * var(q); se <- sqrt(total)
  data.frame(cohort = cohort, analysis = paste0("multiple_imputation_m", m), exposure = exposure,
             n = nrow(dat), events = sum(dat$event), rr = exp(qbar), ci_low = exp(qbar - 1.96 * se),
             ci_high = exp(qbar + 1.96 * se), p = 2 * pnorm(abs(qbar / se), lower.tail = FALSE))
}

existing_results_file <- file.path(in_dir, "sensitivity_results.csv")
cached_mi <- NULL
if (file.exists(existing_results_file)) {
  old_results <- read.csv(existing_results_file, check.names = FALSE)
  cached_mi <- old_results[grepl("^multiple_imputation", old_results$analysis), ]
}
if (!is.null(cached_mi) && nrow(cached_mi) == 4) {
  results <- rbind(results, cached_mi)
} else {
  for (cohort in c("CHARLS", "HRS")) for (x in c("low_pef_adj", "pef_lower_1sd")) {
    results <- rbind(results, mi_effect(cohort, x, m = 10))
  }
}

missingness <- data.frame()
for (cohort in c("CHARLS", "HRS")) {
  file <- if (cohort == "HRS") "hrs_formal_dataset.csv" else "charls_formal_dataset.csv"
  d <- read.csv(file.path(formal_dir, file), check.names = FALSE)
  vars <- if (cohort == "HRS") c("race", "education_years", "smoking", "bmi", "hypertension", "diabetes", "heart", "stroke", "cancer", "cesd") else
    c("race_context", "education", "smoking", "bmi", "hypertension", "diabetes", "heart", "stroke", "cancer", "cesd")
  for (v in vars) missingness <- rbind(missingness, data.frame(cohort = cohort, variable = v, n = nrow(d),
                                                               missing_n = sum(is.na(d[[v]])), missing_pct = 100 * mean(is.na(d[[v]]))))
}

write.csv(results, file.path(in_dir, "sensitivity_results.csv"), row.names = FALSE)
write.csv(attrition, file.path(in_dir, "attrition_status.csv"), row.names = FALSE)
write.csv(missingness, file.path(in_dir, "covariate_missingness.csv"), row.names = FALSE)
cat("Sensitivity results\n"); print(results)
cat("\nAttrition\n"); print(attrition)
