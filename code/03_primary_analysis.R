options(survey.lonely.psu = "adjust")
suppressPackageStartupMessages(library(survey))

root <- normalizePath(Sys.getenv("PEF_PROJECT_ROOT", unset = "."), mustWork = FALSE)
hrs_path <- file.path(root, "outputs", "hrs_preparation", "hrs_first_eligible_analysis_dataset.csv")
charls_path <- file.path(root, "outputs", "charls_preparation", "charls_2015_2018_harmonized_dataset.csv")
out_dir <- file.path(root, "outputs", "formal_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rebuild_pef <- function(d, include_wave = FALSE) {
  d$pef_residual <- NA_real_
  d$low_pef_adj <- NA_integer_
  d$pef_lower_1sd <- NA_real_
  for (s in sort(unique(d$sex))) {
    ix <- which(d$sex == s)
    if (include_wave) {
      fit <- lm(log(pef) ~ age + I(age^2) + height_m + factor(wave), data = d[ix, ])
    } else {
      fit <- lm(log(pef) ~ age + I(age^2) + height_m, data = d[ix, ])
    }
    rr <- residuals(fit)
    cut <- unname(quantile(rr, 0.20, na.rm = TRUE))
    sdv <- sd(rr, na.rm = TRUE)
    used <- ix[as.integer(names(rr))]
    # lm residual names preserve the original row names, which are sequential
    # after filtering below; direct model-frame row names are safest.
    used <- as.integer(rownames(model.frame(fit)))
    d$pef_residual[used] <- rr
    d$low_pef_adj[used] <- as.integer(rr <= cut)
    d$pef_lower_1sd[used] <- -rr / sdv
  }
  d
}

hrs <- read.csv(hrs_path, check.names = FALSE)
hrs <- hrs[hrs$age >= 50 & hrs$age <= 80 & hrs$respondent_weight > 0, ]
rownames(hrs) <- seq_len(nrow(hrs))
hrs <- rebuild_pef(hrs, include_wave = TRUE)
hrs$sex_f <- factor(hrs$sex, levels = c(1, 2), labels = c("Male", "Female"))
hrs$age_group <- factor(ifelse(hrs$age < 65, "50-64", "65-80"), levels = c("50-64", "65-80"))
hrs$smoking_f <- factor(hrs$smoking, levels = c("never", "former", "current"), labels = c("Never", "Former", "Current"))
hrs$race_f <- factor(hrs$race)
hrs$wave_f <- factor(hrs$wave)
hrs$depression_high <- as.integer(hrs$cesd >= 4)

charls <- read.csv(charls_path, check.names = FALSE)
charls <- charls[charls$age >= 50 & charls$age <= 80 & charls$weight > 0, ]
rownames(charls) <- seq_len(nrow(charls))
charls <- rebuild_pef(charls, include_wave = FALSE)
charls$sex_f <- factor(charls$sex, levels = c(1, 2), labels = c("Male", "Female"))
charls$age_group <- factor(ifelse(charls$age < 65, "50-64", "65-80"), levels = c("50-64", "65-80"))
charls$smoking_f <- factor(charls$smoking, levels = c("never", "former", "current"), labels = c("Never", "Former", "Current"))
charls$rural_f <- factor(charls$race_context)
charls$education_f <- factor(charls$education)
charls$depression_high <- as.integer(charls$cesd >= 10)

hrs_design <- svydesign(
  ids = ~survey_half_sample, strata = ~survey_stratum,
  weights = ~respondent_weight, nest = TRUE, data = hrs
)
charls_design <- svydesign(ids = ~communityID, weights = ~weight, data = charls)

hrs_cov <- c(
  "age", "I(age^2)", "sex_f", "height_m", "wave_f", "race_f",
  "education_years", "smoking_f", "bmi", "hypertension", "diabetes",
  "heart", "stroke", "cancer", "cesd"
)
charls_cov <- c(
  "age_c10", "age_c10_sq", "sex_f", "height_c10", "rural_f",
  "education_f", "smoking_f", "bmi_c5", "hypertension", "diabetes",
  "heart", "stroke", "cancer", "cesd_per5"
)

fit_one <- function(design, exposure, covariates, cohort, analysis, model) {
  rhs <- paste(c(exposure, covariates), collapse = " + ")
  f <- as.formula(paste("event ~", rhs))
  fit <- svyglm(f, design = design, family = quasipoisson(link = "log"), na.action = na.omit)
  b <- coef(fit)[exposure]
  se <- sqrt(vcov(fit)[exposure, exposure])
  mf <- model.frame(fit)
  data.frame(
    cohort = cohort, analysis = analysis, exposure = exposure, model = model,
    n = nrow(mf), events = sum(model.response(mf)),
    rr = exp(b), ci_low = exp(b - 1.96 * se), ci_high = exp(b + 1.96 * se),
    p = 2 * pnorm(abs(b / se), lower.tail = FALSE), stringsAsFactors = FALSE
  )
}

main_results <- data.frame()
for (exposure in c("low_pef_adj", "pef_lower_1sd")) {
  main_results <- rbind(
    main_results,
    fit_one(hrs_design, exposure, character(0), "HRS", "survey_weighted", "M0"),
    fit_one(hrs_design, exposure, hrs_cov[1:6], "HRS", "survey_weighted", "M1"),
    fit_one(hrs_design, exposure, hrs_cov[1:10], "HRS", "survey_weighted", "M2"),
    fit_one(hrs_design, exposure, hrs_cov, "HRS", "survey_weighted", "M3"),
    fit_one(charls_design, exposure, character(0), "CHARLS", "survey_weighted", "M0"),
    fit_one(charls_design, exposure, charls_cov[1:5], "CHARLS", "survey_weighted", "M1"),
    fit_one(charls_design, exposure, charls_cov[1:9], "CHARLS", "survey_weighted", "M2"),
    fit_one(charls_design, exposure, charls_cov, "CHARLS", "survey_weighted", "M3")
  )
}

subgroup_fit <- function(design, data, cohort, exposure, modifier, covariates) {
  levels_mod <- levels(data[[modifier]])
  out <- data.frame()
  remove_terms <- switch(
    modifier,
    sex_f = "sex_f",
    smoking_f = "smoking_f",
    age_group = character(0),
    character(0)
  )
  cov_sub <- setdiff(covariates, remove_terms)

  int_formula <- as.formula(paste(
    "event ~", paste(c(paste0(exposure, " * ", modifier), cov_sub), collapse = " + ")
  ))
  int_fit <- svyglm(int_formula, design = design, family = quasipoisson(link = "log"), na.action = na.omit)
  int_test <- tryCatch(
    regTermTest(int_fit, as.formula(paste0("~", exposure, ":", modifier))),
    error = function(e) NULL
  )
  int_p <- if (is.null(int_test)) NA_real_ else as.numeric(int_test$p)

  for (lev in levels_mod) {
    sub_design <- subset(design, get(modifier) == lev)
    fit <- svyglm(
      as.formula(paste("event ~", paste(c(exposure, cov_sub), collapse = " + "))),
      design = sub_design, family = quasipoisson(link = "log"), na.action = na.omit
    )
    b <- coef(fit)[exposure]
    se <- sqrt(vcov(fit)[exposure, exposure])
    mf <- model.frame(fit)
    out <- rbind(out, data.frame(
      cohort = cohort, exposure = exposure, modifier = modifier, level = lev,
      n = nrow(mf), events = sum(model.response(mf)), rr = exp(b),
      ci_low = exp(b - 1.96 * se), ci_high = exp(b + 1.96 * se),
      p = 2 * pnorm(abs(b / se), lower.tail = FALSE), interaction_p = int_p,
      stringsAsFactors = FALSE
    ))
  }
  out
}

subgroup_results <- data.frame()
for (exposure in c("low_pef_adj", "pef_lower_1sd")) {
  for (modifier in c("sex_f", "age_group", "smoking_f")) {
    subgroup_results <- rbind(
      subgroup_results,
      subgroup_fit(hrs_design, hrs, "HRS", exposure, modifier, hrs_cov),
      subgroup_fit(charls_design, charls, "CHARLS", exposure, modifier, charls_cov)
    )
  }
}

fmt_mean_sd <- function(des, var) {
  m <- as.numeric(svymean(as.formula(paste0("~", var)), des, na.rm = TRUE))
  s <- sqrt(as.numeric(svyvar(as.formula(paste0("~", var)), des, na.rm = TRUE)))
  sprintf("%.1f (%.1f)", m, s)
}

fmt_prop <- function(des, expression_text) {
  # Coerce the logical expression to 0/1.  Without this coercion svymean()
  # treats it as a two-level factor and the first coefficient is P(FALSE).
  m <- as.numeric(svymean(as.formula(paste0("~I(as.numeric(", expression_text, "))")), des, na.rm = TRUE))
  sprintf("%.1f%%", 100 * m)
}

table1_one <- function(design, cohort) {
  complete_design <- subset(design, !is.na(low_pef_adj))
  groups <- list(Overall = complete_design, `Non-low PEF` = subset(complete_design, low_pef_adj == 0), `Low PEF` = subset(complete_design, low_pef_adj == 1))
  rows <- list(
    list("Unweighted N", "n", NULL),
    list("Age, years", "mean", "age"),
    list("Female", "prop", "sex == 2"),
    list("Measured height, m", "mean", "height_m"),
    list("PEF, L/min", "mean", "pef"),
    list("BMI, kg/m2", "mean", "bmi"),
    list("Never smoker", "prop", "smoking_f == 'Never'"),
    list("Former smoker", "prop", "smoking_f == 'Former'"),
    list("Current smoker", "prop", "smoking_f == 'Current'"),
    list("Hypertension", "prop", "hypertension == 1"),
    list("Diabetes", "prop", "diabetes == 1"),
    list("Heart disease", "prop", "heart == 1"),
    list("Stroke", "prop", "stroke == 1"),
    list("Cancer", "prop", "cancer == 1"),
    list("Elevated depressive symptoms", "prop", "depression_high == 1")
  )
  out <- data.frame(Cohort = cohort, Characteristic = vapply(rows, `[[`, "", 1), stringsAsFactors = FALSE)
  for (nm in names(groups)) {
    des <- groups[[nm]]
    vals <- character(length(rows))
    for (i in seq_along(rows)) {
      typ <- rows[[i]][[2]]
      val <- rows[[i]][[3]]
      vals[i] <- if (typ == "n") {
        as.character(nrow(des$variables))
      } else if (typ == "mean") {
        fmt_mean_sd(des, val)
      } else {
        fmt_prop(des, val)
      }
    }
    out[[nm]] <- vals
  }
  out
}

table1 <- rbind(table1_one(charls_design, "CHARLS"), table1_one(hrs_design, "HRS"))

write.csv(main_results, file.path(out_dir, "main_results.csv"), row.names = FALSE)
write.csv(subgroup_results, file.path(out_dir, "subgroup_results.csv"), row.names = FALSE)
write.csv(table1, file.path(out_dir, "table1.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(charls, file.path(out_dir, "charls_formal_dataset.csv"), row.names = FALSE)
write.csv(hrs, file.path(out_dir, "hrs_formal_dataset.csv"), row.names = FALSE)

cat("Main results\n")
print(main_results)
cat("\nSubgroup results\n")
print(subgroup_results)
