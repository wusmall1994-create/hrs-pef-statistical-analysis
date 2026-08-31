options(survey.lonely.psu = "adjust")
suppressPackageStartupMessages({
  library(survey)
})

root <- normalizePath(Sys.getenv("PEF_PROJECT_ROOT", unset = "."), mustWork = FALSE)
formal_dir <- file.path(root, "outputs", "formal_analysis")
out_dir <- file.path(root, "outputs", "clinical_value")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

prep_data <- function(cohort) {
  fn <- if (cohort == "HRS") "hrs_formal_dataset.csv" else "charls_formal_dataset.csv"
  d <- read.csv(file.path(formal_dir, fn), check.names = FALSE)
  d$sex_f <- factor(d$sex, levels = c(1, 2), labels = c("Male", "Female"))
  d$smoking_f <- factor(d$smoking, levels = c("never", "former", "current"), labels = c("Never", "Former", "Current"))
  d$pef_resid_z <- -d$pef_lower_1sd  # Higher values indicate higher-than-expected PEF.
  if (cohort == "HRS") {
    d$race_f <- factor(d$race)
    d$wave_f <- factor(d$wave)
  } else {
    d$rural_f <- factor(d$race_context)
    d$education_f <- factor(d$education)
  }
  d
}

covariates <- function(cohort) {
  if (cohort == "HRS") {
    c("age", "I(age^2)", "sex_f", "height_m", "wave_f", "race_f",
      "education_years", "smoking_f", "bmi", "hypertension", "diabetes",
      "heart", "stroke", "cancer", "cesd")
  } else {
    c("age_c10", "age_c10_sq", "sex_f", "height_c10", "rural_f",
      "education_f", "smoking_f", "bmi_c5", "hypertension", "diabetes",
      "heart", "stroke", "cancer", "cesd_per5")
  }
}

make_design <- function(d, cohort) {
  if (cohort == "HRS") {
    svydesign(ids = ~survey_half_sample, strata = ~survey_stratum,
              weights = ~respondent_weight, nest = TRUE, data = d)
  } else {
    svydesign(ids = ~communityID, weights = ~weight, data = d)
  }
}

analytic_data <- function(fit, design) {
  rn <- rownames(model.frame(fit))
  pos <- match(rn, rownames(design$variables))
  if (anyNA(pos)) stop("Could not map model-frame rows to survey-design rows")
  design$variables[pos, , drop = FALSE]
}

weighted_group_risk <- function(design, group_value) {
  des <- subset(design, low_pef_adj == group_value)
  est <- svymean(~event, des, na.rm = TRUE)
  risk <- as.numeric(coef(est)[1]); se <- as.numeric(SE(est)[1])
  c(risk = risk, low = max(0, risk - 1.96 * se), high = min(1, risk + 1.96 * se))
}

standardized_risk <- function(fit, analytic, weight_name, exposure_value) {
  nd <- analytic
  nd$low_pef_adj <- exposure_value
  tt <- delete.response(terms(fit))
  X <- model.matrix(tt, nd, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  beta <- coef(fit); V <- vcov(fit)
  X <- X[, names(beta), drop = FALSE]
  w <- analytic[[weight_name]]
  w <- w / sum(w)
  pred <- as.numeric(exp(X %*% beta))
  risk <- sum(w * pred)
  grad <- colSums(X * as.numeric(w * pred))
  se <- sqrt(as.numeric(t(grad) %*% V %*% grad))
  c(risk = risk, low = max(0, risk - 1.96 * se), high = min(1, risk + 1.96 * se))
}

rcs_basis <- function(x, knots) {
  k <- length(knots)
  if (k < 4) stop("At least four knots are required")
  tp3 <- function(z) pmax(z, 0)^3
  den <- (knots[k] - knots[1])^2
  out <- sapply(seq_len(k - 2), function(j) {
    (tp3(x - knots[j]) -
       ((knots[k] - knots[j]) / (knots[k] - knots[k - 1])) * tp3(x - knots[k - 1]) +
       ((knots[k - 1] - knots[j]) / (knots[k] - knots[k - 1])) * tp3(x - knots[k])) / den
  })
  if (is.null(dim(out))) out <- matrix(out, nrow = length(x), byrow = TRUE)
  colnames(out) <- paste0("pef_rcs", seq_len(ncol(out)))
  out
}

abs_rows <- list(); rd_rows <- list(); spline_tests <- list(); spline_curves <- list(); knot_rows <- list()

for (cohort in c("CHARLS", "HRS")) {
  d <- prep_data(cohort)
  des <- make_design(d, cohort)
  covs <- covariates(cohort)
  wname <- if (cohort == "HRS") "respondent_weight" else "weight"

  primary_formula <- as.formula(paste("event ~", paste(c("low_pef_adj", covs), collapse = " + ")))
  primary_fit <- svyglm(primary_formula, design = des, family = quasipoisson(link = "log"), na.action = na.omit)
  analytic <- analytic_data(primary_fit, des)
  analytic_des <- make_design(analytic, cohort)

  for (g in c(0, 1)) {
    crude <- weighted_group_risk(analytic_des, g)
    adj <- standardized_risk(primary_fit, analytic, wname, g)
    abs_rows[[length(abs_rows) + 1]] <- data.frame(
      cohort = cohort, pef_group = ifelse(g == 1, "Low PEF", "Non-low PEF"),
      n = sum(analytic$low_pef_adj == g), events = sum(analytic$event[analytic$low_pef_adj == g]),
      crude_risk = crude["risk"], crude_ci_low = crude["low"], crude_ci_high = crude["high"],
      adjusted_risk = adj["risk"], adjusted_ci_low = adj["low"], adjusted_ci_high = adj["high"],
      stringsAsFactors = FALSE
    )
  }

  crude_rd_fit <- svyglm(event ~ low_pef_adj, design = analytic_des,
                         family = quasibinomial(link = "identity"))
  crude_b <- coef(crude_rd_fit)["low_pef_adj"]
  crude_se <- sqrt(vcov(crude_rd_fit)["low_pef_adj", "low_pef_adj"])
  a0 <- standardized_risk(primary_fit, analytic, wname, 0)
  a1 <- standardized_risk(primary_fit, analytic, wname, 1)

  # Delta-method covariance for the adjusted marginal risk difference.
  tt <- delete.response(terms(primary_fit)); beta <- coef(primary_fit); V <- vcov(primary_fit)
  nd0 <- analytic; nd0$low_pef_adj <- 0
  nd1 <- analytic; nd1$low_pef_adj <- 1
  X0 <- model.matrix(tt, nd0, contrasts.arg = primary_fit$contrasts, xlev = primary_fit$xlevels)[, names(beta), drop = FALSE]
  X1 <- model.matrix(tt, nd1, contrasts.arg = primary_fit$contrasts, xlev = primary_fit$xlevels)[, names(beta), drop = FALSE]
  ww <- analytic[[wname]] / sum(analytic[[wname]])
  p0 <- as.numeric(exp(X0 %*% beta)); p1 <- as.numeric(exp(X1 %*% beta))
  g0 <- colSums(X0 * as.numeric(ww * p0)); g1 <- colSums(X1 * as.numeric(ww * p1))
  grad_rd <- g1 - g0
  adj_rd <- a1["risk"] - a0["risk"]
  adj_rd_se <- sqrt(as.numeric(t(grad_rd) %*% V %*% grad_rd))
  rr_b <- coef(primary_fit)["low_pef_adj"]
  rr_se <- sqrt(vcov(primary_fit)["low_pef_adj", "low_pef_adj"])

  rd_rows[[length(rd_rows) + 1]] <- data.frame(
    cohort = cohort, n = nrow(analytic), events = sum(analytic$event),
    crude_risk_difference = crude_b,
    crude_rd_ci_low = crude_b - 1.96 * crude_se,
    crude_rd_ci_high = crude_b + 1.96 * crude_se,
    adjusted_risk_difference = adj_rd,
    adjusted_rd_ci_low = adj_rd - 1.96 * adj_rd_se,
    adjusted_rd_ci_high = adj_rd + 1.96 * adj_rd_se,
    adjusted_rr = exp(rr_b), adjusted_rr_ci_low = exp(rr_b - 1.96 * rr_se),
    adjusted_rr_ci_high = exp(rr_b + 1.96 * rr_se),
    stringsAsFactors = FALSE
  )

  # Four-knot restricted cubic spline at survey-weighted 5th, 35th, 65th and 95th percentiles.
  knots <- weighted_quantile(d$pef_resid_z, d[[wname]], c(.05, .35, .65, .95))
  B <- rcs_basis(d$pef_resid_z, knots)
  d$pef_rcs1 <- B[, 1]; d$pef_rcs2 <- B[, 2]
  spline_des <- make_design(d, cohort)
  spline_formula <- as.formula(paste("event ~", paste(c("pef_resid_z", "pef_rcs1", "pef_rcs2", covs), collapse = " + ")))
  spline_fit <- svyglm(spline_formula, design = spline_des, family = quasipoisson(link = "log"), na.action = na.omit)
  spline_analytic <- analytic_data(spline_fit, spline_des)
  overall_test <- regTermTest(spline_fit, ~ pef_resid_z + pef_rcs1 + pef_rcs2)
  nonlinear_test <- regTermTest(spline_fit, ~ pef_rcs1 + pef_rcs2)

  spline_tests[[length(spline_tests) + 1]] <- data.frame(
    cohort = cohort, n = nrow(spline_analytic), events = sum(spline_analytic$event),
    knot_5 = knots[1], knot_35 = knots[2], knot_65 = knots[3], knot_95 = knots[4],
    p_overall = as.numeric(overall_test$p), p_nonlinear = as.numeric(nonlinear_test$p),
    stringsAsFactors = FALSE
  )
  knot_rows[[length(knot_rows) + 1]] <- data.frame(
    cohort = cohort, percentile = c(5, 35, 65, 95), pef_residual_sd = knots,
    stringsAsFactors = FALSE
  )

  bounds <- weighted_quantile(spline_analytic$pef_resid_z, spline_analytic[[wname]], c(.01, .99))
  grid <- seq(bounds[1], bounds[2], length.out = 241)
  Bg <- rcs_basis(grid, knots); Bref <- rcs_basis(0, knots)
  beta_s <- coef(spline_fit); V_s <- vcov(spline_fit)
  term_names <- c("pef_resid_z", "pef_rcs1", "pef_rcs2")
  for (i in seq_along(grid)) {
    delta_small <- c(grid[i], Bg[i, 1], Bg[i, 2]) - c(0, Bref[1, 1], Bref[1, 2])
    delta <- setNames(rep(0, length(beta_s)), names(beta_s))
    delta[term_names] <- delta_small
    est <- sum(delta * beta_s)
    se <- sqrt(as.numeric(t(delta) %*% V_s %*% delta))
    spline_curves[[length(spline_curves) + 1]] <- data.frame(
      cohort = cohort, pef_residual_sd = grid[i], rr = exp(est),
      ci_low = exp(est - 1.96 * se), ci_high = exp(est + 1.96 * se), reference = 0,
      stringsAsFactors = FALSE
    )
  }
}

absolute_risk <- do.call(rbind, abs_rows)
risk_difference <- do.call(rbind, rd_rows)
spline_test <- do.call(rbind, spline_tests)
spline_curve <- do.call(rbind, spline_curves)
spline_knots <- do.call(rbind, knot_rows)

write.csv(absolute_risk, file.path(out_dir, "absolute_risk_results.csv"), row.names = FALSE)
write.csv(risk_difference, file.path(out_dir, "risk_difference_results.csv"), row.names = FALSE)
write.csv(spline_test, file.path(out_dir, "spline_tests.csv"), row.names = FALSE)
write.csv(spline_curve, file.path(out_dir, "spline_curve_source.csv"), row.names = FALSE)
write.csv(spline_knots, file.path(out_dir, "spline_knots.csv"), row.names = FALSE)

cat("Absolute risk results\n"); print(absolute_risk)
cat("\nRisk-difference results\n"); print(risk_difference)
cat("\nSpline tests\n"); print(spline_test)
