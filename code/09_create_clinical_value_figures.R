suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

root <- normalizePath(Sys.getenv("PEF_PROJECT_ROOT", unset = "."), mustWork = FALSE)
in_dir <- file.path(root, "outputs", "clinical_value")
fig_dir <- file.path(in_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c(navy = "#17365D", blue = "#4E79A7", orange = "#D55E00",
             light_blue = "#D9EAF7", light_orange = "#FCE4D6", grey = "#6B7280")

theme_pub <- function(base_size = 8.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      axis.title = element_text(colour = "#222222"),
      strip.background = element_rect(fill = "#EAF0F7", colour = NA),
      strip.text = element_text(face = "bold", size = base_size + 0.2),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1.2, margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "#4B5563", margin = margin(b = 6)),
      plot.caption = element_text(size = base_size - 1.0, colour = "#555555", hjust = 0, margin = margin(t = 6)),
      legend.title = element_blank(),
      legend.position = "top"
    )
}

save_pub <- function(plot, stem, width_mm = 183, height_mm = 105, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = "Arial")
  print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, "_600dpi.tiff"), width = w, height = h, units = "in", res = dpi)
  print(plot); dev.off()
  ragg::agg_png(paste0(stem, "_preview.png"), width = w, height = h, units = "in", res = 180)
  print(plot); dev.off()
}

# Figure 1: parallel design and harmonization.
box_df <- data.frame(
  cohort = rep(c("CHARLS", "HRS"), each = 3),
  x = rep(c(1.7, 5.0, 8.3), 2),
  y = rep(c(2.7, 1.25), each = 3),
  label = c(
    "Wave 3 (2015)\nAge 50-80; PEF measured\nNo reported chronic lung disease",
    "Wave 4 (2018)\nFirst reported physician\ndiagnosis; ~3-year interval",
    "N=8,087; 604 events\nM3 RR 1.63 (1.37-1.93)",
    "First eligible enhanced\nface-to-face wave (2006-2018)\nAge 50-80; PEF measured\nNo reported CLD",
    "Next biennial interview\nFirst reported physician\ndiagnosis; ~2-year interval",
    "N=17,030; 280 events\nM3 RR 3.15 (2.32-4.28)"
  ), stringsAsFactors = FALSE
)
box_df$fill <- ifelse(box_df$cohort == "CHARLS", palette["light_blue"], palette["light_orange"])

p_design <- ggplot() +
  geom_segment(data = data.frame(y = c(2.7, 1.25)),
               aes(x = 2.75, xend = 3.95, y = y, yend = y),
               arrow = arrow(length = grid::unit(2.5, "mm")), linewidth = 0.55, colour = palette["grey"]) +
  geom_segment(data = data.frame(y = c(2.7, 1.25)),
               aes(x = 6.05, xend = 7.25, y = y, yend = y),
               arrow = arrow(length = grid::unit(2.5, "mm")), linewidth = 0.55, colour = palette["grey"]) +
  geom_rect(data = box_df, aes(xmin = x - 1.05, xmax = x + 1.05, ymin = y - .55, ymax = y + .55, fill = fill),
            colour = palette["navy"], linewidth = 0.35) +
  geom_text(data = box_df, aes(x = x, y = y, label = label), size = 2.4, lineheight = .92, colour = "#1F2937") +
  annotate("text", x = .45, y = 2.7, label = "CHARLS", hjust = 1, fontface = "bold", size = 3.0, colour = palette["blue"]) +
  annotate("text", x = .45, y = 1.25, label = "HRS", hjust = 1, fontface = "bold", size = 3.0, colour = palette["orange"]) +
  annotate("rect", xmin = .55, xmax = 9.45, ymin = .10, ymax = .68, fill = "#F3F6F9", colour = "#A6B9D0", linewidth = .3) +
  annotate("text", x = 5, y = .39,
           label = "Harmonized logic: sex-specific relative low PEF | age 50-80 | aligned M3 covariates\nsurvey-weighted modified Poisson | cohort-specific inference",
           size = 2.45, lineheight = .95, colour = "#374151") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 3.45), clip = "off") +
  labs(title = "Parallel longitudinal analyses across two healthcare settings",
       subtitle = "The exposure and analysis logic were harmonized; follow-up interval and diagnosis context were intentionally retained as cohort-specific features",
       caption = "CLD, chronic lung disease; PEF, peak expiratory flow; M3, fully adjusted model. Outcome denotes subsequent clinical recognition, not spirometry-confirmed incident COPD.") +
  theme_void(base_family = "Arial", base_size = 8.2) +
  theme(plot.title = element_text(face = "bold", size = 9.4, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 7.7, colour = "#4B5563", margin = margin(b = 6)),
        plot.caption = element_text(size = 6.7, colour = "#555555", hjust = 0, margin = margin(t = 5)),
        plot.margin = margin(7, 7, 7, 7))

save_pub(p_design, file.path(fig_dir, "Figure1_parallel_cohort_design"), 183, 105)

# Figure 3: absolute risks and risk differences.
absrisk <- read.csv(file.path(in_dir, "absolute_risk_results.csv"), check.names = FALSE)
rd <- read.csv(file.path(in_dir, "risk_difference_results.csv"), check.names = FALSE)
absrisk$pef_group <- factor(absrisk$pef_group, levels = c("Non-low PEF", "Low PEF"))
absrisk$cohort <- factor(absrisk$cohort, levels = c("CHARLS", "HRS"))
abs_long <- rbind(
  data.frame(absrisk[, c("cohort", "pef_group")], estimate = "Crude weighted", risk = absrisk$crude_risk, low = absrisk$crude_ci_low, high = absrisk$crude_ci_high),
  data.frame(absrisk[, c("cohort", "pef_group")], estimate = "M3 standardized", risk = absrisk$adjusted_risk, low = absrisk$adjusted_ci_low, high = absrisk$adjusted_ci_high)
)
abs_long$estimate <- factor(abs_long$estimate, levels = c("Crude weighted", "M3 standardized"))

p_risk <- ggplot(abs_long, aes(pef_group, risk, colour = estimate, shape = estimate)) +
  geom_errorbar(aes(ymin = low, ymax = high), position = position_dodge(width = .38), width = .10, linewidth = .45) +
  geom_point(position = position_dodge(width = .38), size = 2.1) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c("Crude weighted" = unname(palette["grey"]), "M3 standardized" = unname(palette["orange"]))) +
  scale_shape_manual(values = c("Crude weighted" = 1, "M3 standardized" = 16)) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x), limits = c(0, .14), expand = expansion(mult = c(0, .04))) +
  labs(x = NULL, y = "Cumulative risk (95% CI)", title = "a  Cohort-specific absolute risks") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 7.2), plot.title = element_text(size = 8.5))

rd$cohort_label <- factor(ifelse(rd$cohort == "CHARLS", "CHARLS (~3 years)", "HRS (~2 years)"),
                          levels = c("HRS (~2 years)", "CHARLS (~3 years)"))
p_rd <- ggplot(rd, aes(adjusted_risk_difference, cohort_label)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .4, colour = "#666666") +
  geom_errorbar(aes(xmin = adjusted_rd_ci_low, xmax = adjusted_rd_ci_high), orientation = "y", width = .16, linewidth = .55, colour = palette["blue"]) +
  geom_point(size = 2.6, colour = palette["orange"]) +
  geom_text(aes(x = adjusted_rd_ci_high + .0015, label = sprintf("%+.1f pp", 100 * adjusted_risk_difference)), hjust = 0, size = 2.6, colour = "#333333") +
  scale_x_continuous(labels = function(x) sprintf("%.0f", 100 * x), limits = c(0, .065), breaks = seq(0, .06, .02), expand = expansion(mult = c(.02, .18))) +
  labs(x = "Adjusted risk difference (percentage points)", y = NULL,
       title = "b  M3-standardized excess risk") +
  theme_pub() +
  theme(legend.position = "none", plot.title = element_text(size = 8.5))

fig_abs <- p_risk + p_rd + plot_layout(widths = c(1.55, 1)) +
  plot_annotation(
    title = "Low PEF identifies higher subsequent clinical-recognition risk in both cohorts",
    caption = "Risks are survey-weighted cumulative risks over each cohort's follow-up interval, not incidence rates.\nAdjusted estimates use M3 marginal standardization; bars show 95% CIs. pp, percentage points."
  ) & theme(plot.title = element_text(face = "bold", size = 9.4),
            plot.caption = element_text(size = 6.7, colour = "#555555", hjust = 0))

save_pub(fig_abs, file.path(fig_dir, "Figure3_absolute_risk_and_difference"), 183, 105)
write.csv(abs_long, file.path(fig_dir, "Figure3_source_absolute_risks.csv"), row.names = FALSE)
write.csv(rd, file.path(fig_dir, "Figure3_source_risk_difference.csv"), row.names = FALSE)

# Figure 4: continuous dose-response.
curve <- read.csv(file.path(in_dir, "spline_curve_source.csv"), check.names = FALSE)
tests <- read.csv(file.path(in_dir, "spline_tests.csv"), check.names = FALSE)
curve$cohort <- factor(curve$cohort, levels = c("CHARLS", "HRS"))
tests$label <- sprintf("Overall P %s\nNonlinear P %s",
                       ifelse(tests$p_overall < .001, "<0.001", sprintf("=%.3f", tests$p_overall)),
                       ifelse(tests$p_nonlinear < .001, "<0.001", sprintf("=%.3f", tests$p_nonlinear)))
tests$x <- 1.58; tests$y <- 7.1

p_spline <- ggplot(curve, aes(pef_residual_sd, rr)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = palette["light_blue"], alpha = .85) +
  geom_line(linewidth = .75, colour = palette["blue"]) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = .4, colour = "#666666") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = .35, colour = "#777777") +
  geom_text(data = tests, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 1, vjust = 1, size = 2.65, lineheight = .95, colour = "#333333") +
  facet_wrap(~cohort, nrow = 1) +
  scale_x_continuous(breaks = -3:2) +
  scale_y_log10(breaks = c(.5, 1, 2, 4, 8), labels = c("0.5", "1", "2", "4", "8")) +
  coord_cartesian(xlim = c(-3.1, 1.75), ylim = c(.35, 8)) +
  labs(title = "Continuous PEF dose-response across CHARLS and HRS",
       subtitle = "Restricted cubic spline of the sex-specific PEF residual; zero denotes the expected PEF level",
       x = "PEF residual (SD; lower values indicate lower-than-expected PEF)",
       y = "Adjusted risk ratio (95% CI; log scale)",
       caption = "Survey-weighted modified Poisson M3 models. Knots: weighted 5th, 35th, 65th and 95th percentiles.\nCurves span the 1st to 99th weighted percentiles and are referenced to residual z=0.") +
  theme_pub() +
  theme(legend.position = "none")

save_pub(p_spline, file.path(fig_dir, "Figure4_PEF_restricted_cubic_spline"), 183, 105)

cat("Created publication figures in", fig_dir, "\n")
