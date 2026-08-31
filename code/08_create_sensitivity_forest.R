suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# Figure contract
# Core conclusion: the association between low PEF and later clinical recognition
# of chronic lung disease persists across selection-bias, exposure-definition,
# measurement-quality, delayed-incidence and missing-data analyses.
# Archetype: quantitative robustness forest plot, faceted by cohort.
# Target/output: supplementary double-column figure, 183 x 138 mm; editable SVG/PDF,
# 600-dpi TIFF and PNG preview; source data exported as CSV.
# Reviewer risk: sensitivity estimates use different subsets and definitions, so
# the figure emphasizes direction and interval stability rather than rank ordering.

root <- normalizePath(Sys.getenv("PEF_PROJECT_ROOT", unset = "."), mustWork = FALSE)
sens_dir <- file.path(root, "outputs", "sensitivity")
formal_dir <- file.path(root, "outputs", "formal_analysis")
fig_dir <- file.path(sens_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sens <- read.csv(file.path(sens_dir, "sensitivity_results.csv"), check.names = FALSE)
main <- read.csv(file.path(formal_dir, "main_results.csv"), check.names = FALSE)
main <- subset(main, model == "M3" & exposure == "low_pef_adj")
main$analysis <- "formal_primary_M3"
main <- main[, c("cohort", "analysis", "exposure", "n", "events", "rr", "ci_low", "ci_high", "p")]

keep_labels <- c(
  formal_primary_M3 = "Formal primary M3",
  loss_to_followup_IPCW = "Loss-to-follow-up IPCW",
  alternative_exposure_low_pef_weighted = "Survey-weighted residual q20",
  alternative_exposure_low_pef_q10 = "Stricter residual q10",
  alternative_exposure_low_pef_raw_weighted = "Sex-specific raw PEF q20",
  PEF_50_to_800 = "PEF restricted to 50-800 L/min",
  full_effort_only = "Full effort only",
  repeatable_top2_within_40 = "Top two attempts within 40 L/min",
  repeatable_top2_within_10pct = "Top two attempts within 10%",
  exclude_baseline_lung_medication = "Exclude baseline lung medication",
  four_year_delayed_incidence = "Four-year delayed incidence",
  five_year_delayed_incidence = "Five-year delayed incidence",
  persistent_two_consecutive_reports = "Persistent report at two interviews",
  multiple_imputation_m10 = "Multiple imputation (m=10)"
)

d <- rbind(main, subset(sens, exposure %in% c("low_pef_adj", "low_pef_weighted", "low_pef_q10", "low_pef_raw_weighted")))
d <- d[d$analysis %in% names(keep_labels), ]
d$label <- unname(keep_labels[d$analysis])
order_levels <- rev(unname(keep_labels))
d$label <- factor(d$label, levels = order_levels)
d$cohort <- factor(d$cohort, levels = c("CHARLS", "HRS"))

write.csv(d, file.path(fig_dir, "sensitivity_forest_source.csv"), row.names = FALSE)

p <- ggplot(d, aes(x = rr, y = label)) +
  geom_vline(xintercept = 1, colour = "#6B6B6B", linewidth = .4, linetype = 2) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = .18,
                linewidth = .55, colour = "#4C78A8") +
  geom_point(size = 2.15, shape = 16, colour = "#D55E00") +
  facet_wrap(~cohort, nrow = 1, scales = "free_y") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4.5, 6), limits = c(.9, 6.2)) +
  labs(x = "Adjusted risk ratio (log scale)", y = NULL,
       title = "Robustness of low PEF associations across sensitivity analyses",
       caption = "Survey-weighted modified Poisson models. Points show adjusted risk ratios; bars show 95% CIs.\nIPCW, inverse-probability-of-censoring weighting; PEF, peak expiratory flow; q10/q20, lowest weighted decile/quintile.") +
  theme_classic(base_size = 8.2, base_family = "Arial") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = "black", size = 7.3),
        strip.background = element_rect(fill = "#EAF0F7", colour = NA),
        strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold", size = 10),
        plot.caption = element_text(hjust = 0, size = 6.7, colour = "#3F3F3F"),
        panel.spacing.x = unit(8, "mm"), plot.margin = margin(5, 6, 4, 5))

base <- file.path(fig_dir, "Supplementary_Figure_sensitivity_forest")
svglite::svglite(paste0(base, ".svg"), width = 183 / 25.4, height = 138 / 25.4)
print(p); dev.off()
cairo_pdf(paste0(base, ".pdf"), width = 183 / 25.4, height = 138 / 25.4, family = "Arial")
print(p); dev.off()
ragg::agg_tiff(paste0(base, "_600dpi.tiff"), width = 183 / 25.4, height = 138 / 25.4, units = "in", res = 600,
               background = "white")
print(p); dev.off()
ragg::agg_png(paste0(base, "_preview.png"), width = 183, height = 138, units = "mm", res = 180,
              background = "white")
print(p); dev.off()
cat(base, "\n")
