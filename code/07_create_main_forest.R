suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Figure contract
# Conclusion: lower PEF precedes later clinical recognition of chronic lung
# disease in both cohorts; magnitude is stronger in HRS, while prespecified
# sex, age and smoking interactions are not supported.
# Archetype: two-panel quantitative forest plot.
# Export: 183 mm double-column SVG/PDF/TIFF plus PNG preview.
# Review risks: cross-cohort heterogeneity, different follow-up structures,
# and over-interpretation of within-subgroup significance.

root <- normalizePath(Sys.getenv("PEF_PROJECT_ROOT", unset = "."), mustWork = FALSE)
in_dir <- file.path(root, "outputs", "formal_analysis")
out_dir <- file.path(in_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

main <- read.csv(file.path(in_dir, "main_results.csv"), check.names = FALSE)
sub <- read.csv(file.path(in_dir, "subgroup_results.csv"), check.names = FALSE)

cols <- c(CHARLS = "#0072B2", HRS = "#D55E00")
shapes <- c(CHARLS = 16, HRS = 15)
pd <- position_dodge(width = 0.45)

main3 <- subset(main, model == "M3")
main3$label <- ifelse(main3$exposure == "low_pef_adj", "Low PEF", "PEF per 1-SD lower")
main3$label <- factor(main3$label, levels = c("PEF per 1-SD lower", "Low PEF"))

p_main <- ggplot(main3, aes(x = rr, y = label, colour = cohort, shape = cohort)) +
  geom_vline(xintercept = 1, linewidth = 0.45, colour = "#777777", linetype = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.14,
                 linewidth = 0.65, position = pd) +
  geom_point(size = 2.5, position = pd) +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4.5), limits = c(0.9, 5.0)) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  labs(x = "Adjusted risk ratio (log scale)", y = NULL,
       title = "A  Main associations", colour = NULL, shape = NULL) +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = "black"),
        legend.position = "top", legend.justification = "right",
        plot.title = element_text(face = "bold", size = 10),
        plot.margin = margin(4, 6, 4, 4))

sub_low <- subset(sub, exposure == "low_pef_adj")
sub_low$section <- factor(sub_low$modifier,
                          levels = c("sex_f", "age_group", "smoking_f"),
                          labels = c("Sex", "Age, years", "Smoking"))
sub_low$level_label <- factor(
  sub_low$level,
  levels = c("Current", "Former", "Never", "65-80", "50-64", "Female", "Male")
)

p_sub <- ggplot(sub_low, aes(x = rr, y = level_label, colour = cohort, shape = cohort)) +
  geom_vline(xintercept = 1, linewidth = 0.45, colour = "#777777", linetype = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.14,
                 linewidth = 0.6, position = pd) +
  geom_point(size = 2.25, position = pd) +
  facet_grid(section ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4.5, 6.5), limits = c(0.9, 7.0)) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  labs(x = "Adjusted risk ratio for low PEF (log scale)", y = NULL,
       title = "B  Prespecified subgroup analyses", colour = NULL, shape = NULL) +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = "black"),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", colour = "#333333"),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 10),
        plot.margin = margin(4, 6, 4, 4))

caption <- paste0(
  "Survey-weighted modified Poisson models (M3). Points are risk ratios and bars are 95% CIs. ",
  "\nInteraction P values (sex/age/smoking): CHARLS 0.315/0.873/0.270; ",
  "HRS 0.881/0.547/0.699. PEF, peak expiratory flow."
)

fig <- (p_main | p_sub) +
  plot_layout(widths = c(0.82, 1.18)) +
  plot_annotation(caption = caption,
                  theme = theme(plot.caption = element_text(size = 7.3, hjust = 0,
                                                            colour = "#333333",
                                                            margin = margin(t = 5)),
                                plot.background = element_rect(fill = "white", colour = NA)))

svg_file <- file.path(out_dir, "Figure_dual_cohort_forest.svg")
pdf_file <- file.path(out_dir, "Figure_dual_cohort_forest.pdf")
tif_file <- file.path(out_dir, "Figure_dual_cohort_forest_600dpi.tiff")
png_file <- file.path(out_dir, "Figure_dual_cohort_forest_preview.png")

ggsave(svg_file, fig, width = 183, height = 126, units = "mm", device = svglite::svglite)
ggsave(pdf_file, fig, width = 183, height = 126, units = "mm", device = cairo_pdf)
ragg::agg_tiff(tif_file, width = 183, height = 126, units = "mm", res = 600,
               background = "white")
print(fig)
dev.off()
ragg::agg_png(png_file, width = 183, height = 126, units = "mm", res = 180,
              background = "white")
print(fig)
dev.off()

cat(paste(svg_file, pdf_file, tif_file, png_file, sep = "\n"), "\n")
