# plot_agedepth_vignette.R
# -----------------------------------------------------------------------
# Age-depth model figure for the vignette: three sites on one panel.
# Age (cal yr BP) on y-axis, reversed (present at top).
# Depth on x-axis, reversed (surface at right).
# Shared axis limits so accumulation rate slopes are directly comparable.
# -----------------------------------------------------------------------

library(ggplot2)

# ---- User settings -----------------------------------------------------
sites    <- c("CH10", "CO", "SI17")
rds_fmt  <- "CharAnalysis_2_0_R/tests/%s_ensemble_results.rds"
out_file <- "CharAnalysis_2_0_R/tests/agedepth_vignette.png"

site_colours <- c("CH10" = "#2166ac", "CO" = "#1b7837", "SI17" = "black")

# ---- Load chron_ci from each site's RDS --------------------------------
ci_list <- lapply(sites, function(s) {
  rds <- sprintf(rds_fmt, s)
  if (!file.exists(rds)) stop("RDS not found: ", rds)
  ens <- readRDS(rds)
  if (is.null(ens$chron_ci)) stop("chron_ci missing in RDS for site: ", s)
  df <- ens$chron_ci[order(ens$chron_ci$depth_cm), ]
  df$site <- s
  df
})
df_all <- do.call(rbind, ci_list)
df_all$site <- factor(df_all$site, levels = sites)

# ---- Shared axis limits ------------------------------------------------
x_hi <- ceiling(max(df_all$depth_cm, na.rm = TRUE) / 10) * 10
x_lo <- 0
y_hi <- ceiling(max(df_all$ci95_hi,  na.rm = TRUE) / 500) * 500
y_br <- seq(0, y_hi, by = 1000)

# ---- Plot --------------------------------------------------------------
p_ad <- ggplot(df_all, aes(x = depth_cm, y = median_age,
                            colour = site, fill = site)) +
  geom_ribbon(aes(ymin = ci95_lo, ymax = ci95_hi),
              alpha = 0.45, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = site_colours, name = NULL) +
  scale_fill_manual(  values = site_colours, name = NULL) +
  scale_x_reverse(expand = expansion(0)) +
  scale_y_reverse(breaks = y_br, expand = expansion(0)) +
  coord_cartesian(
    xlim = c(x_hi + 20, x_lo),   # deep end has 20 cm breathing room; surface at right edge
    ylim = c(y_hi + 100, -100)   # oldest age at bottom; slight padding top and bottom
  ) +
  labs(
    x = "Depth below mud-water interface (cm)",
    y = "Cal yr BP"
  ) +
  theme_classic() +
  theme(
    aspect.ratio         = 1,
    axis.text            = element_text(size = 12),
    axis.title           = element_text(size = 13),
    panel.grid.major     = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.minor     = element_blank(),
    legend.position      = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background    = element_blank(),
    legend.key.size      = unit(1.0, "lines"),
    legend.text          = element_text(size = 12)
  )

print(p_ad)

ggsave(out_file, plot = p_ad,
       width = 3, height = 3, units = "in", dpi = 300)
message("Saved: ", out_file)
