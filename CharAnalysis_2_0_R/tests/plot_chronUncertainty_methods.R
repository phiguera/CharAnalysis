# plot_chronUncertainty_methods.R
# -----------------------------------------------------------------------
# Methods illustration figure for chronological uncertainty analysis.
#
# Row 1    : Reference (median chronology) CHAR record — full panel-(a)
#            content with "+" at reference peaks (black).
# Rows 2-6 : Five random ensemble iterations — same panel-(a) content.
#            Detected peaks shown as filled circles with horizontal bars
#            spanning the matching window (±match_halfwin_k of nearest
#            reference peak). Black = matched; grey = unmatched.
#            Vertical grey lines mark reference peak ages in all panels.
# Row 7    : Panel-(d) equivalent from the five example iterations only.
#
# Requires in workspace:
#   ensemble        -- CharEnsemble object
#   out             -- reference CharAnalysis run
#   ref_peak_ages   -- from run_ensemble_analysis.R
#   match_halfwin_k -- from run_ensemble_analysis.R
# -----------------------------------------------------------------------

library(ggplot2)
library(patchwork)

# ---- user settings -----------------------------------------------------
set.seed(42)
n_ens    <- 5
sel_iter <- sample(seq_len(ensemble$n_iter), n_ens)

# Optional y-axis ceiling for CHAR panels (all rows).
# Add site-specific overrides below; NA = data-driven 99.9th percentile.
# CH10 has highly skewed CHAR values that compress the y-axis without a cap.
y_ceil_char <- switch(out$site,
  "CH10" = 50,
  NA     # default: computed from data
)

# Optional x-axis crop: oldest age to display (cal yr BP).
# NA = full record. Use e.g. 3000 to zoom into the most recent 3000 years.
x_crop <- switch(out$site,
  "CH10" = 3000,
  NA     # default: full record
)

# ---- shared axis parameters --------------------------------------------
age_grid      <- ensemble$age_grid
bar_width_ens <- mean(diff(age_grid))
bar_width_ref <- mean(diff(out$charcoal$ybpI))
zone_div      <- out$pretreatment$zoneDiv
x_lo          <- min(zone_div)                          # youngest age (right side)
x_hi          <- if (!is.na(x_crop)) x_crop else max(zone_div)  # oldest age (left side)
x_lims        <- c(x_hi, x_lo)
x_ticks       <- seq(0, x_hi, by = 1000)
x_minor       <- seq(0, x_hi, by = 500)
cPeak         <- out$peak_analysis$cPeak

# ---- threshold interpolated onto ensemble age grid --------------------
x_ref  <- out$charcoal$ybpI
y3_ref <- out$char_thresh$pos[, ncol(out$char_thresh$pos)]
y4_ref <- if (is.matrix(out$char_thresh$neg)) {
  out$char_thresh$neg[, 1L]
} else {
  rep(out$char_thresh$neg[1L], length(x_ref))
}
thresh_pos <- approx(x_ref, y3_ref, xout = age_grid, rule = 2)$y
thresh_neg <- approx(x_ref, y4_ref, xout = age_grid, rule = 2)$y

# ---- consistent y ceiling across all panels ----------------------------
y_ceil <- if (!is.na(y_ceil_char)) {
  y_ceil_char
} else {
  max(
    quantile(ensemble$charPeak_matrix + ensemble$charBkg_matrix,
             probs = 0.999, na.rm = TRUE),
    quantile(out$charcoal$accI, probs = 0.999, na.rm = TRUE)
  )
}

# ---- y-axis label (avoid # comment issue) ------------------------------
y_axis_lbl <- expression(CHAR ~ ("#" ~ cm^{-2} ~ yr^{-1}))

# ---- shared x scales and theme for CHAR panels -------------------------
theme_char <- theme_classic(base_size = 9) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.title.y.left  = element_text(size = 8),
    axis.title.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank(),
    axis.text.y.right  = element_blank(),
    plot.title         = element_text(size = 8),
    plot.margin        = margin(1, 5, 1, 5)
  )

scale_x_char <- scale_x_reverse(
  limits       = x_lims,
  breaks       = x_ticks,
  minor_breaks = x_minor,
  expand       = expansion(0),
  labels       = NULL
)

scale_y_char <- function() {
  scale_y_continuous(
    limits = c(0, y_ceil * 1.15),  # extra headroom for zone labels on p_ref
    expand = expansion(0)
  )
}

# Vertical lines at reference peak ages — filtered to visible window
ref_peak_ages_vis <- ref_peak_ages[ref_peak_ages <= x_hi]
vlines_ref <- geom_vline(xintercept = ref_peak_ages_vis,
                         colour = "grey78", linewidth = 0.3)

# Zone boundary vlines (interior boundaries only) for all panels
inner_zones <- zone_div[-c(1L, length(zone_div))]
zone_vlines <- if (length(inner_zones) > 0L) {
  geom_vline(xintercept = inner_zones, colour = "grey50",
             linewidth = 0.8, linetype = "dashed")
} else NULL

# ========================================================================
# Row 1: Reference run
# ========================================================================
sel_col    <- ncol(out$post$CharcoalCharPeaks)
ref_pk_idx <- which(out$post$CharcoalCharPeaks[, sel_col] > 0)

t_pos_ref <- if (cPeak == 1L) out$charcoal$accIS + y3_ref else
               out$charcoal$accIS * y3_ref
t_neg_ref <- if (cPeak == 1L) out$charcoal$accIS + y4_ref else
               out$charcoal$accIS * y4_ref

p_ref <- ggplot(data.frame(age = x_ref, y = out$charcoal$accI),
                aes(x = age, y = y)) +
  geom_col(fill = "black", colour = NA, width = bar_width_ref) +
  geom_line(data = data.frame(age = x_ref, y = out$charcoal$accIS),
            aes(x = age, y = y), colour = "grey50", linewidth = 1.0) +
  geom_line(data = data.frame(age = x_ref, y = t_pos_ref),
            aes(x = age, y = y), colour = "red", linewidth = 0.5) +
  geom_line(data = data.frame(age = x_ref, y = t_neg_ref),
            aes(x = age, y = y), colour = "red", linewidth = 0.5) +
  geom_point(
    data = data.frame(age = x_ref[ref_pk_idx]),
    aes(x = age, y = y_ceil * 0.91),
    shape = 3, size = 2.5, stroke = 1.1, colour = "black"
  ) +
  zone_vlines +
  scale_x_char +
  scale_y_char() +
  labs(x = NULL, y = NULL,
       title = paste0(out$site, ": CharAnalysis results (benchmark)")) +
  theme_char

# Zone labels on p_ref (centred in each zone, above data area)
if (length(inner_zones) > 0L) {
  for (zd in inner_zones) {
    p_ref <- p_ref + annotate("segment",
      x = zd, xend = zd,
      y = y_ceil * 1.01, yend = y_ceil * 1.10,
      colour = "grey50", linewidth = 1.0, linetype = "dashed")
  }
  for (z in seq_along(zone_div[-1L])) {
    mid_x <- mean(zone_div[z:(z + 1L)])
    p_ref <- p_ref + annotate("text",
      x = mid_x, y = y_ceil * 1.055,
      label = paste0("Zone ", z), hjust = 0.5, size = 2.5)
  }
}

# ========================================================================
# Rows 2-6: Ensemble iterations
# ========================================================================

# For each detected peak: classify matched/unmatched and get matching
# half-window. Bars are clamped to plot bounds so the oldest/youngest
# peaks are not silently dropped by scale_x_reverse.
classify_peaks <- function(iter_idx) {
  pk_idx <- which(ensemble$peaks_matrix[, iter_idx] == 1L)
  if (length(pk_idx) == 0L) {
    return(data.frame(age     = numeric(0),
                      matched = logical(0),
                      xmin    = numeric(0),
                      xmax    = numeric(0),
                      y_pos   = numeric(0)))
  }
  pk_ages   <- age_grid[pk_idx]
  nearest_k <- sapply(pk_ages, function(pa) which.min(abs(ref_peak_ages - pa)))
  matched   <- mapply(function(pa, k) abs(pa - ref_peak_ages[k]) <= match_halfwin_k[k],
                      pk_ages, nearest_k)
  halfwin   <- match_halfwin_k[nearest_k]
  data.frame(
    age     = pk_ages,
    matched = matched,
    xmin    = pmax(pk_ages - halfwin, x_lo),
    xmax    = pmin(pk_ages + halfwin, x_hi),
    y_pos   = y_ceil * 0.91
  )
}

panel_plots <- lapply(seq_along(sel_iter), function(j) {
  i     <- sel_iter[j]
  accI  <- ensemble$charBkg_matrix[, i] + ensemble$charPeak_matrix[, i]
  accIS <- ensemble$charBkg_matrix[, i]
  t_pos <- if (cPeak == 1L) accIS + thresh_pos else accIS * thresh_pos
  t_neg <- if (cPeak == 1L) accIS + thresh_neg else accIS * thresh_neg

  df_pk <- classify_peaks(i)

  # y-axis label on row 4 overall (= ensemble row j == 3)
  y_title <- if (j == 3) y_axis_lbl else NULL

  ggplot(data.frame(age = age_grid, y = accI), aes(x = age, y = y)) +
    vlines_ref +
    zone_vlines +
    geom_col(fill = "black", colour = NA, width = bar_width_ens) +
    geom_line(data = data.frame(age = age_grid, y = accIS),
              aes(x = age, y = y), colour = "grey50", linewidth = 1.0) +
    geom_line(data = data.frame(age = age_grid, y = t_pos),
              aes(x = age, y = y), colour = "red", linewidth = 0.5) +
    geom_line(data = data.frame(age = age_grid, y = t_neg),
              aes(x = age, y = y), colour = "red", linewidth = 0.5) +
    geom_errorbarh(
      data = df_pk,
      aes(xmin = xmin, xmax = xmax, y = y_pos, colour = matched),
      height = y_ceil * 0.06, linewidth = 0.7
    ) +
    geom_point(
      data = df_pk,
      aes(x = age, y = y_pos, colour = matched),
      shape = 16, size = 2.2
    ) +
    scale_colour_manual(
      values = c("TRUE" = "black", "FALSE" = "grey65"),
      guide  = "none"
    ) +
    scale_x_char +
    scale_y_char() +
    labs(x = NULL, y = y_title, title = paste0("Iteration i = ", i)) +
    theme_char
})

# ========================================================================
# Bottom panel: panel-(d) equivalent from 5 iterations only
# ========================================================================
sub_peaks <- ensemble$peaks_matrix[, sel_iter]
n_sub     <- length(sel_iter)

# ---- matched peaks (reference-aligned) --------------------------------
assigned_sub <- lapply(seq_along(ref_peak_ages), function(k) {
  win_k  <- match_halfwin_k[k]
  ages_k <- numeric(0)
  for (j in seq_len(n_sub)) {
    det_idx <- which(sub_peaks[, j] == 1L)
    if (length(det_idx) == 0L) next
    det_ages <- age_grid[det_idx]
    dists    <- abs(det_ages - ref_peak_ages[k])
    best     <- which.min(dists)
    if (dists[best] <= win_k) ages_k <- c(ages_k, det_ages[best])
  }
  ages_k
})

df_sub <- data.frame(
  ref_age  = ref_peak_ages,
  det_freq = sapply(assigned_sub, function(v) length(v) / n_sub * 100),
  ci95_lo  = sapply(assigned_sub, function(v)
    if (length(v) > 1) quantile(v, 0.025) else NA_real_),
  ci95_hi  = sapply(assigned_sub, function(v)
    if (length(v) > 1) quantile(v, 0.975) else NA_real_)
)

# ---- ensemble-only peaks (unmatched, using same logic as iteration panels) ----
# Use classify_peaks() directly so "unmatched" means the same thing here
# as in the grey circles plotted above (nearest-reference-peak rule).
ens_only_min_freq <- 40   # minimum detection frequency (%) to show in bottom panel
clust_win         <- mean(match_halfwin_k)

# Collect unmatched peak ages with iteration index so we can cap
# contributions at one detection per iteration per cluster.
ens_only <- data.frame(age = numeric(0), iter = integer(0))
for (j in seq_along(sel_iter)) {
  df_pk_j      <- classify_peaks(sel_iter[j])
  unmatched_ages <- df_pk_j$age[!df_pk_j$matched]
  if (length(unmatched_ages) > 0L)
    ens_only <- rbind(ens_only,
                      data.frame(age = unmatched_ages, iter = j))
}

# Greedy clustering: group unmatched peaks within clust_win of each other.
# Detection frequency counts unique iterations per cluster (not total peaks),
# so one iteration contributing two nearby peaks counts only once.
cluster_unmatched <- function(ens_df, win, n_iter) {
  if (nrow(ens_df) == 0L)
    return(data.frame(clust_age = numeric(0), det_freq  = numeric(0),
                      ci95_lo   = numeric(0), ci95_hi   = numeric(0)))
  ord       <- order(ens_df$age)
  ages      <- ens_df$age[ord]
  iters     <- ens_df$iter[ord]
  clust     <- integer(length(ages))
  k         <- 1L
  clust[1L] <- k
  for (i in seq_along(ages)[-1L]) {
    if (ages[i] - ages[i - 1L] <= win) {
      clust[i] <- k
    } else {
      k <- k + 1L
      clust[i] <- k
    }
  }
  clust_ids <- unique(clust)
  data.frame(
    clust_age = sapply(clust_ids, function(id) median(ages[clust == id])),
    # count unique iterations, not total peaks
    det_freq  = sapply(clust_ids,
                       function(id) length(unique(iters[clust == id])) / n_iter * 100),
    ci95_lo   = sapply(clust_ids, function(id) {
      # one representative age per iteration (median within iteration)
      v <- tapply(ages[clust == id], iters[clust == id], median)
      if (length(v) > 1L) quantile(v, 0.025) else NA_real_
    }),
    ci95_hi   = sapply(clust_ids, function(id) {
      v <- tapply(ages[clust == id], iters[clust == id], median)
      if (length(v) > 1L) quantile(v, 0.975) else NA_real_
    })
  )
}

df_ens_only <- cluster_unmatched(ens_only, clust_win, n_sub)
df_ens_only <- df_ens_only[df_ens_only$det_freq >= ens_only_min_freq, ]

# ---- shared y limits (include both matched and ensemble-only points) ---
all_freqs  <- c(df_sub$det_freq, df_ens_only$det_freq)
y_ceil_bot <- max(c(all_freqs, 0), na.rm = TRUE)
y_floor    <- 10                  # fixed: 1/5 = 20% is min possible; 10 gives room below
y_top      <- y_ceil_bot * 1.10  # 10% headroom above max point
y_breaks   <- seq(20, floor(y_ceil_bot / 20) * 20, by = 20)
cap_ht     <- (y_top - y_floor) * 0.02

p_bottom <- ggplot(df_sub, aes(x = ref_age, y = det_freq)) +
  vlines_ref +
  zone_vlines +
  geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi),
                width = cap_ht, linewidth = 0.7, colour = "black",
                na.rm = TRUE) +
  geom_point(shape = 21, fill = "black", colour = "black", size = 2.2) +
  # Ensemble-only clusters: timing CI bars + open grey circles
  {if (nrow(df_ens_only) > 0)
    geom_errorbarh(data = df_ens_only,
                   aes(xmin = ci95_lo, xmax = ci95_hi, y = det_freq),
                   height = cap_ht, linewidth = 0.7, colour = "grey50",
                   na.rm = TRUE, inherit.aes = FALSE)
  } +
  {if (nrow(df_ens_only) > 0)
    geom_point(data = df_ens_only, aes(x = clust_age, y = det_freq),
               shape = 21, fill = "white", colour = "grey50",
               size = 2.2, stroke = 0.8, inherit.aes = FALSE)
  } +
  geom_hline(yintercept = 50, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  scale_x_reverse(
    limits       = x_lims,
    breaks       = x_ticks,
    minor_breaks = x_minor,
    labels       = x_ticks / 1000,
    expand       = expansion(0)
  ) +
  scale_y_continuous(
    limits = c(y_floor, y_top),
    breaks = y_breaks,
    expand = expansion(0)
  ) +
  guides(x = guide_axis(minor.ticks = TRUE)) +
  labs(
    x        = "cal yr BP (× 1000)",
    y        = "Detection\nfreq. (%)",
    subtitle = sprintf(
      "Detection frequency from %d example iterations  |  filled = ref-matched, open = ensemble-only",
      n_sub)
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.subtitle = element_text(size = 8),
    plot.margin   = margin(3, 5, 5, 5)
  )

# ========================================================================
# Assemble
# ========================================================================
fig_methods <- wrap_plots(
  c(list(p_ref), panel_plots, list(p_bottom)),
  ncol    = 1,
  heights = c(rep(1, n_ens + 1), 0.5)
)

print(fig_methods)

ggsave(
  filename = sprintf("CharAnalysis_2_0_R/tests/%s_chronUncertainty_methods.png",
                     out$site),
  plot   = fig_methods,
  width  = 6, height = 10, units = "in", dpi = 300
)
message("Saved: ", sprintf("CharAnalysis_2_0_R/tests/%s_chronUncertainty_methods.png",
                            out$site))
