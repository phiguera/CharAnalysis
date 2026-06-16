# plot_prob_peak.R
# Four-panel figure:
#   (a) C_interp, C_bkg, peak ID  -- panel (a) from char_plot_sni()
#   (b) SNI time series           -- panel (c) from char_plot_sni()
#   (c) Bubbles scaled by P(peak) -- per-peak summary
#   (d) P(peak) grey bars         -- ensemble probability of peak
#
# Requires:
#   ensemble  -- object in workspace from char_run_ensemble.R
#   out       -- CharAnalysis output from reference run (see below)

library(CharAnalysis)
library(ggplot2)
library(patchwork)

# ---- reference run (if `out` not already in workspace) ------------------
params_file <- "CharAnalysis_2_0_R/inst/validation/CH10_charParams.csv"
devtools::load_all("CharAnalysis_2_0_R")   # load dev version with plots = FALSE
out <- CharAnalysis(params_file, plots = FALSE)

# ---- replicate char_plot_sni panel (a) ----------------------------------
# Variable names and logic match char_plot_results.R exactly.
charcoal      <- out$charcoal
char_thresh   <- out$char_thresh
peak_analysis <- out$peak_analysis
pretreatment  <- out$pretreatment
post          <- out$post
site          <- out$site

zone_div    <- pretreatment$zoneDiv
transform   <- pretreatment$transform
cPeak       <- peak_analysis$cPeak
thresh_type <- peak_analysis$threshType

x  <- charcoal$ybpI
y  <- charcoal$accI
y2 <- charcoal$accIS
y3 <- char_thresh$pos[, ncol(char_thresh$pos)]
y4 <- if (is.matrix(char_thresh$neg)) char_thresh$neg[, 1L] else rep(char_thresh$neg[1L], length(x))

t_pos_line <- if (cPeak == 1L) y2 + y3 else y2 * y3
t_neg_line <- if (cPeak == 1L) y2 + y4 else y2 * y4

ccp      <- post$CharcoalCharPeaks
T_thresh <- ncol(ccp)
y_max    <- max(y, na.rm = TRUE)

lev <- c(0.78, 0.85, 0.92)
if (T_thresh < 3L) lev <- lev[(4L - T_thresh):3L]
sel_col <- T_thresh

has_ggtext  <- requireNamespace("ggtext", quietly = TRUE)
y_lbl_axis  <- if (has_ggtext) ggtext::element_markdown() else element_text()

y_lbl <- if (has_ggtext) {
  switch(as.character(transform),
         "1" = "log CHAR (# cm<sup>-2</sup> yr<sup>-1</sup>)",
         "2" = "ln CHAR (# cm<sup>-2</sup> yr<sup>-1</sup>)",
              "CHAR (# cm<sup>-2</sup> yr<sup>-1</sup>)")
} else {
  "CHAR (# cm^-2 yr^-1)"
}

x_ticks       <- seq(0, max(zone_div, na.rm = TRUE), by = 1000)
x_minor_ticks <- seq(0, max(zone_div, na.rm = TRUE), by = 500)

df_bar  <- data.frame(x = x, y = y)
df_back <- data.frame(x = x, y2 = y2)
df_tpos <- data.frame(x = x, t = t_pos_line)
df_tneg <- data.frame(x = x, t = t_neg_line)

pa <- ggplot(df_bar, aes(x = x, y = y)) +
  geom_col(fill = "black", colour = "black", width = mean(diff(x), na.rm = TRUE)) +
  geom_line(data = df_back, aes(x = x, y = y2),
            colour = "grey50", linewidth = 1.2) +
  geom_line(data = df_tpos, aes(x = x, y = t),
            colour = "red", linewidth = 0.6) +
  geom_line(data = df_tneg, aes(x = x, y = t),
            colour = "red", linewidth = 0.6)

# Zone dividers
if (length(zone_div) > 2L) {
  for (zd in zone_div[-c(1L, length(zone_div))]) {
    pa <- pa + annotate("segment",
                        x = zd, xend = zd,
                        y = y_max * 1.01, yend = y_max * 1.10,
                        colour = "grey50", linewidth = 1.0)
  }
  for (z in seq_along(zone_div[-1L])) {
    mid_x <- mean(zone_div[z:(z + 1L)])
    pa <- pa + annotate("text",
                        x = mid_x, y = y_max * 1.05,
                        label = paste0("Zone ", z),
                        hjust = 0.5, size = 2.5)
  }
}

# Peak markers: grey dots for alternative thresholds, black + for selected
tv          <- peak_analysis$threshValues
non_sel_idx <- setdiff(seq_len(T_thresh), sel_col)
non_sel_tv  <- tv[non_sel_idx]
show_cols   <- utils::head(non_sel_idx[order(non_sel_tv)], 3L)

for (rank in seq_along(show_cols)) {
  j      <- show_cols[rank]
  lev_j  <- lev[rank]
  pk_idx <- which(ccp[, j] > 0)
  if (length(pk_idx) > 0L) {
    df_pk <- data.frame(x = x[pk_idx], y = y_max * lev_j)
    pa <- pa + geom_point(data = df_pk, aes(x = x, y = y),
                          shape = 16L, size = 1, colour = "grey50")
  }
}

sel_tv     <- tv[sel_col]
sorted_tv  <- non_sel_tv[order(non_sel_tv)][seq_along(show_cols)]
match_rank <- which(sorted_tv == sel_tv)
if (length(match_rank) == 0L) match_rank <- ceiling(length(show_cols) / 2L)
mIndex     <- lev[match_rank[1L]]

pk_idx_sel <- which(ccp[, sel_col] > 0)
if (length(pk_idx_sel) > 0L) {
  df_pk_sel <- data.frame(x = x[pk_idx_sel], y = y_max * mIndex)
  pa <- pa +
    geom_point(data = df_pk_sel, aes(x = x, y = y),
               shape = 3L, size = 2.5, colour = "white", stroke = 1.5) +
    geom_point(data = df_pk_sel, aes(x = x, y = y),
               shape = 3L, size = 2.5, colour = "black", stroke = 0.9)
}

pa <- pa +
  scale_x_reverse(limits = c(max(zone_div), min(zone_div)),
                  breaks = x_ticks, minor_breaks = x_minor_ticks,
                  labels = x_ticks / 1000,
                  expand = expansion(0)) +
  guides(x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(limits = c(0, y_max * 1.15), expand = expansion(0)) +
  theme_classic() +
  theme(axis.title.y = y_lbl_axis,
        axis.title.x = element_blank(),
        axis.text.x  = element_blank()) +
  labs(title = paste0(site,
                      " (a) C_interpolated, C_background, and peak ID (+)"),
       y = y_lbl)

if (has_ggtext)
  pa <- pa + theme(plot.title = ggtext::element_markdown())

# ---- panel (b): SNI time series (panel c from char_plot_sni) -------------
sni_val <- char_thresh$SNI
sni_series <- if (length(sni_val) == 1L) rep(sni_val, length(x)) else sni_val
df_sni <- data.frame(x = x, sni = sni_series)

y_tick_sni <- pretty(c(0, max(sni_series, na.rm = TRUE)), n = 4)
y_lim_sni  <- c(0, max(y_tick_sni))

p_sni <- ggplot(df_sni, aes(x = x, y = sni)) +
  geom_line(colour = "black", linewidth = 0.7) +
  geom_hline(yintercept = 3, linetype = "dashed") +
  scale_x_reverse(limits = c(max(zone_div), min(zone_div)),
                  breaks = x_ticks, minor_breaks = x_minor_ticks,
                  labels = x_ticks / 1000,
                  expand = expansion(0)) +
  guides(x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(limits = y_lim_sni, breaks = y_tick_sni,
                     expand = expansion(0)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank()) +
  labs(y = "SNI",
       title = "(b) Local signal-to-noise index")

# ---- panel (c): P(peak) as grey bars + peak markers ----------------------
age       <- ensemble$age_grid
prob      <- ensemble$prob_peak
y_max_prob <- max(prob, na.rm = TRUE)

df_prob <- data.frame(age = age, prob = prob)

pb <- ggplot(df_prob, aes(x = age, y = prob)) +
  geom_col(fill = "grey50", colour = "grey50",
           width = mean(diff(age), na.rm = TRUE))

# Grey dots for alternative threshold peaks (same logic as panel a)
for (rank in seq_along(show_cols)) {
  j      <- show_cols[rank]
  lev_j  <- lev[rank]
  pk_idx <- which(ccp[, j] > 0)
  if (length(pk_idx) > 0L) {
    df_pk <- data.frame(age = x[pk_idx], y = y_max_prob * lev_j)
    pb <- pb + geom_point(data = df_pk, aes(x = age, y = y),
                          shape = 16L, size = 1, colour = "grey50")
  }
}

# Black + for selected (final) threshold peaks
if (length(pk_idx_sel) > 0L) {
  df_pk_sel_b <- data.frame(age = x[pk_idx_sel], y = y_max_prob * mIndex)
  pb <- pb +
    geom_point(data = df_pk_sel_b, aes(x = age, y = y),
               shape = 3L, size = 2.5, colour = "white", stroke = 1.5) +
    geom_point(data = df_pk_sel_b, aes(x = age, y = y),
               shape = 3L, size = 2.5, colour = "black", stroke = 0.9)
}

pb <- pb +
  scale_x_reverse(limits = c(max(zone_div), min(zone_div)),
                  breaks = x_ticks, minor_breaks = x_minor_ticks,
                  labels = x_ticks / 1000,
                  expand = expansion(0)) +
  guides(x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(limits = c(0, y_max_prob * 1.15), expand = expansion(0)) +
  theme_classic() +
  labs(x = "cal. yr BP (x 1000)",
       y = "P(peak)",
       title = sprintf("(d) Probability of peak  [n = %d chronologies]",
                       ensemble$n_iter))

# ---- panel (d): peaksFinal bubbles scaled by P(peak) ---------------------
# Reference peak ages from out
peak_ages <- x[pk_idx_sel]

# Look up ensemble P(peak) at each reference peak age (nearest grid point)
peak_prob <- sapply(peak_ages, function(a) {
  idx <- which.min(abs(ensemble$age_grid - a))
  ensemble$prob_peak[idx]
})

# Normalise so the highest-probability peak always maps to size = 1
peak_prob_norm <- peak_prob / max(peak_prob)

df_peaks <- data.frame(age = peak_ages, y = 1, prob_norm = peak_prob_norm)

pc <- ggplot(df_peaks, aes(x = age, y = y, size = prob_norm)) +
  geom_point(shape = 21, fill = "grey50", colour = "black") +
  scale_size_continuous(limits = c(0, 1), range = c(1, 8), guide = "none") +
  scale_x_reverse(limits = c(max(zone_div), min(zone_div)),
                  breaks = x_ticks, minor_breaks = x_minor_ticks,
                  labels = x_ticks / 1000,
                  expand = expansion(0)) +
  guides(x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(limits = c(0.5, 1.5), breaks = NULL,
                     expand = expansion(0)) +
  theme_classic() +
  theme(axis.ticks.y  = element_blank(),
        axis.title.y  = element_blank(),
        axis.line.y   = element_blank(),
        axis.title.x  = element_blank(),
        axis.text.x   = element_blank()) +
  labs(title = "(c) Reference peaks scaled by P(peak)")

# ---- combine with patchwork ---------------------------------------------
pa / p_sni / pc / pb + plot_layout(heights = c(3, 1, 0.5, 1.5))
