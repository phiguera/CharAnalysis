#' CharAnalysis output figures
#'
#' Nine output figures mirroring the MATLAB CharAnalysis v2.0 plots.
#' Function names follow R snake_case conventions.
#'
#' \describe{
#'   \item{[char_plot_raw()]}{Figure 1 (allFigures only):
#'     C_raw, C_interpolated, and C_background smoothing options.}
#'   \item{[char_plot_thresh_diag()]}{Figure 2 (allFigures only):
#'     Local threshold determination diagnostics (5x5 window grid).}
#'   \item{[char_plot_peaks()]}{Figure 3: Resampled CHAR with
#'     background trend (top) and C_peak with thresholds and peak markers (bottom).}
#'   \item{[char_plot_sni()]}{Figure 4: Sensitivity to alternative
#'     thresholds and signal-to-noise index.}
#'   \item{[char_plot_cumulative()]}{Figure 5: Cumulative peaks through time.}
#'   \item{[char_plot_fri()]}{Figure 6: FRI distributions by
#'     zone with Weibull model fits.}
#'   \item{[char_plot_fire_history()]}{Figure 7: Peak magnitude (top),
#'     FRIs through time with smoothed FRI and CI ribbon (middle), and smoothed
#'     fire frequency (bottom).}
#'   \item{[char_plot_zones()]}{Figure 8: Between-zone comparisons
#'     of raw CHAR distributions (CDF and box plots).}
#'   \item{[char_plot_all()]}{Convenience wrapper: produces all figures and
#'     optionally saves them as PDF files.}
#' }
#'
#' @name char_plot
#' @aliases char_plot_raw char_plot_thresh_diag
#'   char_plot_peaks char_plot_sni
#'   char_plot_cumulative char_plot_fri
#'   char_plot_fire_history char_plot_zones char_plot_all
#'
#' @param out Named list returned by [CharAnalysis()].  Must contain
#'   \code{charcoal}, \code{pretreatment}, \code{peak_analysis},
#'   \code{char_thresh}, \code{post}, and \code{site}.
#' @param save   Logical.  If \code{TRUE}, each figure is saved as a PDF in
#'   \code{out_dir}.  Default \code{FALSE}.
#' @param out_dir Directory for saved PDFs.  Default: current working directory.
#' @param width,height PDF dimensions in inches.  Defaults: 11 x 8.5.
#'
#' @return
#'   Individual figure functions each return a \pkg{patchwork} / \pkg{ggplot2}
#'   object.  \code{char_plot_all()} returns a named list of all figure objects.
#'
#' @details
#'   Requires the \pkg{ggplot2} package.  Multi-panel layout uses
#'   \pkg{patchwork} if available; otherwise panels are printed separately
#'   with a message.
#'
#' @seealso [CharAnalysis()], [char_post_process()]
#'
#' @examples
#' \donttest{
#'   # Run pipeline on the bundled example dataset, then plot:
#'   params_file <- system.file("validation", "CO_charParams.csv",
#'                              package = "CharAnalysis")
#'   out <- CharAnalysis(params_file)
#'   char_plot_peaks(out)
#'   char_plot_fire_history(out)
#'   # Individual figures can also be called directly:
#'   char_plot_sni(out)
#'   char_plot_fri(out)
#'   # Save all figures to PDF in a temporary directory:
#'   char_plot_all(out, save = TRUE, out_dir = tempdir())
#' }
NULL

# Suppress R CMD check NOTEs for ggplot2 aes() column names used via NSE.
utils::globalVariables(c(
  "x", "y", "y_lo", "y_hi",           # generic aesthetics
  "accI", "bkg", "peak", "pos", "neg", # char_plot_peaks
  "age", "acc", "smooth", "method",    # char_plot_raw
  "mag",                               # char_plot_fire_history
  "ff",                                # char_plot_fire_history (fire freq)
  "lbl",                               # char_plot_fri (zone labels)
  "xpos", "lo", "hi", "mfri", "sni",  # char_plot_sni
  "zone"                               # char_plot_zones
))

# =============================================================================
# Internal helpers
# =============================================================================

# Use ggtext::element_markdown() for titles when available so that
# HTML <sub>...</sub> tags render as proper subscripts.  Falls back to
# plain element_text() if ggtext is not installed (tags are stripped first).
.has_ggtext <- function() requireNamespace("ggtext", quietly = TRUE)

.char_theme <- function(base_size = 10) {
  title_el <- if (.has_ggtext()) {
    ggtext::element_markdown(size = base_size, face = "bold")
  } else {
    ggplot2::element_text(size = base_size, face = "bold")
  }
  # Axis-title element: use element_markdown() when ggtext is available so
  # that HTML <sup>/<sub> tags render as proper super/subscripts.  Plain text
  # passes through unchanged, so this is safe for all axis labels.
  axis_title_el <- if (.has_ggtext()) {
    ggtext::element_markdown(size = base_size)
  } else {
    ggplot2::element_text(size = base_size)
  }
  ggplot2::theme_classic(base_size = base_size) +
  ggplot2::theme(
    axis.line         = ggplot2::element_line(colour = "black"),
    panel.grid        = ggplot2::element_blank(),
    plot.title        = title_el,
    axis.title.x      = axis_title_el,
    axis.title.y      = axis_title_el,
    axis.ticks        = ggplot2::element_line(colour = "black"),
    axis.ticks.length = ggplot2::unit(3, "pt")
  )
}

# Format a title string: keep HTML when ggtext is available, strip tags otherwise.
.title <- function(html) {
  if (.has_ggtext()) html else gsub("<[^>]+>", "", html)
}

# Return a title theme element of the correct class (must match .char_theme).
# Always call this instead of element_text() when overriding plot.title so
# that patchwork can merge elements across panels without a class conflict.
.title_el <- function(size = 10, face = "bold", hjust = 0) {
  if (.has_ggtext()) {
    ggtext::element_markdown(size = size, face = face, hjust = hjust)
  } else {
    ggplot2::element_text(size = size, face = face, hjust = hjust)
  }
}

.char_xscale <- function(zone_div) {
  x_max   <- max(zone_div)
  x_break <- seq(0, x_max, by = 1000)
  x_label <- x_break / 1000
  list(
    ggplot2::scale_x_reverse(
      limits = c(x_max, zone_div[1L]),
      breaks = x_break,
      labels = x_label
    ),
    ggplot2::xlab("cal. yr BP (x 1000)")
  )
}

.zone_lines <- function(zone_div, y_lo, y_hi) {
  # Vertical dashed lines at interior zone boundaries
  inner <- zone_div[-c(1L, length(zone_div))]
  if (length(inner) == 0L) return(NULL)
  lapply(inner, function(z) {
    ggplot2::annotate("segment",
                      x = z, xend = z,
                      y = y_lo, yend = y_hi,
                      colour = "grey50", linewidth = 0.8, linetype = "dashed")
  })
}

.zone_labels <- function(zone_div, y_pos, base_size = 10) {
  # Text labels centred in each zone (e.g. "Zone 1", "Zone 2")
  n_zones <- length(zone_div) - 1L
  if (n_zones <= 1L) return(NULL)
  mids   <- (zone_div[-length(zone_div)] + zone_div[-1L]) / 2
  labels <- paste0("Zone ", seq_len(n_zones))
  lapply(seq_len(n_zones), function(z) {
    ggplot2::annotate("text",
                      x = mids[z], y = y_pos,
                      label = labels[z],
                      size  = base_size / ggplot2::.pt,
                      hjust = 0.5)
  })
}

.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required: install.packages('ggplot2')")
}

.combine_panels <- function(panels, ncol = 1L, heights = NULL) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- Reduce(`/`, panels)
    if (!is.null(heights))
      p <- p + patchwork::plot_layout(heights = heights)
    return(p)
  }
  # Fallback: print each panel in sequence
  message("Install 'patchwork' for combined multi-panel layout: ",
          "install.packages('patchwork')")
  invisible(panels[[length(panels)]])
}

# =============================================================================
# char_plot_peaks  --  Figure 3: C_interp / C_background / C_peak
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_peaks <- function(out) {

  .require_ggplot2()

  charcoal     <- out$charcoal
  pretreatment <- out$pretreatment
  peak_analysis <- out$peak_analysis
  char_thresh  <- out$char_thresh
  post         <- out$post
  site         <- out$site

  zone_div  <- pretreatment$zoneDiv
  yr_interp <- pretreatment$yrInterp
  x         <- charcoal$ybpI
  r         <- yr_interp

  # ---------- Panel (a): C_interpolated and C_background -------------------
  df_a <- data.frame(
    x    = x,
    accI = charcoal$accI,
    bkg  = charcoal$accIS
  )

  y_max_a <- max(df_a$accI, na.rm = TRUE)
  y_lim_a <- c(0, 1.1 * y_max_a)

  y_label_a <- if (!is.null(pretreatment$transform) &&
                   pretreatment$transform > 0L) {
    expression(paste("C"["interp"], " (transformed pieces cm"^{-2}*" yr"^{-1}*")"))
  } else {
    expression(paste("C"["interp"], " (pieces cm"^{-2}*" yr"^{-1}*")"))
  }

  # Smoothing window label (yr): take last element of smoothing$yr vector
  smooth_yr <- out$smoothing$yr
  smooth_yr_label <- if (!is.null(smooth_yr)) {
    paste0(smooth_yr[length(smooth_yr)], "-yr")
  } else {
    "low-frequency"
  }

  pa <- ggplot2::ggplot(df_a) +
    ggplot2::geom_col(ggplot2::aes(x = x, y = accI),
                      fill = "black", colour = "black", width = r) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = bkg),
                       colour = "grey50", linewidth = 1.2) +
    .zone_lines(zone_div, y_max_a * 1.01, y_max_a * 1.09) +
    .zone_labels(zone_div, y_max_a * 1.05) +
    .char_xscale(zone_div) +
    ggplot2::scale_y_continuous(limits = y_lim_a, expand = c(0, 0)) +
    ggplot2::ylab(y_label_a) +
    ggplot2::labs(
      title = .title(paste0(
        "(a) ", site,
        ": C<sub>interpolated</sub> (", yr_interp, " yr) and",
        " C<sub>background</sub> defined by ", smooth_yr_label, "-yr trends"))
    ) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +
    .char_theme()

  # ---------- Panel (b): C_peak and thresholds ------------------------------
  df_b <- data.frame(
    x     = x,
    peak  = charcoal$peak,
    pos   = char_thresh$pos[, ncol(char_thresh$pos)],
    neg   = char_thresh$neg[, ncol(char_thresh$neg)]
  )

  base_val <- if (!is.null(peak_analysis$cPeak) &&
                  peak_analysis$cPeak == 2L) 1 else 0

  peak_rows   <- x[post$peakIn]
  screen_rows <- x[post$peakScreenIn]
  y_max_b <- max(df_b$peak, na.rm = TRUE)
  y_min_b <- min(df_b$neg,  na.rm = TRUE)

  y_label_b <- if (!is.null(peak_analysis$cPeak) &&
                   peak_analysis$cPeak == 2L) {
    expression(paste("C"["peak"], " (C"["interp"]*" / C"["bkg"]*")"))
  } else {
    expression(paste("C"["peak"], " (pieces cm"^{-2}*" yr"^{-1}*")"))
  }

  pb <- ggplot2::ggplot(df_b) +
    ggplot2::geom_col(ggplot2::aes(x = x, y = peak),
                      fill = "black", colour = "black", width = r) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = pos), colour = "red") +
    ggplot2::geom_line(ggplot2::aes(x = x, y = neg), colour = "red") +
    ggplot2::geom_hline(yintercept = base_val, colour = "black", linewidth = 0.4)

  if (length(screen_rows) > 0L) {
    df_screen <- data.frame(x = screen_rows,
                             y = rep(0.8 * y_max_b, length(screen_rows)))
    pb <- pb + ggplot2::geom_point(data = df_screen,
                                    ggplot2::aes(x = x, y = y),
                                    colour = "grey60", size = 1.5, shape = 16)
  }
  if (length(peak_rows) > 0L) {
    df_peaks <- data.frame(x = peak_rows,
                            y = rep(0.8 * y_max_b, length(peak_rows)))
    pb <- pb + ggplot2::geom_point(data = df_peaks,
                                    ggplot2::aes(x = x, y = y),
                                    colour = "black", size = 2, shape = 3,
                                    stroke = 1)
  }

  peak_type <- if (!is.null(peak_analysis$cPeak) &&
                   peak_analysis$cPeak == 2L) "ratio" else "residual"
  # Two-line title matching MATLAB: line break after "peaks"
  cpeak_label <- if (.has_ggtext()) {
    if (peak_type == "ratio") {
      .title("(b) C<sub>peak</sub> (C<sub>interpolated</sub> / C<sub>background</sub>), thresholds defining C<sub>noise</sub>, and peaks<br>identified (gray peaks fail to pass peak-magnitude test)")
    } else {
      .title("(b) C<sub>peak</sub> (C<sub>interpolated</sub> - C<sub>background</sub>), thresholds defining C<sub>noise</sub>, and peaks<br>identified (gray dots fail to pass peak-magnitude test)")
    }
  } else {
    if (peak_type == "ratio") {
      "(b) C_peak (C_interpolated / C_background), thresholds defining C_noise, and peaks\nidentified (gray peaks fail to pass peak-magnitude test)"
    } else {
      "(b) C_peak (C_interpolated - C_background), thresholds defining C_noise, and peaks\nidentified (gray dots fail to pass peak-magnitude test)"
    }
  }

  pb <- pb +
    .zone_lines(zone_div, y_min_b, y_max_b * 1.0) +
    .char_xscale(zone_div) +
    ggplot2::scale_y_continuous(limits = c(y_min_b, 1.1 * y_max_b),
                                 expand = c(0, 0)) +
    ggplot2::ylab(y_label_b) +
    ggplot2::labs(title = cpeak_label) +
    .char_theme()

  .combine_panels(list(pa, pb), heights = c(1.2, 1))
}

# =============================================================================
# char_plot_fire_history  --  Figure 7: peak magnitude / FRIs / fire frequency
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_fire_history <- function(out) {

  .require_ggplot2()

  charcoal      <- out$charcoal
  pretreatment  <- out$pretreatment
  peak_analysis <- out$peak_analysis
  post          <- out$post
  site          <- out$site

  zone_div  <- pretreatment$zoneDiv
  yr_interp <- pretreatment$yrInterp
  x         <- charcoal$ybpI
  r         <- yr_interp

  peak_rows   <- x[post$peakIn]
  screen_rows <- x[post$peakScreenIn]

  # ---------- Panel 1: Peak magnitude ---------------------------------------
  df_mag <- data.frame(x   = x,
                        mag = charcoal$peakMagnitude)
  y_max_mag <- max(df_mag$mag, na.rm = TRUE)
  if (y_max_mag == 0) y_max_mag <- 1  # avoid degenerate limits

  p1 <- ggplot2::ggplot(df_mag) +
    ggplot2::geom_col(ggplot2::aes(x = x, y = mag),
                      fill = "black", colour = "black", width = r)

  if (length(screen_rows) > 0L) {
    df_s <- data.frame(x = screen_rows,
                        y = rep(0.8 * y_max_mag, length(screen_rows)))
    p1 <- p1 + ggplot2::geom_point(data = df_s, ggplot2::aes(x = x, y = y),
                                    colour = "grey60", size = 1.5, shape = 16)
  }
  if (length(peak_rows) > 0L) {
    df_pk <- data.frame(x = peak_rows,
                         y = rep(0.8 * y_max_mag, length(peak_rows)))
    p1 <- p1 + ggplot2::geom_point(data = df_pk, ggplot2::aes(x = x, y = y),
                                    colour = "red", size = 2, shape = 3,
                                    stroke = 1.2)
  }

  p1 <- p1 +
    .zone_lines(zone_div, y_max_mag * 0.90, y_max_mag * 1.0) +
    .zone_labels(zone_div, y_max_mag * 0.95) +
    .char_xscale(zone_div) +
    ggplot2::scale_y_continuous(limits = c(0, 1.1 * y_max_mag), expand = c(0, 0)) +
    ggplot2::ylab(expression(paste("peak mag. (pieces cm"^{-2}*" peak"^{-1}*")"))) +
    ggplot2::labs(title = .title(paste0(
        "Peak magnitude, FRIs, and fire frequ.\n\n", site))) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +
    .char_theme()

  # ---------- Panel 2: FRIs through time + smoothed FRI ---------------------
  FRIyr  <- post$FRIyr
  FRI    <- post$FRI
  smFRIyr <- post$smFRIyr
  smFRI   <- post$smFRI
  smFRIci <- post$smFRIci   # [upper, lower] columns (MATLAB convention: [2.5%, 97.5%])
  alpha   <- post$alpha %||% 0.05

  has_fri <- !is.null(FRIyr) && length(FRIyr) > 1L &&
             !all(is.na(FRIyr)) && length(smFRI) > 2L

  if (has_fri) {
    y_max_fri <- 1.1 * max(FRI, na.rm = TRUE)
    df_fri  <- data.frame(x = FRIyr, y = FRI)
    df_sm   <- data.frame(x = smFRIyr, y = smFRI)
    # CI: MATLAB smFRIci column 1 = 2.5th quantile (lower), column 2 = 97.5th (upper)
    df_ci   <- data.frame(x    = smFRIyr,
                           y_lo = smFRIci[, 1L],
                           y_hi = smFRIci[, 2L])

    p2 <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = df_ci,
                            ggplot2::aes(x = x, ymin = y_lo, ymax = y_hi),
                            fill = "grey80", alpha = 0.8) +
      ggplot2::geom_point(data = df_fri, ggplot2::aes(x = x, y = y),
                           shape = 22, fill = "grey75", colour = "black",
                           size  = 2) +
      ggplot2::geom_line(data = df_sm, ggplot2::aes(x = x, y = y),
                          colour = "black", linewidth = 1) +
      .zone_lines(zone_div, 0, y_max_fri) +
      .char_xscale(zone_div) +
      ggplot2::scale_y_continuous(limits = c(0, y_max_fri), expand = c(0, 0)) +
      ggplot2::ylab(
        if (.has_ggtext()) {
          paste0("FRI (yr fire<sup>-1</sup>)<br>",
                 peak_analysis$peakFrequ, "-yr mean<br>",
                 round((1 - alpha) * 100), "% CI")
        } else {
          paste0("FRI (yr fire^-1)\n",
                 peak_analysis$peakFrequ, "-yr mean\n",
                 round((1 - alpha) * 100), "% CI")
        }) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x  = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) +
      .char_theme()
  } else {
    p2 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                         label = "Fewer than 3 FRIs - FRI plot not available",
                         hjust = 0.5) +
      ggplot2::theme_void()
  }

  # ---------- Panel 3: Smoothed fire frequency ------------------------------
  df_ff <- data.frame(x  = x,
                       ff = charcoal$smoothedFireFrequ)
  y_max_ff <- max(1.1 * max(df_ff$ff, na.rm = TRUE), .Machine$double.eps)

  p3 <- ggplot2::ggplot(df_ff) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = ff),
                        colour = "black", linewidth = 0.8) +
    .zone_lines(zone_div, 0, y_max_ff) +
    .char_xscale(zone_div) +
    ggplot2::scale_y_continuous(limits = c(0, y_max_ff), expand = c(0, 0)) +
    ggplot2::ylab(
      if (.has_ggtext()) {
        paste0("fire freq.<br>(fires ",
               peak_analysis$peakFrequ, " yr<sup>-1</sup>)")
      } else {
        paste0("fire freq.\n(fires ",
               peak_analysis$peakFrequ, " yr^-1)")
      }) +
    .char_theme()

  .combine_panels(list(p1, p2, p3), heights = c(1.2, 1.2, 1))
}

# =============================================================================
# char_plot_cumulative  --  Figure 5: cumulative peaks through time
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_cumulative <- function(out) {

  .require_ggplot2()

  charcoal     <- out$charcoal
  pretreatment <- out$pretreatment
  post         <- out$post
  site         <- out$site

  zone_div <- pretreatment$zoneDiv
  x_all    <- charcoal$ybpI
  peaks_col <- post$CharcoalCharPeaks[, ncol(post$CharcoalCharPeaks)]

  peak_idx <- which(peaks_col > 0)

  if (length(peak_idx) == 0L) {
    message("char_plot_cumulative: no peaks to plot.")
    return(invisible(NULL))
  }

  # Cumulative count ordered young-to-old (x reversed), so the youngest
  # peak gets the highest cumulative number -- matching MATLAB's flipud(cumsum).
  x_plot <- x_all[peak_idx]
  y_plot <- rev(seq_along(peak_idx))   # equals flipud(cumsum(1:n))

  df <- data.frame(x = x_plot, y = y_plot)

  y_max <- max(y_plot)

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(colour = "black", size = 1.5) +
    .zone_lines(zone_div, 0, y_max) +
    .char_xscale(zone_div) +
    ggplot2::scale_y_continuous(limits = c(0, y_max * 1.05), expand = c(0, 0)) +
    ggplot2::ylab("cumulative number of peaks") +
    ggplot2::labs(
      title = .title(paste0(site, ": Cumulative fires as a function of time"))
    ) +
    .char_theme() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey90",
                                                             linewidth = 0.4))
}

# =============================================================================
# char_plot_fri  --  Figure 6: FRI histograms with Weibull fits by zone
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_fri <- function(out) {

  .require_ggplot2()

  pretreatment <- out$pretreatment
  post         <- out$post
  site         <- out$site

  zone_div <- pretreatment$zoneDiv
  n_zones  <- length(zone_div) - 1L
  zones    <- post$zone   # list of per-zone results from char_post_process

  zone_labels <- paste0("Zone ", seq_len(n_zones))

  # KS pass criterion: p > 0.10 if n < 30, p > 0.05 if n >= 30
  .ks_pass <- function(p, n) {
    if (is.na(p) || is.null(p)) return(FALSE)
    if (n < 30L) p > 0.10 else p > 0.05
  }

  x_lim <- c(0, 800)

  panel_list <- vector("list", n_zones)

  for (z in seq_len(n_zones)) {
    zd <- zones[[z]]

    # ---- Empty / insufficient data ----------------------------------------
    if (is.null(zd) || length(zd$FRI) <= 1L) {
      panel_list[[z]] <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                           label = "insufficient data",
                           hjust = 0.5, vjust = 0.5) +
        ggplot2::labs(title = zone_labels[z]) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = .title_el(face = "bold", hjust = 0.5))
      next
    }

    fri_z    <- zd$FRI
    bin_c    <- zd$bin_centers               # 30, 50, ..., 990
    freq     <- zd$freq                      # counts per bin (length 49)
    prop     <- freq / sum(freq)             # normalise to proportions
    p_ks     <- zd$p_ks
    n_fri    <- length(fri_z)
    m_fri    <- zd$mean_mFRI
    m_ci     <- zd$mean_mFRI_ci              # [2.5%, 97.5%]
    wbl_b    <- zd$wbl_scale                 # scale
    wbl_c    <- zd$wbl_shape                 # shape
    wbl_b_ci <- zd$wbl_scale_ci             # [2.5%, 97.5%]
    wbl_c_ci <- zd$wbl_shape_ci

    df_hist <- data.frame(x = bin_c, prop = prop)

    p <- ggplot2::ggplot(df_hist) +
      ggplot2::geom_col(ggplot2::aes(x = x, y = prop),
                         fill = "grey75", colour = "grey50", width = 20)

    # ---- Weibull overlay (when fit converged and n > 4) --------------------
    # The KS goodness-of-fit test is computed and shown, but the Weibull curve
    # is displayed whenever the model converged (n > 4 and wbl parameters are
    # available).  MATLAB gates on passKS; here we always draw the curve and
    # report the KS p-value in the annotation so the user has full information.
    # Rationale: floating-point differences in the GMM between R and MATLAB
    # can shift the local threshold slightly, changing the peak count and
    # therefore the KS p-value for borderline datasets (e.g. CH10).  Hiding
    # the Weibull in those cases obscures a legitimate and visually good fit.
    pass_ks     <- n_fri > 4L && .ks_pass(p_ks, n_fri)
    show_weibull <- n_fri > 4L && !is.null(wbl_b) && !is.null(wbl_c)
    annot_lines <- character(0L)

    if (show_weibull) {
      wbl_x   <- seq(1, 1000, by = 1)
      wbl_pdf <- stats::dweibull(wbl_x, shape = wbl_c, scale = wbl_b)
      # Scale PDF to match histogram proportions (bin width 20)
      wbl_pdf_scaled <- wbl_pdf * 20
      df_wbl <- data.frame(x = wbl_x, y = wbl_pdf_scaled)
      p <- p + ggplot2::geom_line(data = df_wbl,
                                    ggplot2::aes(x = x, y = y),
                                    colour = "black", linewidth = 1)

      # KS p-value line: show passing result normally; flag a borderline/failing
      # result with "(n.s.)" so the user can see the statistical quality.
      ks_label <- if (is.na(p_ks)) {
        "KS p = NA"
      } else if (pass_ks) {
        sprintf("KS p = %.3f", p_ks)
      } else {
        sprintf("KS p = %.3f (n.s.)", p_ks)
      }

      annot_lines <- c(
        sprintf("Wbl *b* = %d (%d-%d)",
                round(wbl_b), round(wbl_b_ci[1L]), round(wbl_b_ci[2L])),
        sprintf("Wbl *c* = %.2f (%.2f-%.2f)",
                wbl_c, wbl_c_ci[1L], wbl_c_ci[2L]),
        sprintf("mFRI = %d (%d-%d)",
                round(m_fri), round(m_ci[1L]), round(m_ci[2L])),
        sprintf("N = %d", n_fri),
        ks_label
      )
    } else {
      annot_lines <- sprintf("N = %d", n_fri)
    }

    # ---- Text annotation --------------------------------------------------
    # geom_richtext renders HTML: use <br> for line breaks.
    # annotate("text") renders plain text: use \n for line breaks.
    y_max_hist <- max(prop, na.rm = TRUE)

    if (.has_ggtext()) {
      annot_html <- paste(annot_lines, collapse = "<br>")
      # annot_html must live in the data frame, NOT referenced via aes().
      # If it is referenced by name in aes(), ggplot2 lazy-evaluates it and
      # picks up whatever value the variable holds when the panel is finally
      # rendered (after the loop finishes) -- i.e. always the last zone.
      # Putting it in the data frame captures the value immediately.
      p <- p + ggtext::geom_richtext(
        data  = data.frame(x = x_lim[2L], y = y_max_hist, lbl = annot_html),
        ggplot2::aes(x = x, y = y, label = lbl),
        hjust = 1, vjust = 1, size = 2.8,
        fill  = NA, label.colour = NA, lineheight = 1.3
      )
    } else {
      annot_plain <- paste(gsub("\\*", "", annot_lines), collapse = "\n")
      p <- p + ggplot2::annotate(
        "text",
        x = x_lim[2L], y = y_max_hist,
        label = annot_plain,
        hjust = 1, vjust = 1, size = 2.8, lineheight = 1.3
      )
    }

    # ---- Axis and theme ---------------------------------------------------
    p <- p +
      ggplot2::scale_x_continuous(limits = x_lim,
                                   breaks = seq(0, x_lim[2L], by = 200)) +
      ggplot2::scale_y_continuous(limits = c(0, max(y_max_hist * 1.3, 0.35)),
                                   expand = c(0, 0)) +
      ggplot2::xlab("FRI (yr)") +
      ggplot2::ylab(if (z == n_zones) "proportion" else NULL) +
      ggplot2::labs(title = zone_labels[z]) +
      .char_theme() +
      ggplot2::theme(plot.title = .title_el(face = "bold", hjust = 0.5))

    panel_list[[z]] <- p
  }

  # Combine panels left-to-right (patchwork), or print individually.
  # Panels are reversed so Zone 1 is always on the far right, matching
  # MATLAB's inPlot = fliplr(1:nZones) subplot ordering.
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- Reduce(`+`, rev(panel_list)) +
      patchwork::plot_annotation(
        title = .title(paste0(site,
                               ": FRI distributions by zone with Weibull models")),
        theme = ggplot2::theme(plot.title = .title_el(face = "bold"))
      ) +
      patchwork::plot_layout(nrow = 1L)
    return(combined)
  }

  message("Install 'patchwork' for combined multi-panel layout.")
  invisible(panel_list[[n_zones]])
}

# =============================================================================
# char_plot_zones  --  Figure 8: between-zone raw CHAR comparisons
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_zones <- function(out) {

  .require_ggplot2()

  charcoal     <- out$charcoal
  pretreatment <- out$pretreatment
  site         <- out$site

  zone_div    <- pretreatment$zoneDiv
  n_zones     <- length(zone_div) - 1L
  zone_labels <- paste0("Zone ", seq_len(n_zones))

  if (is.null(charcoal$acc) || is.null(charcoal$ybp)) {
    message("char_plot_zones: raw CHAR (charcoal$acc / charcoal$ybp) not available.")
    return(invisible(NULL))
  }

  # Assign each raw sample to a zone
  zone_id <- cut(charcoal$ybp,
                  breaks = zone_div,
                  labels = zone_labels,
                  include.lowest = TRUE, right = FALSE)
  df_raw <- data.frame(
    acc  = charcoal$acc,
    ybp  = charcoal$ybp,
    zone = zone_id
  )
  df_raw <- df_raw[!is.na(df_raw$zone), ]

  # Zone colours: black + grey scale cycling for additional zones
  zone_cols <- c("black", "grey50", "grey70",
                  "#1b7837", "#762a83", "#d6604d",
                  "#f46d43", "#006837")[seq_len(n_zones)]

  # ---- Left panel: empirical CDFs ------------------------------------------
  x_max_cdf <- max(df_raw$acc, na.rm = TRUE)

  p_cdf <- ggplot2::ggplot(df_raw,
                             ggplot2::aes(x = acc, colour = zone,
                                          group = zone)) +
    ggplot2::stat_ecdf(linewidth = 1) +
    ggplot2::scale_colour_manual(values = zone_cols,
                                  name   = NULL,
                                  labels = zone_labels) +
    ggplot2::scale_x_continuous(limits = c(0, x_max_cdf), expand = c(0, 0)) +
    ggplot2::xlab(expression(paste("CHAR (pieces cm"^{-2}*" yr"^{-1}*")"))) +
    ggplot2::ylab("cumulative proportion") +
    ggplot2::labs(title = .title("CDFs of zone-specific raw CHAR, with KS test results")) +
    .char_theme() +
    ggplot2::theme(legend.position        = c(0.85, 0.15),
                   legend.justification   = c("right", "bottom"),
                   legend.background      = ggplot2::element_rect(fill  = "white",
                                                                    colour = "grey80",
                                                                    linewidth = 0.3))

  # ---- KS test table (annotated on CDF plot) --------------------------------
  if (n_zones > 1L) {
    zone_levels <- levels(df_raw$zone)
    ks_mat <- matrix(NA_real_, n_zones, n_zones,
                      dimnames = list(zone_levels, zone_levels))
    for (i in seq_len(n_zones - 1L)) {
      for (j in (i + 1L):n_zones) {
        xi <- df_raw$acc[df_raw$zone == zone_levels[i]]
        xj <- df_raw$acc[df_raw$zone == zone_levels[j]]
        if (length(xi) > 1L && length(xj) > 1L) {
          res <- suppressWarnings(stats::ks.test(xi, xj))
          ks_mat[i, j] <- round(res$p.value, 3L)
        }
      }
    }
    # Build annotation string
    ks_lines <- "KS p-values:"
    for (i in seq_len(n_zones - 1L)) {
      for (j in (i + 1L):n_zones) {
        if (!is.na(ks_mat[i, j])) {
          ks_lines <- c(ks_lines,
                         sprintf("  %s vs %s: p = %.3f",
                                 zone_levels[i], zone_levels[j],
                                 ks_mat[i, j]))
        }
      }
    }
    p_cdf <- p_cdf +
      ggplot2::annotate("text",
                         x    = max(df_raw$acc, na.rm = TRUE) * 0.5,
                         y    = 0.45,
                         label = paste(ks_lines, collapse = "\n"),
                         hjust = 0, vjust = 1, size = 3)
  }

  # ---- Right panel: box plots ----------------------------------------------
  p_box <- ggplot2::ggplot(df_raw,
                             ggplot2::aes(x = zone, y = acc,
                                          fill = zone, colour = zone)) +
    ggplot2::geom_boxplot(
      outlier.shape  = 16, outlier.size = 1,
      coef           = 1.5   # standard whiskers at 1.5 * IQR
    ) +
    ggplot2::scale_fill_manual(
      values = ggplot2::alpha(zone_cols, 0.4),
      guide  = "none") +
    ggplot2::scale_colour_manual(values = zone_cols, guide = "none") +
    ggplot2::xlab("zone") +
    ggplot2::ylab(expression(paste("CHAR (pieces cm"^{-2}*" yr"^{-1}*")"))) +
    ggplot2::labs(title = .title("Box plots of raw CHAR per zone")) +
    .char_theme()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- (p_cdf | p_box) +
      patchwork::plot_annotation(
        title = .title(paste0(site,
                               ": Between-zone comparisons of raw CHAR distributions")),
        theme = ggplot2::theme(plot.title = .title_el(face = "bold"))
      ) +
      patchwork::plot_layout(widths = c(2, 1))
    return(combined)
  }

  message("Install 'patchwork' for combined multi-panel layout.")
  print(p_cdf)
  invisible(p_box)
}

# =============================================================================
# char_plot_raw  --  Figure 1: C_raw / C_resampled / C_background options
# Only produced when out$results$allFigures == 1.
# Mirrors MATLAB CharPretreatment.m (subplot 1) + CharSmooth.m (subplot 2).
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_raw <- function(out) {

  .require_ggplot2()

  charcoal  <- out$charcoal
  smoothing <- out$smoothing
  pretreat  <- out$pretreatment
  site      <- out$site %||% ""
  yr_interp <- pretreat$yrInterp
  has_gt    <- .has_ggtext()

  # Helper: subscript label, rendered in HTML when ggtext is available
  sub_lbl <- function(base, sub)
    if (has_gt) paste0(base, "<sub>", sub, "</sub>") else paste0(base, "_", sub)
  # Helper: axis title element that renders markdown/HTML
  md_axis <- function()
    if (has_gt) ggtext::element_markdown() else ggplot2::element_text()

  char_lbl <- sub_lbl("CHAR (# cm", "-2")
  char_lbl <- if (has_gt)
    "CHAR (# cm<sup>-2</sup> yr<sup>-1</sup>)" else "CHAR (# cm^-2 yr^-1)"

  x_lim <- c(max(charcoal$ybp, na.rm = TRUE) + 100,
              min(charcoal$ybp, na.rm = TRUE))  # reversed (old -> young right)

  # Legend labels
  lbl_raw  <- sub_lbl("C", "raw")
  lbl_int  <- paste0(sub_lbl("C", "interpolated"), ": ", yr_interp, " yr")

  # -- Subplot 1: C_raw (bars) + C_interpolated (step) + legend --------------
  df_raw <- data.frame(age = charcoal$ybp, acc = charcoal$acc)
  df_int <- data.frame(age = charcoal$ybpI - 0.5 * yr_interp,
                       acc = charcoal$accI)

  p1 <- ggplot2::ggplot() +
    ggplot2::geom_bar(data = df_raw,
                      ggplot2::aes(x = age, y = acc, fill = "raw"),
                      stat = "identity",
                      width = diff(range(charcoal$ybp)) / length(charcoal$ybp),
                      colour = NA) +
    ggplot2::geom_step(data = df_int,
                       ggplot2::aes(x = age, y = acc, colour = "int"),
                       linewidth = 0.8) +
    ggplot2::scale_fill_manual(
      name   = NULL,
      values = c("raw" = "grey50"),
      labels = c("raw" = lbl_raw),
      guide  = ggplot2::guide_legend(
        override.aes = list(colour = NA, fill = "grey50",
                            linetype = 0, size = 4))
    ) +
    ggplot2::scale_colour_manual(
      name   = NULL,
      values = c("int" = "black"),
      labels = c("int" = lbl_int),
      guide  = ggplot2::guide_legend(
        override.aes = list(fill = NA, linetype = 1, linewidth = 0.8))
    ) +
    ggplot2::scale_x_reverse(limits = x_lim) +
    ggplot2::coord_cartesian(
      ylim = c(0, stats::quantile(charcoal$acc, 0.99, na.rm = TRUE))) +
    .char_theme() +
    ggplot2::theme(
      axis.title.y = md_axis(),
      legend.text  = md_axis()
    ) +
    ggplot2::labs(
      title = .title(paste0(site, " (a) ", sub_lbl("C", "raw"), " and ",
                            sub_lbl("C", "interpolated"), ": ", yr_interp, " yr")),
      x = "time (cal. yr BP)",
      y = char_lbl
    )

  # -- Subplot 2: C_interpolated (bars) + all 5 smoothing curves -------------
  smooth_names <- c("Lowess", "Robust Lowess",
                    "Moving Average", "Moving Median", "Moving Mode")
  all_curves   <- charcoal$accIS_all   # N x 5 matrix

  df_int2 <- data.frame(age = charcoal$ybpI, acc = charcoal$accI)

  df_smooth <- do.call(rbind, lapply(seq_len(5L), function(k) {
    data.frame(age    = charcoal$ybpI,
               smooth = all_curves[, k],
               method = smooth_names[k])
  }))
  df_smooth$method <- factor(df_smooth$method, levels = smooth_names)

  selected_name <- smooth_names[smoothing$method]
  df_sel  <- df_smooth[df_smooth$method == selected_name, ]
  df_rest <- df_smooth[df_smooth$method != selected_name, ]

  # Build a named colour vector: selected method = black, others = Set1 palette
  set1_4 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  non_sel_names <- smooth_names[smooth_names != selected_name]
  col_map <- stats::setNames(
    c("black", set1_4[seq_along(non_sel_names)]),
    c(selected_name, non_sel_names)
  )
  # Legend linewidth override in factor-level order (selected = thick)
  lw_override <- ifelse(smooth_names == selected_name, 1.2, 0.6)

  p2_title <- paste0(
    "(b) ", sub_lbl("C", "interpolated"),
    " and options for a ", smoothing$yr, " yr ",
    sub_lbl("C", "background"),
    " (selected: ", selected_name, ")"
  )

  p2 <- ggplot2::ggplot() +
    ggplot2::geom_bar(data = df_int2,
                      ggplot2::aes(x = age, y = acc),
                      stat = "identity",
                      width = diff(range(charcoal$ybpI, na.rm = TRUE)) /
                        length(charcoal$ybpI),
                      fill = "grey50", colour = NA) +
    ggplot2::geom_line(data = df_rest,
                       ggplot2::aes(x = age, y = smooth, colour = method),
                       linewidth = 0.6, alpha = 0.7) +
    ggplot2::geom_line(data = df_sel,
                       ggplot2::aes(x = age, y = smooth, colour = method),
                       linewidth = 1.2) +
    ggplot2::scale_x_reverse(limits = x_lim) +
    ggplot2::coord_cartesian(ylim = c(0, max(charcoal$accI, na.rm = TRUE))) +
    ggplot2::scale_colour_manual(
      values = col_map,
      name   = "Method",
      guide  = ggplot2::guide_legend(
        override.aes = list(linewidth = lw_override)
      )
    ) +
    .char_theme() +
    ggplot2::theme(axis.title.y = md_axis()) +
    ggplot2::labs(
      title = .title(p2_title),
      x     = "time (cal. yr BP)",
      y     = char_lbl
    )

  fig <- patchwork::wrap_plots(p1, p2, ncol = 1L)
  print(fig)
  invisible(fig)
}


# =============================================================================
# char_plot_thresh_diag  --  Figure 2: threshold determination diagnostics
# Only produced when out$results$allFigures == 1.
# Global threshold: single panel histogram + noise PDF + threshold lines.
# Local threshold:  5x5 grid of per-sample window distributions.
# Mirrors MATLAB CharThreshGlobal.m (lines 201-268) and
#         CharThreshLocal.m (lines 225-291).
#
# NOTE -- R vs. MATLAB difference (intentional, documented):
#   When threshMethod == 3 (Gaussian mixture model), the MATLAB version draws
#   both the noise and signal Gaussian components plus their weighted mixture
#   on each subplot.  The R version plots only the noise component (the
#   component with the smaller mean), which is the distribution used to define
#   the threshold.  This is an intentional simplification: the noise component
#   is the only part of the GMM that directly determines the threshold, and
#   plotting only it avoids visual confusion in windows where the two components
#   have very similar means.  This difference is noted in the package
#   documentation and release notes.
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_thresh_diag <- function(out) {

  .require_ggplot2()

  charcoal     <- out$charcoal
  char_thresh  <- out$char_thresh
  peak_analysis <- out$peak_analysis
  pretreat     <- out$pretreatment
  site         <- out$site %||% ""

  thresh_type <- peak_analysis$threshType  # 1 = global, 2 = local

  # ============================================================
  # GLOBAL THRESHOLD (Figure 2, single panel)
  # ============================================================
  if (thresh_type == 1L) {

    possible  <- char_thresh$possible
    noise_pdf <- char_thresh$noise_pdf
    thresh_pos <- char_thresh$pos[1L, ]   # 4 candidate levels

    # Histogram of C_peak proportions at the candidate bins
    bin_w  <- mean(diff(possible))
    counts <- graphics::hist(charcoal$peak, breaks = possible,
                             plot = FALSE)$counts
    prop   <- counts / sum(counts)
    df_hist <- data.frame(x = possible[-length(possible)] + bin_w / 2,
                          y = prop)

    p <- ggplot2::ggplot(df_hist, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_col(fill = "grey50", colour = "grey50", width = bin_w)

    # Noise PDF overlay (when data-defined threshold)
    if (!is.null(noise_pdf) && length(noise_pdf) == length(possible) &&
        !identical(noise_pdf, -99)) {
      df_pdf <- data.frame(x = possible,
                           y = noise_pdf * bin_w)
      p <- p + ggplot2::geom_line(data = df_pdf,
                                   ggplot2::aes(x = x, y = y),
                                   colour = "black", linewidth = 1.2)
    }

    y_max <- max(prop, na.rm = TRUE)

    # Dashed lines for all 4 candidate thresholds
    for (tv in thresh_pos) {
      p <- p + ggplot2::geom_vline(xintercept = tv,
                                   linetype = "dashed", colour = "black")
    }
    # Solid red line for selected (4th) threshold
    p <- p + ggplot2::geom_vline(xintercept = thresh_pos[4L],
                                 colour = "red", linewidth = 1.0)

    # Annotation: threshold value and mean SNI
    xannot <- thresh_pos[4L]
    sni_mean <- mean(char_thresh$SNI, na.rm = TRUE)
    p <- p +
      ggplot2::annotate("text", x = xannot * 1.5, y = y_max * 0.75,
                        label = paste0("Threshold = ",
                                       round(xannot, 4)),
                        hjust = 0, size = 3) +
      ggplot2::annotate("text", x = xannot * 1.5, y = y_max * 0.60,
                        label = paste0("SNI = ", round(sni_mean, 2)),
                        hjust = 0, size = 3)

    p <- p +
      .char_theme() +
      ggplot2::labs(
        title = paste0(site, ": ", pretreat$zoneDiv[1L], " to ",
                       pretreat$zoneDiv[length(pretreat$zoneDiv)],
                       " cal. yr BP"),
        x     = "peak CHAR (# cm^-2 yr^-1)",
        y     = "proportion or scaled density"
      )

    print(p)
    return(invisible(p))
  }

  # ============================================================
  # LOCAL THRESHOLD (Figure 2, 5x5 diagnostic grid)
  # ============================================================
  diag <- char_thresh$diag
  if (is.null(diag) || length(diag) == 0L) {
    message("char_plot_thresh_diag: no local threshold diagnostic data available.")
    return(invisible(NULL))
  }

  n_panels <- length(diag)
  n_col    <- 5L
  n_row    <- ceiling(n_panels / n_col)

  panels <- lapply(seq_len(n_panels), function(k) {
    d   <- diag[[k]]
    X   <- d$X
    bin_w_k <- diff(range(X, na.rm = TRUE)) / 50
    if (bin_w_k <= 0) bin_w_k <- 0.01
    breaks_k <- seq(min(X, na.rm = TRUE),
                    max(X, na.rm = TRUE) + bin_w_k,
                    by = bin_w_k)
    h <- graphics::hist(X, breaks = breaks_k, plot = FALSE)
    df_h <- data.frame(x = h$mids, y = h$counts / sum(h$counts))

    x_seq <- seq(min(X, na.rm = TRUE), max(X, na.rm = TRUE),
                 length.out = 200L)

    # Base panel
    pk <- ggplot2::ggplot(df_h, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_col(fill = "grey75", colour = "grey75",
                        width = bin_w_k)

    # PDF overlay: two-component GMM or single Gaussian
    if (peak_analysis$threshMethod == 3L && d$prop2 > 0) {
      df_pdf1 <- data.frame(
        x = x_seq,
        y = stats::dnorm(x_seq, d$mu1, d$sig1) * bin_w_k * d$prop1)
      df_pdf2 <- data.frame(
        x = x_seq,
        y = stats::dnorm(x_seq, d$mu2, d$sig2) * bin_w_k * d$prop2)
      df_mix  <- data.frame(
        x = x_seq,
        y = df_pdf1$y + df_pdf2$y)
      pk <- pk +
        ggplot2::geom_line(data = df_pdf1,
                           ggplot2::aes(x = x, y = y),
                           colour = "black", linewidth = 0.6) +
        ggplot2::geom_line(data = df_pdf2,
                           ggplot2::aes(x = x, y = y),
                           colour = "black", linewidth = 0.6) +
        ggplot2::geom_line(data = df_mix,
                           ggplot2::aes(x = x, y = y),
                           colour = "blue", linewidth = 0.8)
    } else {
      df_pdf <- data.frame(
        x = x_seq,
        y = stats::dnorm(x_seq, d$mu1, d$sig1) * bin_w_k)
      pk <- pk +
        ggplot2::geom_line(data = df_pdf,
                           ggplot2::aes(x = x, y = y),
                           colour = "black", linewidth = 0.9)
    }

    # Threshold line
    pk <- pk +
      ggplot2::geom_vline(xintercept = d$t_pos,
                          colour = "red", linewidth = 0.6) +
      ggplot2::coord_cartesian(ylim = c(0, 0.25)) +
      .char_theme() +
      ggplot2::theme(
        axis.text    = ggplot2::element_text(size = 6),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        plot.title   = .title_el(size = 7, face = "plain", hjust = 0)
      )

    # Shared axis labels: mirror MATLAB CharThreshLocal.m lines 281-286,
    # which add ylabel to plotIn==11 (middle-left) and xlabel to plotIn==23
    # (bottom-center) only.  Fall back to k==n_panels-2 for xlabel when
    # fewer than 25 panels are plotted (record shorter than 25 windows).
    if (k == 11L) {
      pk <- pk +
        ggplot2::labs(y = "proportion or density (scaled)") +
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = 7, angle = 90,
                                                            vjust = 0.5))
    }
    x_lbl_k <- if (n_panels >= 23L) 23L else max(n_panels - 2L, 1L)
    if (k == x_lbl_k) {
      x_txt <- if (.has_ggtext()) {
        "CHAR (# cm<sup>-2</sup> yr<sup>-1</sup>)"
      } else {
        "CHAR (# cm^-2 yr^-1)"
      }
      pk <- pk +
        ggplot2::labs(x = x_txt) +
        ggplot2::theme(
          axis.title.x = if (.has_ggtext())
            ggtext::element_markdown(size = 7)
          else
            ggplot2::element_text(size = 7)
        )
    }

    # Title and annotation
    yr_label <- round(d$yr_bp)
    panel_title <- if (k == 1L) paste0(site, ": ", yr_label, " yr BP")
                   else paste0(yr_label, " yr BP")
    sni_label <- paste0("SNI = ",  round(d$sni, 2), "\n",
                        "KS p = ", round(d$gof, 2), "\n",
                        "t = ",    round(d$t_pos, 3))
    pk <- pk +
      ggplot2::labs(title = panel_title) +
      ggplot2::annotate("text",
                        x    = max(X, na.rm = TRUE),
                        y    = 0.24,
                        label = sni_label,
                        hjust = 1, vjust = 1, size = 2.5)

    pk
  })

  fig <- patchwork::wrap_plots(panels, ncol = n_col)
  print(fig)
  invisible(fig)
}


# =============================================================================
# char_plot_sni  --  Figure 4: Sensitivity to alternative thresholds and SNI
# Mirrors MATLAB char_plot_sni_ThresholdSNI.m
#
# Layout (mimics MATLAB 3x5 subplot grid):
#   Panel (a): C_interp bar chart with C_back and final threshold; peaks at
#              each threshold marked as dots, selected threshold as +.
#   Panel (b): Global -- C_peak histogram, noise PDF, peak-count curve.
#              Local  -- mean FRI +/- 95% CI by zone for each threshold.
#   Panel (c): SNI time series with dashed reference at SNI = 3.
#   Panel (d): SNI distribution boxplot.
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_sni <- function(out) {

  .require_ggplot2()

  charcoal     <- out$charcoal
  char_thresh  <- out$char_thresh
  peak_analysis <- out$peak_analysis
  pretreatment <- out$pretreatment
  post         <- out$post
  site         <- out$site

  zone_div  <- pretreatment$zoneDiv
  transform <- pretreatment$transform
  cPeak     <- peak_analysis$cPeak
  thresh_type <- peak_analysis$threshType   # 1 = global, 2 = local

  x   <- charcoal$ybpI
  y   <- charcoal$accI
  y2  <- charcoal$accIS
  y3  <- char_thresh$pos[, ncol(char_thresh$pos)]   # final threshold (+)
  y4  <- if (is.matrix(char_thresh$neg)) char_thresh$neg[, 1L]
         else rep(char_thresh$neg[1L], length(x))

  # Positive and negative threshold lines on the C_interp scale
  t_pos_line <- if (cPeak == 1L) y2 + y3 else y2 * y3
  t_neg_line <- if (cPeak == 1L) y2 + y4 else y2 * y4

  # CharcoalCharPeaks: N x T matrix (T = number of thresholds)
  ccp        <- post$CharcoalCharPeaks
  T_thresh   <- ncol(ccp)
  y_max      <- max(y, na.rm = TRUE)

  # y-levels for peak markers (mirroring MATLAB's 0.78 / 0.85 / 0.92)
  lev <- c(0.78, 0.85, 0.92)
  if (T_thresh < 3L) lev <- lev[(4L - T_thresh):3L]

  # Identify which column is the "selected" (Final) threshold (always last column)
  sel_col <- T_thresh

  # y-axis label: mirrors charYLabel(transform)
  y_lbl <- if (.has_ggtext()) {
    switch(as.character(transform),
           "1" = "log CHAR (# cm<sup>-2</sup> yr<sup>-1</sup>)",
           "2" = "ln CHAR (# cm<sup>-2</sup> yr<sup>-1</sup>)",
               "CHAR (# cm<sup>-2</sup> yr<sup>-1</sup>)")
  } else {
    switch(as.character(transform),
           "1" = "log CHAR (# cm^-2 yr^-1)",
           "2" = "ln CHAR (# cm^-2 yr^-1)",
               "CHAR (# cm^-2 yr^-1)")
  }
  y_lbl_axis <- if (.has_ggtext()) ggtext::element_markdown() else ggplot2::element_text()

  # Helper: x-axis ticks in 1000s of years
  x_ticks <- seq(0, max(zone_div, na.rm = TRUE), by = 1000)

  # ============================================================
  # PANEL (a): C_interp + C_back + threshold + peak markers
  # ============================================================
  df_bar  <- data.frame(x = x, y = y)
  df_back <- data.frame(x = x, y2 = y2)
  df_tpos <- data.frame(x = x, t = t_pos_line)
  df_tneg <- data.frame(x = x, t = t_neg_line)

  pa <- ggplot2::ggplot(df_bar, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_col(fill = "black", colour = "black", width = mean(diff(x), na.rm = TRUE)) +
    ggplot2::geom_line(data = df_back, ggplot2::aes(x = x, y = y2),
                       colour = "grey50", linewidth = 1.2) +
    ggplot2::geom_line(data = df_tpos, ggplot2::aes(x = x, y = t),
                       colour = "red", linewidth = 0.6) +
    ggplot2::geom_line(data = df_tneg, ggplot2::aes(x = x, y = t),
                       colour = "red", linewidth = 0.6)

  # Zone dividers
  if (length(zone_div) > 2L) {
    for (zd in zone_div[-c(1L, length(zone_div))]) {
      pa <- pa + ggplot2::annotate("segment",
                                   x = zd, xend = zd,
                                   y = y_max * 1.01, yend = y_max * 1.10,
                                   colour = "grey50", linewidth = 1.0)
    }
    for (z in seq_along(zone_div[-1L])) {
      mid_x <- mean(zone_div[z:(z + 1L)])
      pa <- pa + ggplot2::annotate("text",
                                   x = mid_x, y = y_max * 1.05,
                                   label = paste0("Zone ", z),
                                   hjust = 0.5, size = 2.5)
    }
  }

  # Peak markers: plot the first (T_thresh - 1) threshold columns as grey dots,
  # sorted by threshold value lowest → highest, mapped to y-levels 0.78 → 0.92
  # (bottom to top).  Then overlay black + at whichever y-level corresponds to
  # the Final threshold value, mirroring MATLAB char_plot_sni_ThresholdSNI.m.
  #
  # The 4th threshValue always duplicates one of the first three; the matching
  # y-level is found by comparing threshold values, not peak vectors.
  tv          <- peak_analysis$threshValues          # all T_thresh values
  non_sel_idx <- setdiff(seq_len(T_thresh), sel_col) # columns 1:(T_thresh-1)
  non_sel_tv  <- tv[non_sel_idx]                     # their threshold values
  show_cols   <- utils::head(non_sel_idx[order(non_sel_tv)], 3L) # sorted asc, ≤3

  for (rank in seq_along(show_cols)) {
    j      <- show_cols[rank]
    lev_j  <- lev[rank]
    pk_idx <- which(ccp[, j] > 0)
    if (length(pk_idx) > 0L) {
      df_pk <- data.frame(x = x[pk_idx], y = y_max * lev_j)
      pa <- pa +
        ggplot2::geom_point(data = df_pk, ggplot2::aes(x = x, y = y),
                            shape = 16L, size = 1, colour = "grey50")
    }
  }

  # Determine mIndex: find the rank (in the sorted non-selected columns) of the
  # column whose threshValue matches the Final threshValue, then read its y-level.
  sel_tv        <- tv[sel_col]
  sorted_tv     <- non_sel_tv[order(non_sel_tv)][seq_along(show_cols)]
  match_rank    <- which(sorted_tv == sel_tv)
  if (length(match_rank) == 0L) match_rank <- ceiling(length(show_cols) / 2L)
  mIndex        <- lev[match_rank[1L]]

  # Selected (final) peaks: white layer erases the grey dot, then black +
  pk_idx_sel <- which(ccp[, sel_col] > 0)
  if (length(pk_idx_sel) > 0L) {
    df_pk_sel <- data.frame(x = x[pk_idx_sel], y = y_max * mIndex)
    pa <- pa +
      ggplot2::geom_point(data = df_pk_sel, ggplot2::aes(x = x, y = y),
                          shape = 3L, size = 2.5, colour = "white", stroke = 1.5) +
      ggplot2::geom_point(data = df_pk_sel, ggplot2::aes(x = x, y = y),
                          shape = 3L, size = 2.5, colour = "black", stroke = 0.9)
  }

  pa <- pa +
    ggplot2::scale_x_reverse(limits = c(max(zone_div), min(zone_div)),
                              breaks = x_ticks,
                              labels = x_ticks / 1000,
                              expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(limits = c(0, y_max * 1.15),
                                expand = ggplot2::expansion(0)) +
    .char_theme() +
    ggplot2::theme(
      axis.title.y = y_lbl_axis,
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = if (.has_ggtext())
        paste0(site, " (a) C<sub>interpolated</sub>, C<sub>background</sub>, and peak ID (+)")
      else
        paste0(site, " (a) C_interpolated, C_background, and peak ID (+)"),
      y = y_lbl
    )

  if (.has_ggtext()) {
    pa <- pa + ggplot2::theme(plot.title = ggtext::element_markdown())
  }

  # ============================================================
  # PANEL (b): Threshold sensitivity
  # ============================================================
  if (thresh_type == 1L) {
    # ---- Global: C_peak histogram + noise PDF + peak-count curve ------------
    cpk      <- charcoal$peak
    possible <- char_thresh$possible
    noise_pdf <- char_thresh$noise_pdf

    # Histogram (probability)
    bk   <- possible
    cnts <- graphics::hist(cpk, breaks = bk, plot = FALSE)
    n    <- cnts$counts / sum(cnts$counts)
    xh_c <- cnts$mids

    # Noise PDF scaled to histogram bin width (method > 1 only)
    bw <- mean(diff(bk), na.rm = TRUE)
    has_noise <- !is.null(noise_pdf) && length(noise_pdf) > 1L &&
                   !all(noise_pdf == -99)

    # Right axis: number of peaks as function of threshold
    y3tot <- vapply(possible, function(tv) sum(cpk > tv, na.rm = TRUE),
                    integer(1L))

    # Filter x range for peak-count curve (mirrors MATLAB cPeak logic)
    if (cPeak == 1L) {
      xplot_thresh <- xh_c[xh_c > 0]
    } else {
      xplot_thresh <- xh_c[xh_c >= 1]
    }
    yplot_thresh <- stats::approx(possible, y3tot, xout = xplot_thresh,
                                   method = "linear", rule = 2)$y

    # Scale factor to overlay right axis on left axis
    scale_fac <- max(n, na.rm = TRUE) / max(y3tot[y3tot > 0], na.rm = TRUE)

    df_hist <- data.frame(x = xh_c, n = n)
    df_cnt  <- data.frame(x = xplot_thresh, y = yplot_thresh * scale_fac)

    pb <- ggplot2::ggplot(df_hist, ggplot2::aes(x = x, y = n)) +
      ggplot2::geom_col(fill = "grey75", colour = "grey75",
                        width = bw * 0.98)

    if (has_noise) {
      noise_pdf_sc <- stats::approx(possible, noise_pdf, xout = xh_c,
                                     method = "linear", rule = 2)$y * bw
      df_pdf <- data.frame(x = xh_c, y = noise_pdf_sc)
      pb <- pb + ggplot2::geom_line(data = df_pdf,
                                    ggplot2::aes(x = x, y = y),
                                    linetype = "dashed", linewidth = 1.2)
    }

    pb <- pb +
      ggplot2::geom_line(data = df_cnt, ggplot2::aes(x = x, y = y),
                         linewidth = 1.2)

    # Threshold vertical lines (grey dashed for all but last, solid for last)
    n_tv    <- ncol(char_thresh$pos)
    t_vals  <- char_thresh$pos[1L, ]   # first row (scalar for global)
    n_max   <- max(n, na.rm = TRUE)
    for (ti in seq_len(n_tv)) {
      if (ti < n_tv) {
        pb <- pb + ggplot2::geom_vline(xintercept = t_vals[ti],
                                       linetype = "dashed",
                                       colour = "grey75", linewidth = 0.8)
      } else {
        pb <- pb +
          ggplot2::geom_vline(xintercept = t_vals[ti],
                              colour = "black", linewidth = 0.8) +
          ggplot2::annotate("text", x = t_vals[ti], y = 0,
                            label = "<", angle = 90,
                            fontface = "bold", size = 3.5, vjust = 0.5)
        if (has_noise && !is.null(peak_analysis$threshValues)) {
          tv_pct <- peak_analysis$threshValues[n_tv]
          tv_val <- round(t_vals[ti], 4)
          ann_x  <- t_vals[ti] * 1.25
          ann_y  <- 0.9 * n_max
          ann_lbl <- paste0(tv_pct * 100, "th percentile\n= ", tv_val)
          pb <- pb + ggplot2::annotate("label", x = ann_x, y = ann_y,
                                       label = ann_lbl,
                                       hjust = 0, size = 2.5, fill = "white",
                                       label.size = 0)
        }
      }
    }

    x_lim_b <- c(min(cpk, na.rm = TRUE),
                  0.75 * max(cpk, na.rm = TRUE))
    pb <- pb +
      ggplot2::scale_x_continuous(limits = x_lim_b,
                                   expand = ggplot2::expansion(0)) +
      ggplot2::scale_y_continuous(
        name  = "relative frequency",
        limits = c(0, n_max * 1.01),
        expand = ggplot2::expansion(0),
        sec.axis = ggplot2::sec_axis(
          transform = ~ . / scale_fac,
          name  = "# of peaks identified"
        )
      ) +
      ggplot2::labs(
        x     = if (cPeak == 1L)
          "residual CHAR value"
        else
          "CHAR ratio",
        title = if (has_noise)
          "(b) C_peak dist., noise dist., and threshold"
        else
          "(b) C_peak distribution and threshold"
      ) +
      .char_theme()

  } else {
    # ---- Local: mean FRI +/- 95% CI per zone per threshold ------------------
    # Mirrors MATLAB: always plot columns 1-3 of CharcoalCharPeaks, then
    # find in2 = which of those columns shares its threshValue with the Final
    # (last) threshold.  in2 is the column to mark in black with "+".
    n_zones  <- length(zone_div) - 1L
    ccp_cols <- seq_len(min(T_thresh - 1L, 3L))   # columns 1-3 (non-Final)

    # in2: which plotted column matches the Final threshValue (mirrors MATLAB)
    sel_tv_local <- tv[sel_col]
    in2 <- which(tv[ccp_cols] == sel_tv_local)[1L]
    if (is.na(in2)) in2 <- length(ccp_cols)       # fallback: last plotted col

    zone_data <- lapply(seq_len(n_zones), function(z) {
      lapply(ccp_cols, function(j) {
        pk_yrs <- x[ccp[, j] > 0]
        pk_yrs <- pk_yrs[pk_yrs >= zone_div[z] & pk_yrs < zone_div[z + 1L]]
        fri_vals <- if (length(pk_yrs) > 1L) diff(pk_yrs) else NA_real_
        fri_vals <- fri_vals[fri_vals > 0 & !is.na(fri_vals)]
        mfri <- mean(fri_vals, na.rm = TRUE)
        ci   <- if (length(fri_vals) > 1L) {
          bm <- replicate(1000L, mean(sample(fri_vals, length(fri_vals),
                                             replace = TRUE)))
          stats::quantile(bm, c(0.025, 0.975))
        } else {
          c(NA_real_, NA_real_)
        }
        data.frame(zone = z, thresh = j,
                   mfri = mfri, lo = ci[1L], hi = ci[2L])
      })
    })

    df_zfri <- do.call(rbind, do.call(c, zone_data))
    df_zfri <- df_zfri[!is.na(df_zfri$mfri), ]

    # x positions: fliplr(1:n_zones) with offsets for threshold columns
    offsets <- c(-0.25, 0, 0.25)
    df_zfri$xpos <- (n_zones + 1L - df_zfri$zone) +
                    offsets[df_zfri$thresh]

    # Split on in2 (column index within the plotted 1-3 cols matching Final)
    df_nosel <- df_zfri[df_zfri$thresh != in2, ]
    df_sel   <- df_zfri[df_zfri$thresh == in2, ]

    pb <- ggplot2::ggplot() +
      # Non-selected thresholds: grey dots + grey error bars
      ggplot2::geom_errorbar(data = df_nosel,
                             ggplot2::aes(x = xpos, ymin = lo, ymax = hi),
                             colour = "grey50", width = 0.08) +
      ggplot2::geom_point(data = df_nosel,
                          ggplot2::aes(x = xpos, y = mfri),
                          shape = 16L, size = 2, colour = "grey50") +
      # Selected (final) threshold: black + + black error bars
      ggplot2::geom_errorbar(data = df_sel,
                             ggplot2::aes(x = xpos, ymin = lo, ymax = hi),
                             colour = "black", width = 0.08) +
      ggplot2::geom_point(data = df_sel,
                          ggplot2::aes(x = xpos, y = mfri),
                          shape = 16L, size = 2, colour = "black") +
      ggplot2::scale_x_continuous(
        breaks = seq_len(n_zones),
        labels = rev(seq_len(n_zones))
      ) +
      .char_theme() +
      ggplot2::labs(
        x     = "zone",
        y     = "zone-specific\nmean FRI (years)",
        title = "(b) Sensitivity to\nalternative thresholds"
      )
  }

  # ============================================================
  # PANEL (c): SNI time series
  # ============================================================
  sni_val <- char_thresh$SNI
  if (length(sni_val) == 1L) {
    sni_series <- rep(sni_val, length(x))
  } else {
    sni_series <- sni_val
  }
  df_sni <- data.frame(x = x, sni = sni_series)

  y_lim_c <- c(0, 10)
  y_tick_c <- if (y_lim_c[2L] < 20) seq(0, y_lim_c[2L], by = 2)
              else if (y_lim_c[2L] < 50) seq(0, y_lim_c[2L], by = 5)
              else seq(0, y_lim_c[2L], by = 10)

  pc <- ggplot2::ggplot(df_sni, ggplot2::aes(x = x, y = sni)) +
    ggplot2::geom_line(colour = "black", linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = 3, linetype = "dashed") +
    ggplot2::scale_x_reverse(limits = c(max(zone_div), min(zone_div)),
                              breaks = x_ticks,
                              labels = x_ticks / 1000,
                              expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(limits = y_lim_c,
                                breaks = y_tick_c,
                                expand = ggplot2::expansion(0)) +
    .char_theme() +
    ggplot2::labs(
      x     = "time (cal. yr BP x 1000)",
      y     = "signal-to-noise index",
      title = "(c) Local signal-to-noise index"
    )

  # ============================================================
  # PANEL (d): Global SNI boxplot / distribution
  # ============================================================
  df_sni_d <- data.frame(x = "SNI", y = sni_series[!is.na(sni_series)])

  pd <- ggplot2::ggplot(df_sni_d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_boxplot(outlier.shape = 16L, outlier.size = 1.2,
                          width = 0.5) +
    ggplot2::geom_hline(yintercept = 3, linetype = "dashed") +
    ggplot2::scale_y_continuous(limits = y_lim_c,
                                breaks = y_tick_c,
                                expand = ggplot2::expansion(0)) +
    ggplot2::annotate("text",
                      x = 1.35, y = stats::median(df_sni_d$y, na.rm = TRUE),
                      label = round(stats::median(df_sni_d$y, na.rm = TRUE), 2),
                      hjust = 0, size = 3, fontface = "plain") +
    .char_theme() +
    ggplot2::theme(axis.title.y  = ggplot2::element_blank(),
                   axis.text.y   = ggplot2::element_blank(),
                   axis.ticks.y  = ggplot2::element_blank()) +
    ggplot2::labs(
      x     = "global signal-to-\nnoise distribution",
      title = "(d) Global signal-to-noise index"
    )

  # ============================================================
  # Compose with patchwork (mirrors MATLAB 3x5 layout)
  # Wide panels (a, c) take 4 units; narrow panels (b, d) take 1 unit.
  # ============================================================
  design <- "AAAAB\nCCCCD"
  fig <- patchwork::wrap_plots(pa, pb, pc, pd) +
    patchwork::plot_layout(design = design)
  print(fig)
  invisible(fig)
}


# =============================================================================
# char_plot_all  --  produce and optionally save all figures
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_all <- function(out, save = FALSE, out_dir = ".", width = 11, height = 8.5) {

  .require_ggplot2()

  all_figs <- isTRUE(out$results$allFigures == 1L)

  fig1 <- if (all_figs) char_plot_raw(out) else NULL
  fig2 <- if (all_figs) char_plot_thresh_diag(out) else NULL
  fig3 <- char_plot_peaks(out)
  fig4 <- char_plot_sni(out)
  fig5 <- char_plot_cumulative(out)
  fig6 <- char_plot_fri(out)
  fig7 <- char_plot_fire_history(out)
  fig8 <- char_plot_zones(out)

  for (fig in list(fig3, fig4, fig5, fig6, fig7, fig8)) {
    if (!is.null(fig)) print(fig)
  }

  if (save) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    site <- out$site
    figs <- list(
      list(fig = fig1, name = "01_pretreatment"),
      list(fig = fig2, name = "02_threshold_determination"),
      list(fig = fig3, name = "03_CHAR_analysis"),
      list(fig = fig4, name = "04_threshold_sensitivity_SNI"),
      list(fig = fig5, name = "05_cumulative_peaks"),
      list(fig = fig6, name = "06_FRI_distributions"),
      list(fig = fig7, name = "07_continuous_fire_hx"),
      list(fig = fig8, name = "08_zone_comparisons")
    )
    for (f in figs) {
      if (!is.null(f$fig)) {
        path <- file.path(out_dir, paste0(site, "_", f$name, ".pdf"))
        ggplot2::ggsave(path, plot = f$fig, width = width, height = height,
                        units = "in", device = "pdf")
        message("Saved: ", path)
      }
    }
  }

  invisible(list(
    fig_pretreatment       = fig1,
    fig_threshold          = fig2,
    fig_char               = fig3,
    fig_threshold_sni      = fig4,
    fig_cumulative         = fig5,
    fig_fri_dist           = fig6,
    fig_fire_history       = fig7,
    fig_zones              = fig8
  ))
}

# =============================================================================
# Null-coalescing operator (avoid importing rlang)
# =============================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b
