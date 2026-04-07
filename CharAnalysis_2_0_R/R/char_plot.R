#' CharAnalysis ggplot2 figures
#'
#' Two primary output figures mirroring the MATLAB CharAnalysis v2.0 plots:
#'
#' \describe{
#'   \item{[char_plot_char()]}{Figure 3: Resampled CHAR with background trend
#'     (top panel) and C_peak with thresholds and peak markers (bottom panel).}
#'   \item{[char_plot_fire_history()]}{Figure 7: Peak magnitude (top), FRIs
#'     through time with smoothed FRI and CI ribbon (middle), and smoothed
#'     fire frequency (bottom).}
#'   \item{[char_plot_all()]}{Convenience wrapper: produces both figures and
#'     optionally saves them as PDF files.}
#' }
#'
#' @name char_plot
#' @aliases char_plot_char char_plot_fire_history char_plot_all
#'
#' @param out Named list returned by [CharAnalysis()].  Must contain
#'   \code{charcoal}, \code{pretreatment}, \code{peak_analysis},
#'   \code{char_thresh}, \code{post}, and \code{site}.
#' @param save   Logical.  If \code{TRUE}, each figure is saved as a PDF in
#'   \code{out_dir}.  Default \code{FALSE}.
#' @param out_dir Directory for saved PDFs.  Default: current working directory.
#' @param width,height PDF dimensions in inches.  Defaults: 11 × 8.5.
#'
#' @return
#'   \code{char_plot_char()} and \code{char_plot_fire_history()} each return
#'   a \pkg{patchwork} / \pkg{ggplot2} object.
#'   \code{char_plot_all()} returns a named list with elements
#'   \code{fig_char} and \code{fig_fire_history}.
#'
#' @details
#'   Requires the \pkg{ggplot2} package.  Multi-panel layout uses
#'   \pkg{patchwork} if available; otherwise panels are printed separately
#'   with a message.
#'
#' @seealso [CharAnalysis()], [char_post_process()]
#'
#' @examples
#' \dontrun{
#'   out <- CharAnalysis("CO_charParams.csv")
#'   char_plot_char(out)
#'   char_plot_fire_history(out)
#'   # Save both to PDF:
#'   char_plot_all(out, save = TRUE, out_dir = "Results")
#' }
NULL

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
  ggplot2::theme_classic(base_size = base_size) +
  ggplot2::theme(
    axis.line         = ggplot2::element_line(colour = "black"),
    panel.grid        = ggplot2::element_blank(),
    plot.title        = title_el,
    axis.ticks        = ggplot2::element_line(colour = "black"),
    axis.ticks.length = ggplot2::unit(3, "pt")
  )
}

# Format a title string: keep HTML when ggtext is available, strip tags otherwise.
.title <- function(html) {
  if (.has_ggtext()) html else gsub("<[^>]+>", "", html)
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
    ggplot2::xlab("cal. yr BP (\u00d7 1000)")
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
  # Text labels centred in each zone
  n_zones <- length(zone_div) - 1L
  if (n_zones <= 1L) return(NULL)
  mids   <- (zone_div[-length(zone_div)] + zone_div[-1L]) / 2
  labels <- LETTERS[seq_len(n_zones)]
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
# char_plot_char  —  Figure 3: C_interp / C_background / C_peak
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_char <- function(out) {

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
        ": C<sub>interp</sub> (", yr_interp, " yr) and",
        " C<sub>bkg</sub> defined by ", smooth_yr_label, " trends"))
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
  cpeak_label <- .title(if (peak_type == "ratio") {
    "(b) C<sub>peak</sub> (C<sub>interp</sub> / C<sub>bkg</sub>), thresholds, and identified peaks"
  } else {
    "(b) C<sub>peak</sub> (C<sub>interp</sub> \u2212 C<sub>bkg</sub>), thresholds, and identified peaks"
  })

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
# char_plot_fire_history  —  Figure 7: peak magnitude / FRIs / fire frequency
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
        site, ": peak magnitude, FRIs through time, and fire frequency"))) +
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
      ggplot2::ylab(paste0("FRI (yr fire\u207b\u00b9)\n",
                            peak_analysis$peakFrequ, "-yr mean\n",
                            round((1 - alpha) * 100), "% CI")) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x  = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) +
      .char_theme()
  } else {
    p2 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                         label = "Fewer than 3 FRIs \u2014 FRI plot not available",
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
    ggplot2::ylab(paste0("fire freq.\n(fires ",
                          peak_analysis$peakFrequ, " yr\u207b\u00b9)")) +
    .char_theme()

  .combine_panels(list(p1, p2, p3), heights = c(1.2, 1.2, 1))
}

# =============================================================================
# char_plot_all  —  produce and optionally save both figures
# =============================================================================

#' @rdname char_plot
#' @export
char_plot_all <- function(out, save = FALSE, out_dir = ".", width = 11, height = 8.5) {

  .require_ggplot2()

  fig3 <- char_plot_char(out)
  fig7 <- char_plot_fire_history(out)

  print(fig3)
  print(fig7)

  if (save) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    site <- out$site
    f3 <- file.path(out_dir, paste0(site, "_03_CHAR_analysis.pdf"))
    f7 <- file.path(out_dir, paste0(site, "_07_continuous_fire_hx.pdf"))
    ggplot2::ggsave(f3, plot = fig3, width = width, height = height,
                    units = "in", device = "pdf")
    ggplot2::ggsave(f7, plot = fig7, width = width, height = height,
                    units = "in", device = "pdf")
    message("Saved: ", f3)
    message("Saved: ", f7)
  }

  invisible(list(fig_char = fig3, fig_fire_history = fig7))
}

# =============================================================================
# Null-coalescing operator (avoid importing rlang)
# =============================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b
