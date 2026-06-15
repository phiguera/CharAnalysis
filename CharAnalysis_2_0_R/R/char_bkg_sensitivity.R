# Bare column names used inside aes() below, declared to satisfy R CMD check
# (same convention as char_plot_results.R).
utils::globalVariables(c(
  "threshold", "window", "peaks", "value",
  "fri_mean", "fri_lo", "fri_hi", "peaks_scaled"
))

#' Background-window sensitivity analysis
#'
#' Re-runs the core CharAnalysis peak-detection pipeline across a range of
#' background (low-frequency) smoothing windows and summarises how the results
#' change with that choice.  Mirrors \code{bkgCharSensitivity.m} from the MATLAB
#' v2.0 codebase (Figure 5).
#'
#' For each candidate window width the function repeats the four core stages:
#' [char_smooth()] -> compute C_peak -> [char_thresh_global()] or
#' [char_thresh_local()] -> [char_peak_id()], then stores summary statistics.
#'
#' The analysis branches on \code{peak_analysis$threshType}:
#' \describe{
#'   \item{Global threshold (\code{threshType = 1})}{Produces a filled-contour
#'     surface of the number of peaks identified as a function of the candidate
#'     threshold (C_peak units, x) and the smoothing-window width (yr, y).  The
#'     threshold grid is held fixed across windows (see Details) so the surface
#'     is rectangular.}
#'   \item{Local threshold (\code{threshType = 2})}{Produces a four-panel
#'     diagnostic: (a) noise goodness-of-fit (KS p-value) by window, (b)
#'     signal-to-noise index (SNI) by window with a reference line at 3, (c) the
#'     sum of median SNI and median GOF versus window, and (d) mean
#'     fire-return interval +/- SD with the peak count on a second axis.}
#' }
#'
#' @param out     A \code{CharAnalysis} object returned by [CharAnalysis()].
#'   Supplies the resampled charcoal series and all parameter lists.
#' @param save    Logical.  If \code{TRUE}, the figure is saved as a PDF in
#'   \code{out_dir}.  Default \code{FALSE}.
#' @param out_dir Directory for the saved PDF.  Required when \code{save = TRUE}.
#' @param width,height PDF dimensions in inches (default 11 x 8.5, landscape),
#'   matching the other \code{char_plot_*} figures.
#' @param verbose Logical.  Print per-iteration progress messages (default
#'   \code{TRUE}), mirroring the MATLAB console output.
#'
#' @return Invisibly, a named list:
#'   \describe{
#'     \item{bkSmooth}{Numeric vector of evaluated window widths (yr).}
#'     \item{z}{Peak-count matrix.  Global: \eqn{[nWin \times nBin]} over the
#'       fixed positive threshold grid.  Local: \eqn{[nWin \times nThresh]},
#'       one column per \code{threshValues} entry.}
#'     \item{SNI_i}{Global: vector \eqn{[nWin]} of scalar SNI.  Local: matrix
#'       \eqn{[nWin \times N]} of per-sample SNI.}
#'     \item{GOF_i}{Global: vector \eqn{[nWin]} (sentinel \code{-999}; GOF is
#'       not defined on the global path).  Local: matrix \eqn{[nWin \times N]}
#'       of per-sample KS p-values.}
#'     \item{mFRI}{Local only: matrix \eqn{[nWin \times 2]} of the mean and SD
#'       of fire-return intervals at the final \code{threshValues} entry.
#'       \code{NULL} on the global path.}
#'     \item{thresh_bins}{Global only: the fixed candidate threshold grid used
#'       across all windows.  \code{NULL} on the local path.}
#'     \item{fig}{The assembled \pkg{ggplot2} / \pkg{patchwork} figure (or
#'       \code{NULL} if \pkg{ggplot2} is not installed).}
#'   }
#'
#' @details
#'   ## Window range (mirrors MATLAB)
#'   The maximum window is \code{min(1000, floor(max(age)/100) * 100)} so it
#'   stays shorter than the record.  The minimum window is 100 yr on the global
#'   path, or \code{30 * yrInterp} floored to the nearest 100 yr (clamped to
#'   100-500) on the local path, which keeps at least ~30 samples in each local
#'   threshold window.  Windows are then evaluated in 100-yr steps.  (MATLAB
#'   errors when \code{30 * yrInterp} falls outside 100-600; here it is clamped
#'   and a short record that yields \code{bkgMin > bkgMax} raises a clear error.)
#'
#'   ## Fixed threshold grid (global path)
#'   In MATLAB the global-threshold contour relies on \code{bkgSensIn = 1}, which
#'   makes \code{CharThreshGlobal} reuse Figure 2's axis limits so every window
#'   shares one candidate-threshold grid.  Here that coupling is made explicit:
#'   the grid is computed once from the baseline run's C_peak range
#'   (\code{out$charcoal$peak}) and passed to [char_thresh_global()] via its
#'   \code{thresh_bins} argument on every iteration.  This keeps the peak-count
#'   matrix rectangular and reproduces the MATLAB intent without relying on
#'   plot state.
#'
#' @seealso [CharAnalysis()], [char_smooth()], [char_thresh_global()],
#'   [char_thresh_local()], [char_peak_id()]
#'
#' @examples
#' \donttest{
#'   params_file <- system.file("validation", "CO_charParams.csv",
#'                              package = "CharAnalysis")
#'   out <- CharAnalysis(params_file)
#'   sens <- char_bkg_sensitivity(out)
#'   print(sens$fig)
#' }
#' @export
char_bkg_sensitivity <- function(out, save = FALSE, out_dir = NULL,
                                  width = 11, height = 8.5, verbose = TRUE) {

  if (isTRUE(save)) {
    if (is.null(out_dir))
      stop("When save = TRUE, please supply 'out_dir' for the saved PDF. ",
           "Use tempdir() for a transient location, or a path of your ",
           "choosing.", call. = FALSE)
    if (!is.character(out_dir) || length(out_dir) != 1L || nchar(out_dir) == 0L)
      stop("'out_dir' must be a non-empty character string.", call. = FALSE)
  }

  charcoal      <- out$charcoal
  pretreatment  <- out$pretreatment
  smoothing     <- out$smoothing
  peak_analysis <- out$peak_analysis
  results       <- out$results
  site          <- out$site

  thresh_type <- peak_analysis$threshType
  c_peak_mode <- peak_analysis$cPeak
  yr_interp   <- pretreatment$yrInterp

  # ---------------------------------------------------------------------------
  # Window range (see Details).
  # ---------------------------------------------------------------------------
  max_age <- max(charcoal$ybpI, na.rm = TRUE)
  bkg_max <- if (max_age >= 1000) 1000 else floor(max_age / 100) * 100
  if (bkg_max < 100) bkg_max <- 100

  if (thresh_type == 1L) {
    bkg_min <- 100
  } else {
    bkg_min <- floor((30 * yr_interp) / 100) * 100
    bkg_min <- max(100, min(500, bkg_min))
  }

  if (bkg_min > bkg_max) {
    stop("char_bkg_sensitivity: record is too short for a background ",
         "sensitivity analysis (minimum window ", bkg_min, " yr exceeds ",
         "maximum window ", bkg_max, " yr). At least two 100-yr windows ",
         "below the record length are required.", call. = FALSE)
  }

  bk_smooth <- seq(bkg_min, bkg_max, by = 100)
  n_win     <- length(bk_smooth)
  if (n_win < 2L) {
    stop("char_bkg_sensitivity: only one candidate window (", bkg_min,
         " yr) fits this record; a sensitivity analysis needs at least two.",
         call. = FALSE)
  }

  N <- length(charcoal$peak)

  # ---------------------------------------------------------------------------
  # Allocate outputs and (global path) the fixed threshold grid.
  # ---------------------------------------------------------------------------
  if (thresh_type == 1L) {
    thresh_bins <- seq(min(charcoal$peak, na.rm = TRUE),
                       max(charcoal$peak, na.rm = TRUE),
                       length.out = 251L)
    n_pos_bins  <- sum(thresh_bins > 0)
    z      <- matrix(NA_real_, nrow = n_win, ncol = n_pos_bins)
    SNI_i  <- rep(NA_real_, n_win)
    GOF_i  <- rep(NA_real_, n_win)
    mFRI   <- NULL
  } else {
    thresh_bins <- NULL
    n_tv   <- length(peak_analysis$threshValues)
    z      <- matrix(NA_real_, nrow = n_win, ncol = n_tv)
    SNI_i  <- matrix(NA_real_, nrow = n_win, ncol = N)
    GOF_i  <- matrix(NA_real_, nrow = n_win, ncol = N)
    mFRI   <- matrix(NA_real_, nrow = n_win, ncol = 2)  # [mean, sd]
  }

  # ---------------------------------------------------------------------------
  # Sensitivity loop.
  # ---------------------------------------------------------------------------
  for (i in seq_len(n_win)) {

    if (isTRUE(verbose)) {
      message(sprintf(
        "    C_background sensitivity iteration %d of %d: window = %d yr.",
        i, n_win, bk_smooth[i]))
    }

    smoothing_i    <- smoothing
    smoothing_i$yr <- bk_smooth[i]

    # Fresh copy of the resampled inputs for this iteration.
    cc <- charcoal

    # (a) Smooth to estimate C_background for this window.
    cc <- char_smooth(cc, pretreatment, smoothing_i, results, plot_data = 0L)

    # (b) Compute C_peak.
    if (c_peak_mode == 1L) {
      cc$peak <- cc$accI - cc$accIS                 # residuals
    } else {
      if (any(cc$accIS == 0, na.rm = TRUE)) {
        warning("char_bkg_sensitivity: C_background = 0 at window ",
                bk_smooth[i], " yr; ratio C_peak undefined, skipping window.",
                call. = FALSE)
        next
      }
      cc$peak <- cc$accI / cc$accIS                 # ratios
    }

    # (c) Define thresholds.
    if (thresh_type == 1L) {
      ct <- char_thresh_global(cc, pretreatment, peak_analysis,
                               site, results, plot_data = 0L,
                               bkg_sens_in = 1L, thresh_bins = thresh_bins)
    } else {
      ct <- char_thresh_local(cc, smoothing_i, peak_analysis,
                              site, results, plot_data = 0L)
    }

    # (d) Identify peaks.
    pid <- char_peak_id(cc, pretreatment, peak_analysis, ct)
    cc  <- pid$charcoal
    ct  <- pid$char_thresh

    # (e) Store results.
    if (thresh_type == 1L) {
      z[i, ]   <- colSums(cc$charPeaks, na.rm = TRUE)
      SNI_i[i] <- ct$SNI[1L]
      GOF_i[i] <- ct$GOF[1L]
    } else {
      z[i, ]      <- colSums(cc$charPeaks, na.rm = TRUE)
      SNI_i[i, ]  <- ct$SNI
      GOF_i[i, ]  <- ct$GOF
      # Mean / SD of FRIs at the final threshValues column.
      peak_ages   <- cc$ybpI[cc$charPeaks[, ncol(cc$charPeaks)] > 0]
      if (length(peak_ages) >= 2L) {
        fri        <- diff(sort(peak_ages))
        mFRI[i, 1] <- mean(fri, na.rm = TRUE)
        mFRI[i, 2] <- stats::sd(fri, na.rm = TRUE)
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Build the figure.  The numeric results above do not require ggplot2; only
  # the figure does.  If ggplot2 is unavailable, return the numbers with
  # fig = NULL rather than failing (bkgSens is a param-driven opt-in, so the
  # numeric sensitivity is the part the user asked for).
  # ---------------------------------------------------------------------------
  have_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  if (have_ggplot) {
    fig <- if (thresh_type == 1L) {
      .bkg_sens_plot_global(z, bk_smooth, thresh_bins, site)
    } else {
      .bkg_sens_plot_local(z, SNI_i, GOF_i, mFRI, bk_smooth, site)
    }
  } else {
    message("Install 'ggplot2' to render the sensitivity figure: ",
            "install.packages('ggplot2'). Returning numeric results only.")
    fig <- NULL
  }

  # ---------------------------------------------------------------------------
  # Optionally save (PDF only, matching the other char_plot_* figures).
  # ---------------------------------------------------------------------------
  if (isTRUE(save)) {
    if (is.null(fig)) {
      warning("char_bkg_sensitivity: cannot save figure because 'ggplot2' is ",
              "not installed.", call. = FALSE)
    } else {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      path <- file.path(out_dir,
                        paste0(site, "_05_sensitivity_to_C_background.pdf"))
      ggplot2::ggsave(path, plot = fig, width = width, height = height,
                      units = "in", device = "pdf")
      message("Saved: ", path)
    }
  }

  invisible(list(
    bkSmooth    = bk_smooth,
    z           = z,
    SNI_i       = SNI_i,
    GOF_i       = GOF_i,
    mFRI        = mFRI,
    thresh_bins = thresh_bins,
    fig         = fig
  ))
}

# =============================================================================
# Internal plot helpers
# =============================================================================

# Relabel geom_contour_filled() interval labels (e.g. "(0, 10]", "(10, 20]")
# as clean, non-overlapping integer ranges ("0-10", "11-20", ...).  Counts are
# whole numbers, so the lower edge of each bin (exclusive) is shown as lo + 1,
# except the first bin which starts at 0.
.bkg_bin_labels <- function(lvls) {
  nums <- regmatches(lvls, gregexpr("-?[0-9.]+", lvls))
  vapply(nums, function(p) {
    if (length(p) < 2L) return(paste(p, collapse = "-"))
    lo <- as.numeric(p[1L]); hi <- as.numeric(p[2L])
    lo_disp <- if (lo <= 0) 0 else floor(lo) + 1
    paste0(lo_disp, "-", hi)
  }, character(1L))
}

# Global threshold: filled-contour surface of peaks vs (threshold, window).
# Mirrors MATLAB Figure 10 (R: Figure 5) contourf() with the inverted-gray colormap and the
# x-axis truncated to the lower half of the threshold range.
.bkg_sens_plot_global <- function(z, bk_smooth, thresh_bins, site) {

  x <- thresh_bins[thresh_bins > 0]          # positive candidate thresholds

  df <- expand.grid(threshold = x, window = bk_smooth)
  df$peaks <- as.vector(t(z))                # row i = window i, across x

  lab_x <- if (.has_ggtext()) "threshold (C<sub>peak</sub> units)" else
                              "threshold (C_peak units)"
  ttl   <- .title(paste0(
    site, ": peaks identified across threshold and ",
    "C<sub>background</sub> window"))

  ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = window,
                                   z = peaks)) +
    ggplot2::geom_contour_filled() +
    ggplot2::scale_fill_grey(start = 0.92, end = 0.15,
                             name   = "peaks\nidentified (#)",
                             labels = .bkg_bin_labels) +
    ggplot2::coord_cartesian(xlim = c(min(x), stats::median(x))) +
    ggplot2::scale_y_continuous(breaks = bk_smooth) +
    ggplot2::xlab(lab_x) +
    ggplot2::ylab("smoothing window width (yr)") +
    ggplot2::labs(title = ttl) +
    .char_theme()
}

# Local threshold: 2 x 2 diagnostic panels.  Mirrors MATLAB Figure 10 (R: Figure 5)
# subplots (a)-(d).
.bkg_sens_plot_local <- function(z, SNI_i, GOF_i, mFRI, bk_smooth, site) {

  win_f <- factor(bk_smooth, levels = bk_smooth)

  # ---- long frames for the boxplots ----------------------------------------
  to_long <- function(mat) {
    df <- data.frame(
      window = rep(win_f, times = ncol(mat)),
      value  = as.vector(mat)
    )
    df[is.finite(df$value), , drop = FALSE]
  }
  gof_long <- to_long(GOF_i)
  sni_long <- to_long(SNI_i)

  # (a) Noise goodness of fit (KS p-value) ------------------------------------
  pa <- ggplot2::ggplot(gof_long,
                        ggplot2::aes(x = window, y = value)) +
    ggplot2::geom_boxplot(outlier.size = 0.6, fill = "grey85") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::xlab("smoothing window width (yr)") +
    ggplot2::ylab("KS-test result (p-value)") +
    ggplot2::labs(title = .title("(a) Noise distribution goodness of fit")) +
    .char_theme()

  # (b) Signal-to-noise index --------------------------------------------------
  sni_top <- min(max(sni_long$value, na.rm = TRUE), 10)
  pb <- ggplot2::ggplot(sni_long,
                        ggplot2::aes(x = window, y = value)) +
    ggplot2::geom_boxplot(outlier.size = 0.6, fill = "grey85") +
    ggplot2::geom_hline(yintercept = 3, linetype = "dashed",
                        colour = "black") +
    ggplot2::coord_cartesian(ylim = c(0, sni_top)) +
    ggplot2::xlab("smoothing window width (yr)") +
    ggplot2::ylab("signal-to-noise index") +
    ggplot2::labs(title = .title("(b) Signal-to-noise index")) +
    .char_theme()

  # (c) Sum of median SNI and median GOF --------------------------------------
  med_gof <- apply(GOF_i, 1L, stats::median, na.rm = TRUE)
  med_sni <- apply(SNI_i, 1L, stats::median, na.rm = TRUE)
  df_c <- data.frame(window = bk_smooth, value = med_sni + med_gof)
  x_lim <- range(bk_smooth) + c(-0.5, 0.5) * mean(diff(bk_smooth))
  pc <- ggplot2::ggplot(df_c,
                        ggplot2::aes(x = window, y = value)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::coord_cartesian(xlim = x_lim) +
    ggplot2::scale_x_continuous(breaks = bk_smooth) +
    ggplot2::xlab("smoothing window width (yr)") +
    ggplot2::ylab("sum of median SNI & GOF") +
    ggplot2::labs(title = .title("(c) Sum of median SNI and GOF")) +
    .char_theme()

  # (d) Mean FRI +/- SD (left) and peak count (right) -------------------------
  fri_mean <- mFRI[, 1]
  fri_sd   <- mFRI[, 2]
  peaks    <- z[, ncol(z)]

  fri_lo <- fri_mean - fri_sd
  fri_hi <- fri_mean + fri_sd
  rng_fri <- range(c(fri_lo, fri_hi, fri_mean), na.rm = TRUE)
  rng_pk  <- range(peaks, na.rm = TRUE)
  span_fri <- diff(rng_fri); if (span_fri == 0) span_fri <- 1
  span_pk  <- diff(rng_pk);  if (span_pk  == 0) span_pk  <- 1

  to_fri  <- function(p) (p - rng_pk[1]) / span_pk * span_fri + rng_fri[1]
  to_peak <- function(f) (f - rng_fri[1]) / span_fri * span_pk + rng_pk[1]

  df_d <- data.frame(window = bk_smooth, fri_mean = fri_mean,
                     fri_lo = fri_lo, fri_hi = fri_hi,
                     peaks_scaled = to_fri(peaks))

  pd <- ggplot2::ggplot(df_d, ggplot2::aes(x = window)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = fri_lo,
                                        ymax = fri_hi),
                           width = 0.4 * mean(diff(bk_smooth)),
                           colour = "black") +
    ggplot2::geom_line(ggplot2::aes(y = fri_mean),
                       colour = "black", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = fri_mean),
                        colour = "black", size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = peaks_scaled),
                       colour = "blue", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = peaks_scaled),
                        colour = "blue", size = 2) +
    ggplot2::scale_y_continuous(
      name     = "fire-return interval (yr)",
      sec.axis = ggplot2::sec_axis(~ to_peak(.),
                                   name = "# of peaks identified")) +
    ggplot2::coord_cartesian(xlim = x_lim) +
    ggplot2::scale_x_continuous(breaks = bk_smooth) +
    ggplot2::xlab("smoothing window width (yr)") +
    ggplot2::labs(title = .title(
      "(d) Mean fire-return interval +/- SD")) +
    .char_theme() +
    ggplot2::theme(axis.title.y.right = ggplot2::element_text(colour = "blue"),
                   axis.text.y.right  = ggplot2::element_text(colour = "blue"))

  # ---- assemble 2 x 2 --------------------------------------------------------
  if (requireNamespace("patchwork", quietly = TRUE)) {
    fig <- (pa | pb) / (pc | pd)
    fig <- fig + patchwork::plot_annotation(
      title = paste0(site,
                     ": sensitivity to background smoothing window"))
    return(fig)
  }

  message("Install 'patchwork' for the combined 2x2 layout: ",
          "install.packages('patchwork'). Returning panel (d) only.")
  pd
}
