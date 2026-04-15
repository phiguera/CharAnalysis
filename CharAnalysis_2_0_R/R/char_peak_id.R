#' Identify charcoal peaks and apply minimum-count screening
#'
#' Mirrors \code{CharPeakID.m} from the MATLAB v2.0 codebase.  For each
#' threshold column, samples where C_peak exceeds the threshold are flagged,
#' consecutive exceedances are collapsed to a single event (keeping the
#' \strong{last} sample of each run, i.e. the oldest, matching the MATLAB
#' v1.1 / v2.0 convention), and each event is screened against a
#' minimum-count criterion using the Shuie-Bain (1982) extension of the
#' Detre-White (1970) test for unequal sediment volumes.
#'
#' @param charcoal     Named list containing:
#'   \describe{
#'     \item{peak}{C_peak series (\code{accI - accIS} or \code{accI / accIS}),
#'       length \eqn{N}.}
#'     \item{ybpI}{Age at each interpolated sample (cal yr BP), length \eqn{N}.}
#'     \item{accI}{Interpolated charcoal accumulation rate, length \eqn{N}.}
#'     \item{countI}{Interpolated charcoal count, length \eqn{N}.}
#'     \item{volI}{Interpolated sediment volume (cm^3), length \eqn{N}.}
#'   }
#' @param pretreatment Named list with element \code{yrInterp} (interpolation
#'   resolution, years).
#' @param peak_analysis Named list with elements \code{threshType} (1 = global,
#'   2 = local), \code{threshValues} (numeric vector, length \eqn{T}), and
#'   \code{minCountP} (alpha for minimum-count screen).
#' @param char_thresh  Named list returned by [char_thresh_global()] or
#'   [char_thresh_local()], containing \code{possible} (251-bin grid) and
#'   \code{pos} (\eqn{[N \times T]} threshold matrix).
#'
#' @return A named list with two components:
#'   \describe{
#'     \item{charcoal}{Input \code{charcoal} list augmented with:
#'       \itemize{
#'         \item \code{charPeaks}      -- \eqn{[N \times T]} numeric: 1 at peak
#'           samples, 0 elsewhere.
#'         \item \code{charPeaksThresh} -- \eqn{[N \times T]} numeric: threshold
#'           value at each identified peak, 0 elsewhere.
#'         \item \code{peaksTotal}     -- numeric vector length \eqn{T}: total
#'           peaks per threshold column.
#'         \item \code{threshFRI}      -- numeric matrix (\eqn{\leq N \times T}):
#'           fire-return intervals derived from peak ages per threshold column.
#'       }
#'     }
#'     \item{char_thresh}{Input \code{char_thresh} list augmented with
#'       \code{minCountP} -- \eqn{[N \times T]} matrix of Shuie-Bain p-values
#'       (NaN where not computed).}
#'   }
#'
#' @details
#'   ## Threshold value matrix
#'   For **global** thresholds (\code{threshType == 1}), \code{char_thresh$pos}
#'   is a constant-row \eqn{[N \times T]} matrix reused directly.  For
#'   **local** thresholds (\code{threshType == 2}), \code{char_thresh$pos} is
#'   already \eqn{[N \times T]} (per-sample values).
#'
#'   ## Consecutive-peak removal
#'   After flagging all exceedances, a diff-based pass retains only the
#'   \emph{last} sample of each consecutive run -- the oldest sample within a
#'   group of contiguous above-threshold values.  This matches the MATLAB v1.1
#'   algorithm (which the v2.0 comment documents correctly despite the v1.1
#'   comment being misleading).
#'
#'   ## Minimum-count test
#'   For each identified peak \eqn{i} in column \eqn{j}, a time window of
#'   \eqn{\pm 150} yr is constructed around the peak, then narrowed to the
#'   adjacent peaks when they fall within the window.  The test statistic is
#'   \deqn{
#'     d = \frac{|c_{\min} - (c_{\min}+c_{\max})\,v_{\min}/(v_{\min}+v_{\max})| - 0.5}
#'              {\sqrt{(c_{\min}+c_{\max})\,v_{\min}\,v_{\max}/(v_{\min}+v_{\max})^2}}
#'   }
#'   and the p-value is \eqn{1 - \Phi(d)} (standard normal CDF; equivalent to
#'   MATLAB's \code{1 - tcdf(d, 1e10)} because \eqn{t_{1\times10^{10}} \to z}).
#'   Peaks with \eqn{p > \alpha_{\text{peak}}} are removed.
#'
#' @seealso [char_thresh_local()], [char_thresh_global()], [CharAnalysis()]
char_peak_id <- function(charcoal, pretreatment, peak_analysis, char_thresh) {

  r          <- pretreatment$yrInterp
  N          <- length(charcoal$peak)
  alpha_peak <- peak_analysis$minCountP

  # ============================================================
  # BUILD THRESHOLD VALUE MATRIX  [N x T]
  # ============================================================
  # Global threshold: restrict possible grid to positive values, replicate.
  # Local threshold:  char_thresh$pos is already [N x T].

  if (peak_analysis$threshType == 1L) {
    pos_grid         <- char_thresh$possible[char_thresh$possible > 0]
    n_thresholds     <- length(pos_grid)
    threshold_values <- matrix(rep(pos_grid, each = N),
                               nrow = N, ncol = n_thresholds)
  } else {
    threshold_values <- char_thresh$pos
    n_thresholds     <- ncol(threshold_values)
  }

  # ============================================================
  # PEAK FLAGGING  (vectorised broadcast, matches MATLAB v2.0)
  # ============================================================
  # Exceedances are marked 2 so the consecutive-removal step can
  # distinguish "downgraded" (1) from "absent" (0).

  char_peaks <- 2 * (outer(charcoal$peak, rep(1, n_thresholds)) >
                       threshold_values) * 1.0

  # ============================================================
  # CONSECUTIVE-PEAK REMOVAL  (diff-based, column-wise)
  # ============================================================
  # Retains the LAST sample of each consecutive run of exceedances
  # (= the oldest sample within each event group), matching the
  # implemented behaviour of MATLAB CharPeakID.m v1.1 and v2.0.

  for (j in seq_len(n_thresholds)) {
    pk    <- char_peaks[, j] > 0          # logical length-N
    is_end <- c(diff(pk) < 0, pk[N])      # TRUE at run-end positions
    char_peaks[, j] <- as.numeric(is_end)
  }

  # ============================================================
  # charPeaksThresh: threshold value at each identified peak
  # ============================================================
  # Used downstream for plotting and peak magnitude.
  # MATLAB uses thresholdValues(1,:) -- the first row (same for global;
  # per-column constant for local since the threshold is per-sample but
  # the graph marker only needs a representative value).

  char_peaks_thresh <- sweep(char_peaks, 2, threshold_values[1L, ], `*`)

  # ============================================================
  # MINIMUM-COUNT ANALYSIS
  # ============================================================
  mc_window <- round(150 / r) * r        # search half-window in years

  d_mat           <- matrix(0,         nrow = N, ncol = n_thresholds)
  min_count_p_mat <- matrix(NA_real_,  nrow = N, ncol = n_thresholds)

  for (j in seq_len(n_thresholds)) {

    peak_index <- which(char_peaks[, j] > 0)

    if (length(peak_index) <= 1L) next    # need >= 2 peaks

    for (i in seq_along(peak_index)) {

      peak_yr_i <- charcoal$ybpI[peak_index[i]]

      # ---- Time-based window +/-mcWindow around peak ----------------
      # ybpI is ordered youngest->oldest (ascending yr BP values).
      # windowTime[1]: oldest boundary (max ybpI <= peakYr + mcWindow)
      # windowTime[2]: youngest boundary (min ybpI >= peakYr - mcWindow)
      wt1_candidates <- charcoal$ybpI[charcoal$ybpI <= peak_yr_i + mc_window]
      wt2_candidates <- charcoal$ybpI[charcoal$ybpI >= peak_yr_i - mc_window]
      if (length(wt1_candidates) == 0L || length(wt2_candidates) == 0L) next

      window_time_1 <- max(wt1_candidates)
      window_time_2 <- min(wt2_candidates)
      wt_in_1 <- which(charcoal$ybpI == window_time_1)[1L]   # older end index
      wt_in_2 <- which(charcoal$ybpI == window_time_2)[1L]   # younger end index

      # ---- Narrow window to adjacent peaks ------------------------
      peak_yrs_all <- charcoal$ybpI[peak_index]

      if (i == 1L) {
        # Youngest peak: older side bounded by next peak (i+1)
        wp_in_1 <- which(charcoal$ybpI == peak_yrs_all[i + 1L])[1L]
        wp_in_2 <- wt_in_2
      } else if (i == length(peak_index)) {
        # Oldest peak: younger side bounded by previous peak (i-1)
        wp_in_1 <- wt_in_1
        wp_in_2 <- which(charcoal$ybpI == peak_yrs_all[i - 1L])[1L]
      } else {
        wp_in_1 <- which(charcoal$ybpI == peak_yrs_all[i + 1L])[1L]
        wp_in_2 <- which(charcoal$ybpI == peak_yrs_all[i - 1L])[1L]
      }

      # Apply narrowing: use peak boundary when it falls within time window
      if (!is.na(wp_in_1) && wt_in_1 > wp_in_1) wt_in_1 <- wp_in_1
      if (!is.na(wp_in_2) && wt_in_2 < wp_in_2) wt_in_2 <- wp_in_2

      ws <- c(wt_in_1, wt_in_2)   # ws[1]=older index, ws[2]=younger index

      # ---- Count extremes -----------------------------------------
      # countMax: max count in [ws[2] : peakIndex] (younger side)
      seg_young <- charcoal$countI[ws[2L]:peak_index[i]]
      count_max    <- round(max(seg_young))
      count_max_in <- ws[2L] - 1L +
        which(round(charcoal$countI[ws[2L]:peak_index[i]]) == count_max)[1L]

      # countMin: min count in [peakIndex : ws[1]] (older side)
      seg_old  <- charcoal$countI[peak_index[i]:ws[1L]]
      count_min    <- round(min(seg_old))
      count_min_in <- peak_index[i] - 1L +
        which(round(charcoal$countI[peak_index[i]:ws[1L]]) == count_min)[1L]

      vol_max <- charcoal$volI[count_max_in]
      vol_min <- charcoal$volI[count_min_in]

      # ---- Shuie-Bain (1982) test statistic -----------------------
      # Detre-White / Shuie-Bain test for unequal sediment volumes
      # (Gavin 2006 / Charster).
      v_ratio <- vol_min / (vol_min + vol_max)
      d_val   <- (abs(count_min - (count_min + count_max) * v_ratio) - 0.5) /
                 sqrt((count_min + count_max) * v_ratio * (1 - v_ratio))

      d_mat[peak_index[i], j] <- d_val

      # p-value: t_{1e10} -> standard normal; pnorm(d) is equivalent.
      min_count_p_mat[peak_index[i], j] <- 1 - stats::pnorm(d_val)
    }
  }

  # ============================================================
  # REMOVE PEAKS FAILING MINIMUM-COUNT CRITERION
  # ============================================================

  for (j in seq_len(n_thresholds)) {
    fail_idx <- which(char_peaks[, j] > 0 &
                      !is.na(min_count_p_mat[, j]) &
                      min_count_p_mat[, j] > alpha_peak)
    char_peaks[fail_idx, j]        <- 0
    char_peaks_thresh[fail_idx, j] <- 0
  }

  # ============================================================
  # FIRE-RETURN INTERVALS AND PEAK TOTALS
  # ============================================================

  peaks_total <- numeric(n_thresholds)
  thresh_fri  <- matrix(NA_real_, nrow = N, ncol = n_thresholds)

  for (j in seq_len(n_thresholds)) {
    peaks_total[j] <- sum(char_peaks[, j])
    fri_yrs <- diff(charcoal$ybpI[char_peaks[, j] > 0])
    if (length(fri_yrs) > 0L) {
      thresh_fri[seq_along(fri_yrs), j] <- fri_yrs
    }
  }

  # ============================================================
  # PACK RESULTS
  # ============================================================

  charcoal$charPeaks       <- char_peaks
  charcoal$charPeaksThresh <- char_peaks_thresh
  charcoal$peaksTotal      <- peaks_total
  charcoal$threshFRI       <- thresh_fri

  char_thresh$minCountP    <- min_count_p_mat

  list(charcoal   = charcoal,
       char_thresh = char_thresh)
}
