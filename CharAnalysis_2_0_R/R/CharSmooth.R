#' Smooth charcoal record to estimate low-frequency C_background
#'
#' Applies one of five smoothing methods to the interpolated charcoal
#' accumulation rate series (\code{charcoal$accI}) and stores the result in
#' \code{charcoal$accIS}.  Mirrors \code{CharSmooth.m} from the MATLAB v2.0
#' codebase.
#'
#' @param charcoal     Named list with elements \code{accI} (resampled CHAR,
#'   length N), \code{ybpI} (resampled ages), and \code{ybp} (raw ages).
#' @param pretreatment Named list with element \code{yrInterp}.
#' @param smoothing    Named list with elements:
#'   \describe{
#'     \item{method}{Integer 1–5 selecting the smoothing method.}
#'     \item{yr}{Window width in years.}
#'   }
#' @param results      Named list (not used in R; kept for API symmetry).
#' @param plot_data    0/1 flag; ignored in R (no diagnostic plots).
#'
#' @return The input \code{charcoal} list with an additional element
#'   \code{accIS}: the smoothed C_background series (length N).
#'
#' @details
#'   ## Smoothing methods
#'   | Index | Name | Implementation |
#'   |-------|------|----------------|
#'   | 1 | Lowess | [char_lowess()] with \code{iter = 0} |
#'   | 2 | Robust Lowess | [char_lowess()] with \code{iter = 4} |
#'   | 3 | Moving average | \code{zoo::rollapply(..., mean, partial=TRUE)} |
#'   | 4 | Running median + Lowess | Shifted-window median loop, then [char_lowess()] |
#'   | 5 | Running mode + Lowess | Shifted-window 100-bin mode loop, then [char_lowess()] |
#'
#'   ## Span convention
#'   The smoothing window width in data-point units is
#'   \code{s = smoothing$yr / pretreatment$yrInterp}.  This is passed to
#'   [char_lowess()] as \code{span = s} (number of points), which
#'   converts it to the fraction required by \code{stats::lowess()} via
#'   \code{f = round(s) / N}.
#'
#'   ## NaN bridging
#'   NaN values in \code{accI} (from record gaps) cannot be passed to
#'   [char_lowess()] directly.  They are bridged by linear interpolation
#'   before smoothing and restored afterward, matching the MATLAB fallback
#'   path in \code{CharSmooth.m} (used when the Curve Fitting Toolbox is
#'   absent).  Methods 4 and 5 always use the bridged series.
#'
#'   ## Window selection for methods 4 and 5
#'   The boundary window is SHIFTED (not shrunk) to maintain exactly
#'   \code{round(s)} samples, matching MATLAB's loop logic.  Note that
#'   MATLAB's \code{round()} rounds 0.5 away from zero while R uses banker's
#'   rounding (round half to even); this can cause single-sample differences
#'   in window boundaries for odd half-integers.
#'
#' @seealso [char_lowess()], [CharAnalysis()]
char_smooth <- function(charcoal, pretreatment, smoothing,
                         results = NULL, plot_data = 0L) {

  r <- pretreatment$yrInterp            # yr per sample
  s <- smoothing$yr / r                 # window width in data-point units
  N <- length(charcoal$accI)

  # Preallocate all five columns (N x 5)
  char_acc_IS <- matrix(NA_real_, nrow = N, ncol = 5L)

  # =========================================================================
  # NaN BRIDGING
  # NaN entries in accI arise from record gaps. charLowess / lowess() require
  # a NaN-free input. Bridge NaN values by linear interpolation over the gap,
  # then restore NaN positions after all smoothing is done.
  # =========================================================================
  nan_mask   <- is.na(charcoal$accI)
  acc_clean  <- charcoal$accI

  if (any(nan_mask)) {
    idx_all  <- seq_len(N)
    acc_clean[nan_mask] <- approx(
      x    = idx_all[!nan_mask],
      y    = charcoal$accI[!nan_mask],
      xout = idx_all[nan_mask],
      rule = 2L            # extrapolate at boundaries if gap extends to edge
    )$y
  }

  # =========================================================================
  # METHOD 1: Lowess
  # Locally-weighted linear regression, no robustness iterations.
  # Mirrors: smooth(Charcoal.accI, s, 'lowess')
  # =========================================================================
  char_acc_IS[, 1L] <- char_lowess(acc_clean, span = s, iter = 0L)

  # =========================================================================
  # METHOD 2: Robust Lowess
  # Same as method 1 with bisquare re-weighting.
  # MATLAB charLowess uses nIter = 5 (1 initial + 4 robustness updates);
  # this maps to iter = 4 in R's lowess().
  # Mirrors: smooth(Charcoal.accI, s, 'rlowess')
  # =========================================================================
  char_acc_IS[, 2L] <- char_lowess(acc_clean, span = s, iter = 4L)

  # =========================================================================
  # METHOD 3: Moving average
  # Simple arithmetic mean within window; window shrinks at boundaries.
  # Mirrors: smooth(Charcoal.accI, s, 'moving')
  #          → charLowess(accI_clean, s, 'moving')
  #          → movmean(y, k, 'EndPoints', 'shrink')
  # =========================================================================
  k3 <- max(3L, min(round(s), N))
  char_acc_IS[, 3L] <- as.numeric(
    zoo::rollapply(acc_clean, width = k3, FUN = mean,
                   align = "center", partial = TRUE)
  )

  # =========================================================================
  # METHOD 4: Running median + Lowess
  #
  # Step 1: assign each sample the median accI value in a shifted window of
  #         exactly round(s) samples. Window is SHIFTED at boundaries (not
  #         shrunk), matching MATLAB's loop logic in CharSmooth.m lines 101-116.
  # Step 2: apply char_lowess() to the median series.
  # =========================================================================
  k4  <- max(3L, min(round(s), N))
  hw4 <- round(s / 2)

  for (i in seq_len(N)) {
    if (i <= hw4) {
      win <- acc_clean[seq_len(k4)]
    } else if (i >= N - hw4) {
      win <- acc_clean[max(1L, N - k4 + 1L):N]
    } else {
      win <- acc_clean[round(i - 0.5 * s):round(i + 0.5 * s)]
    }
    char_acc_IS[i, 4L] <- stats::median(win)
  }

  char_acc_IS[, 4L] <- char_lowess(char_acc_IS[, 4L], span = s, iter = 0L)

  # =========================================================================
  # METHOD 5: Running mode + Lowess
  #
  # Step 1: divide accI values within the (shifted) window into 100
  #         equally-spaced bins; assign each sample the centre of the
  #         most-populated bin.  If multiple bins share the maximum count,
  #         their centres are averaged (median of modal centres).
  #         Mirrors MATLAB's charHistCounts(win, nBin) + median(modal_centres).
  # Step 2: apply char_lowess() to the mode series.
  # =========================================================================
  n_bin <- 100L
  k5    <- max(3L, min(round(s), N))
  hw5   <- round(s / 2)

  for (i in seq_len(N)) {
    if (i <= hw5) {
      win <- acc_clean[seq_len(k5)]
    } else if (i >= N - hw5) {
      win <- acc_clean[max(1L, N - k5 + 1L):N]
    } else {
      win <- acc_clean[round(i - 0.5 * s):round(i + 0.5 * s)]
    }

    # 100-bin histogram over the window range
    bin_breaks  <- seq(min(win), max(win), length.out = n_bin + 1L)
    bin_centers <- (bin_breaks[-1L] + bin_breaks[-(n_bin + 1L)]) / 2
    # Assign each value to a bin (rightmost bin is closed on both sides)
    bin_idx     <- findInterval(win, bin_breaks, rightmost.closed = TRUE)
    bin_idx     <- pmax(1L, pmin(bin_idx, n_bin))  # clamp to [1, n_bin]
    n_mode      <- tabulate(bin_idx, nbins = n_bin)
    modal_ctrs  <- bin_centers[n_mode == max(n_mode)]
    char_acc_IS[i, 5L] <- stats::median(modal_ctrs)
  }

  char_acc_IS[, 5L] <- char_lowess(char_acc_IS[, 5L], span = s, iter = 0L)

  # =========================================================================
  # RESTORE NaN POSITIONS
  # Set smoothed values back to NA wherever accI was NA.
  # =========================================================================
  if (any(nan_mask)) {
    char_acc_IS[nan_mask, ] <- NA_real_
  }

  # =========================================================================
  # STORE SELECTED METHOD
  # =========================================================================
  charcoal$accIS <- char_acc_IS[, smoothing$method]

  charcoal
}
