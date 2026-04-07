#' Calculate a global threshold value for charcoal peak identification
#'
#' Determines a single threshold applied uniformly across the entire record.
#' The noise component of C_peak is modelled as a zero-mean Gaussian
#' (residuals), a one-mean Gaussian (ratios), or via a Gaussian mixture
#' model (GMM).  Mirrors \code{CharThreshGlobal.m} from the MATLAB v2.0
#' codebase.
#'
#' @param charcoal     Named list containing \code{peak} (C_peak series) and
#'   \code{ybpI}.
#' @param pretreatment Named list with \code{zoneDiv}.
#' @param peak_analysis Named list with \code{cPeak}, \code{threshMethod},
#'   \code{threshValues}.
#' @param site         Character string (site name; unused in R, kept for API
#'   symmetry).
#' @param results      Named list (unused in R, kept for API symmetry).
#' @param plot_data    0/1 flag; ignored in R.
#' @param bkg_sens_in  0/1 flag; ignored in R (no sensitivity loop).
#'
#' @return Named list \code{char_thresh} with elements:
#'   \describe{
#'     \item{possible}{Numeric vector of 251 candidate threshold bins.}
#'     \item{pos}{Numeric matrix [N x 4]: positive threshold for each of the
#'       four \code{threshValues}.}
#'     \item{neg}{Numeric matrix [N x 4] (method 1) or [N x 1] (methods 2–3):
#'       negative threshold.}
#'     \item{noise_pdf}{Estimated noise PDF evaluated at \code{possible}
#'       (methods 2–3), or scalar \code{-99} (method 1).}
#'     \item{mu_hat}{Fitted noise-component mean.}
#'     \item{sigma_hat}{Fitted noise-component standard deviation.}
#'     \item{SNI}{Signal-to-noise index (scalar).}
#'     \item{GOF}{Sentinel vector (\code{-999}, length N).}
#'   }
#'
#' @details
#'   ## Gaussian assumption (threshMethod = 2)
#'   For residuals (\code{cPeak = 1}): noise is zero-mean; sigma is estimated
#'   from the negative half of C_peak, mirrored and pooled.
#'   For ratios (\code{cPeak = 2}): noise is one-mean; values are shifted to
#'   zero, mirrored, shifted back, and sigma estimated from the pooled set.
#'
#'   ## GMM (threshMethod = 3)
#'   MATLAB's \code{GaussianMixture(X, 2, 2, false)} is replicated by
#'   \code{gaussian_mixture_em(X)} from \code{utils_gaussian_mixture.R}, using
#'   the same first/last-point initialisation and loose convergence criterion as
#'   the original Bowman CLUSTER EM.  The noise component is identified as the
#'   Gaussian with the smaller mean (matching MATLAB's
#'   \code{noiseIdx = find(mu == min(mu), 1)}).
#'
#'   ## Bin-lookup for threshold values
#'   Percentile thresholds are mapped to the nearest bin in
#'   \code{possible} (251 equally-spaced values spanning C_peak range).
#'   The v2.0 bug fix is preserved: both sides of the \code{abs()} comparison
#'   use the CHAR-unit threshold value \code{thresh[i]}, not the raw percentile.
#'
#' @seealso [char_thresh_local()], [char_smooth()], [CharAnalysis()]

# Requires: gaussian_mixture_em() from utils_gaussian_mixture.R
# (source that file before sourcing this one, or ensure it is on the search path)

char_thresh_global <- function(charcoal, pretreatment, peak_analysis,
                                site = NULL, results = NULL,
                                plot_data = 0L, bkg_sens_in = 0L) {

  # Candidate threshold bins (251 values spanning the C_peak range)
  pos_thresh_bins <- seq(min(charcoal$peak, na.rm = TRUE),
                         max(charcoal$peak, na.rm = TRUE),
                         length.out = 251L)

  char_thresh <- list(possible = pos_thresh_bins)
  N           <- length(charcoal$peak)
  n_tv        <- length(peak_analysis$threshValues)

  mu_hat    <- NA_real_
  sigma_hat <- NA_real_

  # =========================================================================
  # METHOD 1: User-defined threshold
  # =========================================================================
  if (peak_analysis$threshMethod == 1L) {

    char_thresh$pos <- matrix(NA_real_, nrow = N, ncol = n_tv)

    for (i in seq_len(n_tv)) {
      tv    <- peak_analysis$threshValues[i]
      in1   <- which(pos_thresh_bins >= tv)[1L]
      in2   <- rev(which(pos_thresh_bins <= tv))[1L]
      in1   <- if (is.na(in1)) length(pos_thresh_bins) else in1
      in2   <- if (is.na(in2)) 1L                       else in2
      in_final <- if (abs(pos_thresh_bins[in1] - tv) <=
                      abs(pos_thresh_bins[in2] - tv)) in1 else in2
      char_thresh$pos[, i] <- pos_thresh_bins[in_final]
    }

    char_thresh$neg       <- matrix(-99, nrow = N, ncol = n_tv)
    char_thresh$noise_pdf <- -99

  }

  # =========================================================================
  # METHODS 2 AND 3: Data-defined threshold
  # =========================================================================
  if (peak_analysis$threshMethod > 1L) {

    # -- Estimate noise distribution ----------------------------------------
    if (peak_analysis$threshMethod == 2L) {

      if (peak_analysis$cPeak == 1L) {
        # Residuals: zero-mean Gaussian
        neg_vals  <- charcoal$peak[charcoal$peak <= 0 & !is.na(charcoal$peak)]
        sigma_hat <- stats::sd(c(neg_vals, abs(neg_vals)))
        mu_hat    <- 0
      } else {
        # Ratios: one-mean Gaussian
        sub_vals  <- charcoal$peak[charcoal$peak <= 1 &
                                     !is.na(charcoal$peak)] - 1
        sigma_hat <- stats::sd(c(sub_vals, abs(sub_vals)) + 1)
        mu_hat    <- 1
      }

      char_thresh$noise_pdf <- stats::dnorm(pos_thresh_bins, mu_hat, sigma_hat)
    }

    if (peak_analysis$threshMethod == 3L) {

      # GMM: replicate MATLAB's GaussianMixture(X, 2, 2, false) using the
      # direct R port gaussian_mixture_em() from utils_gaussian_mixture.R.
      gmm_data <- charcoal$peak[!is.na(charcoal$peak)]
      fit      <- tryCatch(
        gaussian_mixture_em(gmm_data, k = 2L),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        # EM failed: fall back to simple Gaussian (method 2 behaviour)
        if (peak_analysis$cPeak == 1L) {
          neg_vals  <- gmm_data[gmm_data <= 0]
          sigma_hat <- stats::sd(c(neg_vals, abs(neg_vals)))
          mu_hat    <- 0
        } else {
          sub_vals  <- gmm_data[gmm_data <= 1] - 1
          sigma_hat <- stats::sd(c(sub_vals, abs(sub_vals)) + 1)
          mu_hat    <- 1
        }
      } else {
        # gaussian_mixture_em returns mu and sigma sorted ascending
        mu_both    <- fit$mu
        sigma_both <- fit$sigma

        if (mu_both[1L] == mu_both[2L]) {
          # Poor GMM fit — re-fit with K = 3, take two smallest-mean components
          warning("char_thresh_global: poor GMM fit (mu1 == mu2). ",
                  "Re-fitting with K = 3.")
          fit3 <- tryCatch(
            gaussian_mixture_em(gmm_data, k = 3L),
            error = function(e) NULL
          )
          if (!is.null(fit3)) {
            mu_both    <- fit3$mu[1:2]
            sigma_both <- fit3$sigma[1:2]
          }
        }

        # Noise component = Gaussian with the smaller mean (index 1, sorted)
        mu_hat    <- mu_both[1L]
        sigma_hat <- sigma_both[1L]
      }

      char_thresh$noise_pdf <- stats::dnorm(pos_thresh_bins, mu_hat, sigma_hat)
    }

    # -- Map percentile thresholds to nearest bin ----------------------------
    thresh_vals <- stats::qnorm(peak_analysis$threshValues, mu_hat, sigma_hat)

    char_thresh$pos <- matrix(NA_real_, nrow = N, ncol = n_tv)

    for (i in seq_len(n_tv)) {
      tv    <- thresh_vals[i]
      in1   <- which(pos_thresh_bins >= tv)[1L]
      in2   <- rev(which(pos_thresh_bins <= tv))[1L]
      in1   <- if (is.na(in1)) length(pos_thresh_bins) else in1
      in2   <- if (is.na(in2)) 1L                       else in2
      # v2.0 bug fix: both operands use thresh_vals[i], not threshValues[i]
      in_final <- if (abs(pos_thresh_bins[in1] - tv) <=
                      abs(pos_thresh_bins[in2] - tv)) in1 else in2
      char_thresh$pos[, i] <- pos_thresh_bins[in_final]
    }

    # Negative threshold: mirror of the THIRD percentile value
    thresh_neg       <- stats::qnorm(1 - peak_analysis$threshValues[3L],
                                      mu_hat, sigma_hat)
    char_thresh$neg  <- matrix(thresh_neg, nrow = N, ncol = 1L)
  }

  # =========================================================================
  # SIGNAL-TO-NOISE INDEX (Kelly et al.)
  # =========================================================================
  t_pos  <- char_thresh$pos[1L, n_tv]    # threshold using final threshValue
  signal <- charcoal$peak[!is.na(charcoal$peak) & charcoal$peak >  t_pos]
  noise  <- charcoal$peak[!is.na(charcoal$peak) & charcoal$peak <= t_pos]

  char_thresh$SNI <- if (length(signal) > 0 && length(noise) > 2) {
    (1 / length(signal)) *
      sum((signal - mean(noise)) / stats::sd(noise)) *
      ((length(noise) - 2) / length(noise))
  } else {
    0
  }

  # =========================================================================
  # GOODNESS-OF-FIT PLACEHOLDER
  # For the global threshold a single KS test is not applied; fill with the
  # sentinel value -999 for downstream compatibility.
  # =========================================================================
  char_thresh$GOF <- rep(-999, N)

  char_thresh$mu_hat    <- mu_hat
  char_thresh$sigma_hat <- sigma_hat

  char_thresh
}
