#' Calculate a per-sample local threshold for charcoal peak identification
#'
#' Determines a sliding-window threshold value for each sample, based on the
#' distribution of C_peak within a window centred on that sample.  The noise
#' component within each window is modelled as a zero-mean Gaussian (residuals),
#' a one-mean Gaussian (ratios), or via a Gaussian mixture model (GMM).
#' Mirrors \code{CharThreshLocal.m} from the MATLAB v2.0 codebase.
#'
#' @param charcoal      Named list with \code{peak} (C_peak series), \code{ybpI}
#'   (resampled ages), and \code{accI}.
#' @param smoothing     Named list with \code{yr} (window width in years).
#' @param peak_analysis Named list with \code{cPeak}, \code{threshMethod},
#'   \code{threshValues} (length 4).
#' @param site          Character string (site name; unused in R, kept for API
#'   symmetry).
#' @param results       Named list (unused in R, kept for API symmetry).
#' @param plot_data     0/1 flag; ignored in R (no diagnostic plots).
#'
#' @return Named list \code{char_thresh} with elements:
#'   \describe{
#'     \item{pos}{Numeric matrix [N x 4]: per-sample positive threshold for each
#'       of the four \code{threshValues}, after Lowess smoothing.}
#'     \item{neg}{Numeric matrix [N x 4]: per-sample negative threshold, after
#'       Lowess smoothing.}
#'     \item{SNI}{Numeric vector [N]: signal-to-noise index time series, Lowess-
#'       smoothed and clamped to \eqn{\geq 0}.}
#'     \item{GOF}{Numeric vector [N]: KS goodness-of-fit p-values per sample
#'       (NA where fewer than 4 noise samples exist).}
#'   }
#'
#' @details
#'   ## Window selection
#'   The half-window is \code{hw = round(0.5 * smoothing$yr / r)} samples, where
#'   \code{r = mean(diff(charcoal$ybpI))} is the record resolution.  The window
#'   is shifted (not shrunk) at record boundaries, matching MATLAB's loop logic.
#'
#'   ## NaN bridging within windows
#'   NaN entries in \code{charcoal$peak} (from record gaps) are replaced with the
#'   neutral value (0 for residuals, 1 for ratios) before distribution fitting.
#'   This preserves window size while preventing gap samples from biasing the fit,
#'   matching the v2.0 bug fix documented in \code{CharThreshLocal.m}.
#'
#'   ## Small-sample fallback
#'   If fewer than 4 non-neutral samples exist in the window, a simple Gaussian
#'   is fitted (same as method 2) and the KS test is skipped.
#'
#'   ## GMM replacement (threshMethod = 3)
#'   MATLAB's \code{GaussianMixture(X, 2, 2, false)} is replaced by
#'   \code{mclust::Mclust(G = 2, modelNames = "V")}.  The noise component is
#'   identified as the Gaussian with the smaller mean.  Poor-fit fallback
#'   (mu1 == mu2) re-fits with G = 3 and takes the two components with the
#'   smallest means, mirroring the MATLAB v2.0 bug fix (the fallback now passes
#'   the local window X, not the full record).
#'
#'   ## KS goodness-of-fit
#'   MATLAB evaluates the fitted normal CDF at 101 equally-spaced bin centres
#'   and calls \code{kstest()} with a custom CDF table.  R uses
#'   \code{stats::ks.test(noise, "pnorm", mu, sigma)}, which evaluates the
#'   continuous CDF at each observation.  P-values may differ by a small amount
#'   for small samples; the statistic converges as sample size grows.
#'
#'   ## Post-loop smoothing
#'   After the per-sample loop, \code{pos}, \code{neg}, and \code{SNI} are
#'   smoothed with \code{char_lowess(span = smoothing$yr / r, iter = 0)}.
#'   SNI is then clamped to \eqn{\geq 0}.
#'
#' @seealso [char_thresh_global()], [char_smooth()], [CharAnalysis()]

# mclust must be attached (not just loaded) so its internal functions
# (e.g. mclustBIC) are on the search path when Mclust() runs.
if (!requireNamespace("mclust", quietly = TRUE))
  stop("Package 'mclust' is required: install.packages('mclust')")
library(mclust)

# Helper: extract per-component standard deviations from a mclust fit.
# Handles sigmasq vs sigma storage differences across mclust versions:
#   - mclust 5.x / 6.x univariate "V": parameters$variance$sigmasq (numeric vector)
#   - multivariate or alternate storage:  parameters$variance$sigma (3-D array)
#   - degenerate fallback: compute per-component SD from classification
.mclust_sigma <- function(fit) {
  var_s <- fit$parameters$variance
  if (!is.null(var_s$sigmasq) && is.numeric(var_s$sigmasq)) {
    return(sqrt(var_s$sigmasq))
  }
  if (!is.null(var_s$sigma) && is.array(var_s$sigma)) {
    return(sqrt(apply(var_s$sigma, 3L, function(m) m[1L, 1L])))
  }
  # Last resort: per-component SD from classification labels
  cl   <- fit$classification
  vapply(sort(unique(cl)), function(k) {
    vals <- fit$data[cl == k]
    if (length(vals) < 2L) 1e-8 else stats::sd(vals)
  }, numeric(1L))
}

char_thresh_local <- function(charcoal, smoothing, peak_analysis,
                               site = NULL, results = NULL,
                               plot_data = 0L) {

  r   <- mean(diff(charcoal$ybpI))          # yr per sample
  hw  <- round(0.5 * smoothing$yr / r)      # half-window in samples
  N   <- length(charcoal$peak)
  n_tv <- length(peak_analysis$threshValues)
  P   <- peak_analysis$threshValues[4L]     # percentile used for SNI

  # Neutral value replaces NaN gap entries before distribution fitting
  # (0 for residuals, 1 for ratios)
  neutral_val <- if (peak_analysis$cPeak == 1L) 0 else 1

  # Preallocate output arrays (NaN-initialised matching MATLAB)
  char_thresh        <- list()
  char_thresh$pos    <- matrix(NA_real_, nrow = N, ncol = n_tv)
  char_thresh$neg    <- matrix(NA_real_, nrow = N, ncol = n_tv)
  char_thresh$SNI    <- rep(NA_real_, N)
  char_thresh$GOF    <- rep(NA_real_, N)

  # ===========================================================================
  # PER-SAMPLE LOOP
  # ===========================================================================
  for (i in seq_len(N)) {

    # ---- Window selection (shifted at boundaries) ---------------------------
    # Mirrors MATLAB CharThreshLocal.m lines 86-92 (1-indexed)
    if (i <= hw) {
      X <- charcoal$peak[seq_len(hw + i)]
    } else if (i > N - hw) {
      X <- charcoal$peak[(i - hw):N]
    } else {
      X <- charcoal$peak[(i - hw):(i + hw)]
    }

    # Replace NaN (gap) values with neutral value to preserve window size.
    # v2.0 bug fix: stripping NaN changes window size and distribution shape
    # near record gaps; substitution of neutral value avoids this.
    X[is.na(X)] <- neutral_val

    # ---- Small-sample fallback: < 4 non-neutral samples --------------------
    # If the window is nearly all neutral (gap region), fit a simple Gaussian
    # instead of a GMM and skip the KS test.
    if (sum(X != neutral_val) < 4L) {

      if (peak_analysis$cPeak == 1L) {
        neg_vals    <- X[X <= 0]
        sigma_hat_i <- if (length(neg_vals) == 0L) {
          stats::sd(X)
        } else {
          stats::sd(c(neg_vals, abs(neg_vals)))
        }
        mu_hat_i <- 0
      } else {
        sub_vals    <- X[X <= 1] - 1
        sigma_hat_i <- stats::sd(c(sub_vals, abs(sub_vals)) + 1)
        mu_hat_i    <- 1
      }

      char_thresh$pos[i, ] <- stats::qnorm(peak_analysis$threshValues,
                                            mu_hat_i, sigma_hat_i)
      char_thresh$neg[i, ] <- stats::qnorm(1 - peak_analysis$threshValues,
                                            mu_hat_i, sigma_hat_i)

      t_pos   <- stats::qnorm(P, mu_hat_i, sigma_hat_i)
      sig_i   <- X[X >  t_pos]
      noise_i <- X[X <= t_pos]

      char_thresh$SNI[i] <- if (length(sig_i) > 0L && length(noise_i) > 2L) {
        (1 / length(sig_i)) *
          sum((sig_i - mean(noise_i)) / stats::sd(noise_i)) *
          ((length(noise_i) - 2L) / length(noise_i))
      } else {
        0
      }

      next  # KS test skipped for small-sample windows
    }

    # ---- Estimate local noise distribution ----------------------------------

    if (peak_analysis$threshMethod == 2L) {

      if (peak_analysis$cPeak == 1L) {
        # Residuals: zero-mean Gaussian
        neg_vals    <- X[X <= 0]
        sigma_hat_i <- stats::sd(c(neg_vals, abs(neg_vals)))
        mu_hat_i    <- 0
      } else {
        # Ratios: one-mean Gaussian
        sub_vals    <- X[X <= 1] - 1
        sigma_hat_i <- stats::sd(c(sub_vals, abs(sub_vals)) + 1)
        mu_hat_i    <- 1
      }

    } else if (peak_analysis$threshMethod == 3L) {

      if (sum(X) == 0) {
        # All-zero window: fall back to simple Gaussian
        # (mirrors MATLAB lines 157-173)
        if (peak_analysis$cPeak == 1L) {
          mu_hat_i    <- 0
          neg_vals    <- X[X <= 0]
          sigma_hat_i <- if (length(neg_vals) == 0L) stats::sd(X)
                         else stats::sd(c(neg_vals, abs(neg_vals)))
        } else {
          mu_hat_i    <- 1
          sub_vals    <- X[X <= 1] - 1
          sigma_hat_i <- stats::sd(c(sub_vals, abs(sub_vals)) + 1)
        }

      } else {
        # GMM: replace MATLAB GaussianMixture(X, 2, 2, false) with mclust.
        # "V" = variable variance (MATLAB EM does not constrain variances equal).
        # tryCatch guards against degenerate windows where mclust cannot converge.
        fit <- tryCatch(
          Mclust(data = X, G = 2L, modelNames = "V", verbose = FALSE),
          error = function(e) NULL
        )

        if (is.null(fit)) {
          # mclust failed to converge: fall back to simple Gaussian (method 2)
          if (peak_analysis$cPeak == 1L) {
            neg_vals    <- X[X <= 0]
            sigma_hat_i <- if (length(neg_vals) == 0L) stats::sd(X)
                           else stats::sd(c(neg_vals, abs(neg_vals)))
            mu_hat_i    <- 0
          } else {
            sub_vals    <- X[X <= 1] - 1
            sigma_hat_i <- stats::sd(c(sub_vals, abs(sub_vals)) + 1)
            mu_hat_i    <- 1
          }
        } else {
          mu_both    <- fit$parameters$mean
          sigma_both <- .mclust_sigma(fit)   # robust across mclust versions

          if (mu_both[1L] == mu_both[2L]) {
            # Poor GMM fit: re-fit with G = 3, take two smallest-mean components.
            # Mirrors the v2.0 bug fix (re-fits on local window X, not full record).
            warning("char_thresh_local: poor GMM fit at sample ", i,
                    " (mu1 == mu2). Re-fitting with G = 3.")
            fit3 <- tryCatch(
              Mclust(data = X, G = 3L, modelNames = "V", verbose = FALSE),
              error = function(e) NULL
            )
            if (!is.null(fit3)) {
              mu3        <- fit3$parameters$mean
              sigma3     <- .mclust_sigma(fit3)
              ord3       <- order(mu3)
              mu_both    <- mu3[ord3][1:2]
              sigma_both <- sigma3[ord3][1:2]
            }
            # If fit3 also failed, keep mu_both/sigma_both from original fit.
          }

          # Noise component = Gaussian with the smaller mean.
          # (matches MATLAB: noiseIdx = find(muHat == min(muHat), 1))
          noise_idx   <- which.min(mu_both)
          mu_hat_i    <- mu_both[noise_idx]
          sigma_hat_i <- sigma_both[noise_idx]
        }
      }
    }

    # ---- Compute local threshold values from inverse normal CDF ------------
    char_thresh$pos[i, ] <- stats::qnorm(peak_analysis$threshValues,
                                          mu_hat_i, sigma_hat_i)
    char_thresh$neg[i, ] <- stats::qnorm(1 - peak_analysis$threshValues,
                                          mu_hat_i, sigma_hat_i)

    # ---- Signal-to-noise index (Kelly et al.) --------------------------------
    t_pos   <- stats::qnorm(P, mu_hat_i, sigma_hat_i)
    sig_i   <- X[X >  t_pos]
    noise_i <- X[X <= t_pos]

    char_thresh$SNI[i] <- if (length(sig_i) > 0L && length(noise_i) > 2L) {
      (1 / length(sig_i)) *
        sum((sig_i - mean(noise_i)) / stats::sd(noise_i)) *
        ((length(noise_i) - 2L) / length(noise_i))
    } else {
      0
    }

    # ---- Goodness-of-fit: one-sample KS test --------------------------------
    # Noise subsample is all window values at or below the threshold.
    # MATLAB builds a 101-point CDF table then calls kstest(); R uses
    # ks.test() against the continuous normal CDF directly. P-values converge
    # for moderate sample sizes; small-sample differences are expected.
    noise_for_ks <- X[X <= t_pos]

    if (length(noise_for_ks) > 3L) {
      ks_result <- suppressWarnings(
        stats::ks.test(noise_for_ks, "pnorm", mu_hat_i, sigma_hat_i)
      )
      char_thresh$GOF[i] <- ks_result$p.value
    }

  }  # end per-sample loop

  # ===========================================================================
  # SMOOTH THRESHOLD AND SNI SERIES WITH LOWESS
  # Mirrors MATLAB CharThreshLocal.m lines 296-304.
  # span = threshYr / r  (same span used for background smoothing)
  # ===========================================================================
  span <- smoothing$yr / r

  char_thresh$SNI <- char_lowess(char_thresh$SNI, span = span, iter = 0L)
  char_thresh$SNI[char_thresh$SNI < 0] <- 0

  for (i in seq_len(n_tv)) {
    char_thresh$pos[, i] <- char_lowess(char_thresh$pos[, i],
                                         span = span, iter = 0L)
    char_thresh$neg[, i] <- char_lowess(char_thresh$neg[, i],
                                         span = span, iter = 0L)
  }

  char_thresh
}
