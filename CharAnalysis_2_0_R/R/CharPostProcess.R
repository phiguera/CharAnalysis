#' Post-processing: FRIs, fire frequency, Weibull statistics, output matrix
#'
#' Pure-computation post-processing step that follows peak identification.
#' Mirrors \code{CharPostProcess.m} and \code{smoothFRI.m} from the MATLAB
#' v2.0 codebase.  No figures are produced here; all computed values are
#' returned for downstream plotting and file output.
#'
#' @param charcoal     Named list from [char_peak_id()] containing (among
#'   others): \code{peak}, \code{ybpI}, \code{cmI}, \code{accI},
#'   \code{accIS}, \code{countI}, \code{volI}, \code{conI},
#'   \code{charPeaks} (\eqn{[N \times T]}), \code{charPeaksThresh}.
#' @param pretreatment Named list with \code{yrInterp} and \code{zoneDiv}.
#' @param peak_analysis Named list with \code{threshType}, \code{threshValues},
#'   \code{minCountP}, and \code{peakFrequ}.
#' @param char_thresh  Named list from [char_thresh_local()] or
#'   [char_thresh_global()] containing \code{pos}, \code{neg}, \code{SNI},
#'   \code{GOF}, and \code{minCountP}.
#' @param smoothing    Named list with \code{yr} (used for smoothed FRI span).
#'
#' @return Named list with three components:
#'   \describe{
#'     \item{charcoal}{Input \code{charcoal} augmented with:
#'       \code{peakInsig} (0/1, peaks failing min-count screen),
#'       \code{peakMagnitude} (integrated C_peak above threshold, per peak),
#'       \code{smoothedFireFrequ} (Lowess-smoothed fire frequency, peaks/ka),
#'       \code{peaksFrequ} (raw sliding-window fire frequency, peaks/ka).}
#'     \item{post}{Named list of intermediate results for plotting:
#'       \code{peakIn}, \code{peakScreenIn}, \code{CharcoalCharPeaks},
#'       \code{threshIn} (global only), \code{peak_mag}, \code{ff_sm},
#'       \code{FRIyr}, \code{FRI}, \code{smFRIyr}, \code{smFRI},
#'       \code{smFRIci}, \code{yis}, \code{FRI_params_zone}, \code{alpha},
#'       \code{nBoot}.}
#'     \item{char_results}{Numeric matrix (\eqn{N \times 33}): the full
#'       output table written to CSV by [CharWriteResults()].  Column
#'       order matches the MATLAB \code{charResults} matrix exactly.}
#'   }
#'
#' @details
#'   ## Smoothed FRI (\code{smooth_fri} helper)
#'   A sliding window of width \code{peak_analysis$peakFrequ} yr steps every
#'   100 yr from -70 to the maximum age.  Within each window the mean FRI and
#'   bootstrapped (1-\code{alpha})x100 % CIs are computed, then the mean and
#'   CI series are smoothed with \code{char_lowess}.  Mirrors
#'   \code{smoothFRI.m}.
#'
#'   ## Weibull statistics
#'   FRIs within each zone (\code{pretreatment$zoneDiv}) are binned (bin width
#'   20 yr, edges 20-1000 yr) and a two-parameter Weibull is fitted to the bin
#'   centres weighted by bin frequency, matching MATLAB's
#'   \code{wblfit(centres, [], [], freq)}.  Bootstrap CIs use
#'   \code{nBoot = 100} resamples of the raw FRIs.  Zones with \eqn{\leq 1}
#'   FRI or max FRI > 5000 yr are left as -999.
#'
#' @seealso [char_peak_id()], [CharAnalysis()], [CharWriteResults()]
#'
#' Requires: char_lowess() from utils_lowess.R,
#'           MASS package for fitdistr()

# ---------------------------------------------------------------------------
# HELPER: smooth_fri -- R port of smoothFRI.m
# ---------------------------------------------------------------------------
smooth_fri <- function(yr, peaks, win_width,
                        alpha    = 0.05,
                        n_boot   = 100L,
                        fri_param = 1L) {

  step_length <- 100L
  max_fris    <- 100L

  steps   <- seq(-70, max(yr), by = step_length)
  n_steps <- length(steps)
  peak_yr <- yr[peaks == 1]

  if (length(peak_yr) < 2L) {
    return(list(FRIyr   = NA_real_, FRI     = NA_real_,
                smFRIyr = NA_real_, smFRI   = NA_real_,
                smFRIci = NA_real_))
  }

  params_mat <- matrix(NA_real_, 3L, n_steps)   # rows: mean, upper-CI, lower-CI

  for (i in seq_len(n_steps)) {
    start_yr <- steps[i] - 0.5 * win_width
    end_yr   <- steps[i] + 0.5 * win_width
    yr_in    <- which(peak_yr >= start_yr & peak_yr < end_yr)
    fris     <- diff(peak_yr[yr_in])

    if (length(fris) > 1L) {
      if (fri_param == 0L) {
        # Standard-deviation bands
        params_mat[1L, i] <- mean(fris)
        se <- stats::sd(fris) / sqrt(length(fris))
        params_mat[2L, i] <- mean(fris) + se
        params_mat[3L, i] <- mean(fris) - se
      } else {
        # Bootstrapped (1-alpha)*100 % CIs
        params_mat[1L, i] <- mean(fris)
        boot_means <- vapply(seq_len(n_boot), function(b) {
          mean(fris[sample.int(length(fris), length(fris), replace = TRUE)])
        }, numeric(1L))
        qs <- stats::quantile(boot_means,
                              probs = c(alpha / 2, 1 - alpha / 2),
                              names = FALSE)
        params_mat[2L, i] <- qs[2L]   # upper CI
        params_mat[3L, i] <- qs[1L]   # lower CI
      }
    }
  }

  # ---- Smooth FRI parameter series with Lowess ----
  in_idx <- which(!is.na(params_mat[1L, ]))
  if (length(in_idx) < 3L) {
    return(list(FRIyr   = peak_yr[-length(peak_yr)],
                FRI     = diff(peak_yr),
                smFRIyr = NA_real_, smFRI   = NA_real_,
                smFRIci = NA_real_))
  }

  x        <- steps[in_idx]
  smooth_in <- (win_width / step_length) / length(x)  # span as fraction

  y  <- char_lowess(params_mat[1L, in_idx], span = smooth_in, iter = 0L)
  y3 <- char_lowess(params_mat[2L, in_idx], span = smooth_in, iter = 0L) # upper
  y2 <- char_lowess(params_mat[3L, in_idx], span = smooth_in, iter = 0L) # lower

  list(
    FRIyr   = peak_yr[-length(peak_yr)],
    FRI     = diff(peak_yr),
    smFRIyr = x,
    smFRI   = y,
    smFRIci = cbind(y3, y2)   # [upper, lower]  matches MATLAB [y3, y2]
  )
}

# ---------------------------------------------------------------------------
# MAIN FUNCTION
# ---------------------------------------------------------------------------
char_post_process <- function(charcoal, pretreatment, peak_analysis,
                               char_thresh, smoothing) {

  r        <- pretreatment$yrInterp
  zone_div <- pretreatment$zoneDiv
  N        <- length(charcoal$ybpI)
  T_thresh <- ncol(charcoal$charPeaks)
  alpha    <- 0.05
  n_boot   <- 100L

  # =========================================================================
  # 1. Resolve peak / threshold index vectors
  # Mirrors CharPostProcess.m lines 48-64.
  # =========================================================================
  if (peak_analysis$threshType == 1L) {
    # Global threshold: identify which column of charPeaks corresponds to
    # each threshValue by matching mean threshold value at peak positions.
    thresh_in1 <- colSums(charcoal$charPeaksThresh) /
                  pmax(colSums(charcoal$charPeaks), 1)  # avoid div-by-zero
    thresh_in <- vapply(seq_len(T_thresh), function(j) {
      min(which(thresh_in1 >= char_thresh$pos[1L, j]))
    }, integer(1L))
    peak_in           <- which(charcoal$charPeaks[, thresh_in[T_thresh]] > 0)
    charcoal_charpeaks <- charcoal$charPeaks[, thresh_in, drop = FALSE]
    peak_screen_in    <- which(!is.na(char_thresh$minCountP[, thresh_in[T_thresh]]) &
                               char_thresh$minCountP[, thresh_in[T_thresh]] >
                               peak_analysis$minCountP)
  } else {
    # Local threshold: use all columns directly
    thresh_in         <- integer(0L)   # not used for local
    peak_in           <- which(charcoal$charPeaks[, T_thresh] > 0)
    charcoal_charpeaks <- charcoal$charPeaks
    peak_screen_in    <- which(!is.na(char_thresh$minCountP[, T_thresh]) &
                               char_thresh$minCountP[, T_thresh] >
                               peak_analysis$minCountP)
  }

  # =========================================================================
  # 2. Peak magnitude
  # For each contiguous run of samples where C_peak exceeds the final
  # threshold, accumulate (sum x resolution) to get pieces cm^-^2 peak^-^1.
  # Mirrors CharPostProcess.m lines 67-103.
  # =========================================================================
  c_peaks <- charcoal$peak - char_thresh$pos[, T_thresh]
  c_peaks[c_peaks < 0] <- 0
  c_peaks[N]           <- 0    # oldest sample never a peak (index N = oldest)

  peak_in_mat <- matrix(0L, N, 2L)
  peak_mag    <- numeric(N)

  for (idx in seq_len(N)) {
    if (idx == 1L && c_peaks[idx] > 0) {
      peak_in_mat[idx, 1L] <- idx
      step <- 1L
      while (idx + step <= N && c_peaks[idx + step] > 0) {
        peak_in_mat[idx, 2L] <- idx + step
        step <- step + 1L
      }
    } else if (idx > 1L && c_peaks[idx] > 0 && c_peaks[idx - 1L] == 0) {
      peak_in_mat[idx, 1L] <- idx
      step <- 1L
      while (idx + step <= N && c_peaks[idx + step] > 0) {
        peak_in_mat[idx, 2L] <- idx + step
        step <- step + 1L
      }
    }

    if (peak_in_mat[idx, 1L] > 0 && peak_in_mat[idx, 2L] == 0)
      peak_in_mat[idx, 2L] <- peak_in_mat[idx, 1L]

    if (peak_in_mat[idx, 2L] > 0) {
      p1 <- peak_in_mat[idx, 1L]; p2 <- peak_in_mat[idx, 2L]
      peak_mag[p2] <- sum(c_peaks[p1:p2]) *
                      ((p2 - p1) + r)
    }
  }

  # =========================================================================
  # 3. Fire frequency time series
  # Sliding window of peakFrequ years scaled to the window actually used at
  # record edges, then smoothed with Lowess.
  # Mirrors CharPostProcess.m lines 106-131.
  # =========================================================================
  ff_yr    <- peak_analysis$peakFrequ
  half_win <- floor(ff_yr / r) / 2

  peaks_frequ <- numeric(N)
  full_win_samples <- round(ff_yr / r)

  for (i in seq_len(N)) {
    if (i < half_win) {
      # Start of record
      win_len  <- half_win + i
      peaks_frequ[i] <- sum(charcoal_charpeaks[seq_len(floor(win_len)),
                                                T_thresh]) *
                         (full_win_samples / floor(win_len))
    } else if (i > N - half_win) {
      # End of record
      win_slice <- charcoal_charpeaks[(i - floor(half_win)):N, T_thresh]
      peaks_frequ[i] <- sum(win_slice) * (ff_yr / r) / length(win_slice)
    } else {
      # Middle of record
      idx_lo <- max(1L,  ceiling(i - 0.5 * (ff_yr / r)) + 1L)
      idx_hi <- min(N,   ceiling(i + 0.5 * round(ff_yr / r)))
      peaks_frequ[i] <- sum(charcoal_charpeaks[idx_lo:idx_hi, T_thresh])
    }
  }

  ff_sm <- char_lowess(peaks_frequ, span = ff_yr / r, iter = 0L)

  # =========================================================================
  # 4. Smoothed FRI curve
  # Mirrors CharPostProcess.m lines 133-151.
  # =========================================================================
  peak_yrs <- charcoal$ybpI[charcoal_charpeaks[, T_thresh] > 0]

  if (length(peak_yrs) > 2L) {
    fri_out <- smooth_fri(charcoal$ybpI,
                          charcoal_charpeaks[, T_thresh],
                          win_width = peak_analysis$peakFrequ,
                          alpha     = alpha,
                          n_boot    = n_boot,
                          fri_param = 1L)
    FRIyr   <- fri_out$FRIyr
    FRI     <- fri_out$FRI
    smFRIyr <- fri_out$smFRIyr
    smFRI   <- fri_out$smFRI
    smFRIci <- fri_out$smFRIci

    # Interpolate smoothed FRI onto ybpI grid (for output column 23)
    if (!all(is.na(smFRI)) && length(smFRI) > 2L) {
      yis_grid <- charcoal$ybpI[charcoal$ybpI < max(smFRIyr)]
      yis <- stats::approx(smFRIyr, smFRI, xout = yis_grid,
                            method = "linear", rule = 1)$y
    } else {
      yis <- rep(-999, length(smFRI))
      warning("char_post_process: fewer than 3 FRIs -- smoothFRI set to -999.")
    }
  } else {
    FRIyr <- NA_real_;  FRI     <- NA_real_
    smFRIyr <- NA_real_; smFRI  <- NA_real_;  smFRIci <- NA_real_
    yis   <- NA_real_
  }

  # =========================================================================
  # 5. Per-zone Weibull / FRI statistics
  # Mirrors CharPostProcess.m lines 153-222.
  #
  # Weibull fit uses MASS::fitdistr on bin-centre-weighted data (equivalent
  # to MATLAB's wblfit(centres, [], [], freq)), with bin width 20 yr and
  # edges 20:20:1000 yr.
  # =========================================================================
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required: install.packages('MASS')")

  n_zones         <- length(zone_div) - 1L
  FRI_params_zone <- matrix(NA_real_, nrow = n_zones, ncol = 10L)
  bin_width       <- 20L
  # bin_edges: 20, 40, ..., 1000  -- 50 values, 49 bins [20,40),[40,60),...,[980,1000)
  # Matches MATLAB histcounts(FRIz, binWidth:binWidth:1000) which has 49 bins.
  bin_edges       <- seq(bin_width, 1000L, by = bin_width)
  bin_centers     <- bin_edges[-length(bin_edges)] + bin_width / 2  # 30,50,...,990

  zone_list <- vector("list", n_zones)

  for (z in seq_len(n_zones)) {
    x_plot <- charcoal$ybpI[charcoal_charpeaks[, T_thresh] > 0]
    x_plot <- x_plot[x_plot >= zone_div[z] & x_plot < zone_div[z + 1L]]
    fri_z  <- diff(x_plot)
    fri_z  <- fri_z[is.finite(fri_z)]   # drop NA/NaN from gaps in ybpI

    if (length(fri_z) <= 1L || max(fri_z) > 5000) next

    # Bin FRIs matching MATLAB histcounts(FRIz, binWidth:binWidth:1000):
    # values outside [bin_edges[1], bin_edges[end]) are silently ignored.
    # R's hist() errors on out-of-range values, so clip first.
    fri_z_bins <- fri_z[fri_z >= bin_edges[1L] &
                         fri_z <  bin_edges[length(bin_edges)]]
    freq_tab <- if (length(fri_z_bins) > 0L) {
      graphics::hist(fri_z_bins,
                     breaks = bin_edges,   # 49 bins [20,40)...[980,1000)
                     plot   = FALSE,
                     right  = FALSE)$counts
    } else {
      integer(length(bin_centers))
    }
    # freq_tab[j] corresponds to bin_centers[j] -- both length 49
    ok <- freq_tab > 0L
    if (sum(ok) < 2L || sum(freq_tab) == 0L) next

    x_rep <- rep(bin_centers[ok], freq_tab[ok])

    # Fit Weibull via MLE.  MASS::fitdistr() can fail with the default
    # starting values when x_rep is small or the likelihood surface is flat
    # (typically < ~10 points or all-unique bin centres).  On failure, retry
    # with method-of-moments starting values before giving up.
    .wbl_start <- function(x) {
      m <- mean(x);  s <- max(stats::sd(x), 1e-6)
      # Approximation: shape ≈ (m/s)^1.086  (valid for CV < ~1)
      sh <- max(0.1, (m / s)^1.086)
      list(shape = sh, scale = m / gamma(1 + 1 / sh))
    }
    fit_wbl <- tryCatch(
      suppressWarnings(MASS::fitdistr(x_rep, "weibull")),
      error = function(e) {
        tryCatch(
          suppressWarnings(MASS::fitdistr(x_rep, "weibull",
                                          start = .wbl_start(x_rep))),
          error = function(e2) NULL
        )
      }
    )
    if (is.null(fit_wbl)) next

    wbl_scale <- unname(fit_wbl$estimate["scale"])   # = MATLAB param(1) = WBLb
    wbl_shape <- unname(fit_wbl$estimate["shape"])   # = MATLAB param(2) = WBLc

    # KS goodness-of-fit against fitted Weibull CDF.
    # Matches MATLAB's kstest(FRIz, [FRIBinKS', wbl_cdf']) where
    # FRIBinKS = 0:20:5000.  MATLAB evaluates the CDF at each data point
    # by linear interpolation from that discrete table, then uses the
    # asymptotic Kolmogorov distribution for the p-value.
    # We replicate both choices: discrete CDF via approxfun, exact = FALSE.
    fri_bin_ks   <- seq(0L, 5000L, by = 20L)
    wbl_cdf_disc <- stats::pweibull(fri_bin_ks,
                                    shape = wbl_shape, scale = wbl_scale)
    disc_cdf_fn  <- stats::approxfun(fri_bin_ks, wbl_cdf_disc,
                                     method = "linear", rule = 2L)
    ks_res <- tryCatch(
      suppressWarnings(stats::ks.test(fri_z, disc_cdf_fn, exact = FALSE)),
      error = function(e) NULL
    )
    p_ks <- if (is.null(ks_res)) NA_real_ else ks_res$p.value

    # Bootstrap CIs on Weibull parameters and mFRI
    mean_mfri_boot <- numeric(n_boot)
    wbl_scale_boot <- numeric(n_boot)
    wbl_shape_boot <- numeric(n_boot)

    for (b in seq_len(n_boot)) {
      fri_t      <- fri_z[sample.int(length(fri_z), length(fri_z), replace = TRUE)]
      fri_t_bins <- fri_t[fri_t >= bin_edges[1L] &
                           fri_t <  bin_edges[length(bin_edges)]]
      freq_t <- if (length(fri_t_bins) > 0L) {
        graphics::hist(fri_t_bins,
                       breaks = bin_edges,
                       plot   = FALSE,
                       right  = FALSE)$counts
      } else {
        integer(length(bin_centers))
      }
      ok_t <- freq_t > 0L
      if (sum(ok_t) < 2L) next
      x_rep_t <- rep(bin_centers[ok_t], freq_t[ok_t])
      fit_t <- tryCatch(
        suppressWarnings(MASS::fitdistr(x_rep_t, "weibull")),
        error = function(e) NULL
      )
      if (!is.null(fit_t)) {
        wbl_scale_boot[b] <- unname(fit_t$estimate["scale"])
        wbl_shape_boot[b] <- unname(fit_t$estimate["shape"])
      } else {
        wbl_scale_boot[b] <- NA_real_
        wbl_shape_boot[b] <- NA_real_
      }
      mean_mfri_boot[b] <- mean(fri_t)
    }

    # Use na.rm = TRUE: failed Weibull fits store NA_real_, and NA > 0 returns
    # NA (not FALSE) in R, so NAs pass through the > 0 filter without it.
    wbl_scale_ci    <- stats::quantile(wbl_scale_boot[wbl_scale_boot > 0],
                                        probs = c(0.025, 0.975), names = FALSE,
                                        na.rm = TRUE)
    wbl_shape_ci    <- stats::quantile(wbl_shape_boot[wbl_shape_boot > 0],
                                        probs = c(0.025, 0.975), names = FALSE,
                                        na.rm = TRUE)
    mean_mfri_ci    <- stats::quantile(mean_mfri_boot[mean_mfri_boot > 0],
                                        probs = c(0.025, 0.975), names = FALSE,
                                        na.rm = TRUE)

    # NOTE: MATLAB CharPostProcess stores CIs as [quantile(2.5%), quantile(97.5%)]
    # in columns labelled uCI / lCI (i.e. uCI = lower bound, lCI = upper bound).
    # The R code follows the same inverted-label convention for column-exact
    # CSV compatibility with the MATLAB charResults output.
    FRI_params_zone[z, ] <- c(
      length(fri_z),           # 1  nFRI
      mean(fri_z),             # 2  mFRI
      mean_mfri_ci[1L],        # 3  mFRI_uCI  (MATLAB convention: 2.5th percentile)
      mean_mfri_ci[2L],        # 4  mFRI_lCI  (MATLAB convention: 97.5th percentile)
      wbl_scale,               # 5  WBLb  (scale)
      wbl_scale_ci[1L],        # 6  WBLb_uCI  (MATLAB convention: 2.5th percentile)
      wbl_scale_ci[2L],        # 7  WBLb_lCI  (MATLAB convention: 97.5th percentile)
      wbl_shape,               # 8  WBLc  (shape)
      wbl_shape_ci[1L],        # 9  WBLc_uCI  (MATLAB convention: 2.5th percentile)
      wbl_shape_ci[2L]         # 10 WBLc_lCI  (MATLAB convention: 97.5th percentile)
    )

    zone_list[[z]] <- list(
      FRI            = fri_z,
      fri_binned     = fri_z_bins,
      bin_centers    = bin_centers,
      freq           = freq_tab,
      wbl_scale      = wbl_scale,
      wbl_shape      = wbl_shape,
      p_ks           = p_ks,
      wbl_scale_ci   = wbl_scale_ci,
      wbl_shape_ci   = wbl_shape_ci,
      mean_mFRI      = mean(fri_z),
      mean_mFRI_ci   = mean_mfri_ci,
      x_plot         = x_plot
    )
  }

  # =========================================================================
  # 6. peakInsig -- samples where minCountP exceeded alphaPeak before removal
  # Mirrors CharPostProcess.m lines 226-228.
  # =========================================================================
  peak_insig <- integer(N)
  peak_insig[peak_screen_in] <- 1L

  # =========================================================================
  # 7. Augment charcoal struct
  # =========================================================================
  charcoal$peakInsig        <- peak_insig
  charcoal$peakMagnitude    <- peak_mag
  charcoal$smoothedFireFrequ <- ff_sm
  charcoal$peaksFrequ       <- peaks_frequ

  # =========================================================================
  # 8. Assemble output matrix  (N x 33)
  # Column order matches MATLAB charResults exactly:
  #   1  cmI       2  ybpI     3  countI    4  volI      5  conI
  #   6  accI      7  accIS    8  peak      9  pos[,1]  10  pos[,2]
  #  11  pos[,3]  12  pos[,T]  13 neg[,T]  14 SNI      15  GOF
  #  16  peaks[,1] 17 peaks[,2] 18 peaks[,3] 19 peaks[,T]
  #  20  peakInsig 21 peakMag  22 smFireFreq 23 yis (smFRIs)
  #  24-33 FRI_params_zone (1:nZones rows)
  # Mirrors CharPostProcess.m lines 233-251.
  # =========================================================================
  char_results <- matrix(NA_real_, nrow = N, ncol = 33L)

  # Expand scalar SNI to vector (global threshold returns scalar)
  sni_col <- if (length(char_thresh$SNI) == 1L) {
    rep(char_thresh$SNI, N)
  } else {
    char_thresh$SNI
  }

  char_results[, 1:22] <- cbind(
    charcoal$cmI,
    charcoal$ybpI,
    charcoal$countI,
    charcoal$volI,
    charcoal$conI,
    charcoal$accI,
    charcoal$accIS,
    charcoal$peak,
    char_thresh$pos,              # [N x T_thresh] -- 4 columns (9-12)
    char_thresh$neg[, T_thresh],  # final-threshold negative column (13)
    sni_col,                      # 14
    char_thresh$GOF,              # 15  (GOF vector; NA for global)
    charcoal_charpeaks,           # [N x T_thresh] peaks 1-4 (16-19)
    charcoal$peakInsig,           # 20
    charcoal$peakMagnitude,       # 21
    charcoal$smoothedFireFrequ    # 22
  )

  # Column 23: smoothed FRI interpolated onto ybpI grid
  if (!all(is.na(yis))) {
    n_yis <- length(yis)
    char_results[seq_len(n_yis), 23L] <- yis
  }

  # Columns 24-33: per-zone FRI statistics (one row per zone)
  if (n_zones > 0L) {
    char_results[seq_len(n_zones), 24:33] <- FRI_params_zone
  }

  # =========================================================================
  # 9. Pack Post struct
  # =========================================================================
  post <- list(
    peakIn            = peak_in,
    peakScreenIn      = peak_screen_in,
    CharcoalCharPeaks = charcoal_charpeaks,
    threshIn          = thresh_in,
    peak_mag          = peak_mag,
    ff_sm             = ff_sm,
    FRIyr             = FRIyr,
    FRI               = FRI,
    smFRIyr           = smFRIyr,
    smFRI             = smFRI,
    smFRIci           = smFRIci,
    yis               = yis,
    FRI_params_zone   = FRI_params_zone,
    zone              = zone_list,
    alpha             = alpha,
    nBoot             = n_boot,
    char_results      = char_results
  )

  list(
    charcoal     = charcoal,
    post         = post,
    char_results = char_results
  )
}
