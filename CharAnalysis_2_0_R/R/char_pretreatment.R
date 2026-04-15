#' Interpolate and pretreat a charcoal time series
#'
#' Resamples the raw charcoal data to equal time steps using a
#' proportion-weighted scheme, fills missing-value gaps by linear
#' interpolation, computes charcoal accumulation rates (CHAR), and optionally
#' applies a log transformation.
#'
#' Mirrors \code{CharPretreatment.m} from the MATLAB v2.0 codebase.  The
#' vectorised proportion matrix (four broadcast cases) produces results
#' numerically identical to the MATLAB implementation within floating-point
#' tolerance (~1e-14).
#'
#' @param char_data    Numeric matrix (n x 6+): cmTop, cmBot, ageTop, ageBot,
#'   charVol, charCount. Rows sorted youngest-first (ascending ageTop).
#' @param site         Character string; site name used in diagnostic messages.
#' @param pretreatment Named list with elements:
#'   \describe{
#'     \item{zoneDiv}{Numeric vector of zone boundaries (cal. yr BP),
#'       strictly ascending, at least 2 values.}
#'     \item{yrInterp}{Resampling interval (yr).  0 = use median raw
#'       resolution (auto).}
#'     \item{transform}{Integer: 0 = none; 1 = log10(x+1); 2 = ln(x+1).}
#'   }
#' @param results      Named list; only \code{allFigures} is referenced.
#'   Pass \code{NULL} to suppress any plot-related behaviour (no plots are
#'   produced in the R implementation in Phase 1).
#' @param plot_data    0/1 integer flag.  Ignored in R (no diagnostic plots);
#'   included for API symmetry with the MATLAB function signature.
#'
#' @return Named list with three elements:
#'   \describe{
#'     \item{charcoal}{List of raw and resampled series:
#'       \code{cm, count, vol, con, ybp, acc} (raw);
#'       \code{cmI, countI, volI, conI, accI, ybpI} (resampled).}
#'     \item{pretreatment}{Input list returned with \code{yrInterp} updated
#'       when it was 0 (auto), and \code{zoneDiv[end]} corrected if it
#'       exceeded the bottom age of the last raw sample.}
#'     \item{gap_in}{Integer matrix (nGaps x 2) of gap row-index pairs, or
#'       \code{matrix(NA_integer_, 0, 2)} when no gaps exist.}
#'   }
#'
#' @details
#'   ## Proportion matrix
#'   For each resampled interval \eqn{i} and each raw sample \eqn{j},
#'   \code{prop_matrix[i,j]} is the fraction of the raw sample's age span
#'   that falls within the resampled interval.  Four mutually exclusive
#'   overlap geometries (Cases A-D) are evaluated via matrix broadcasting
#'   across the full \eqn{[N_{rs} \times N_{raw}]} grid -- no loops required.
#'
#'   | Case | Geometry                                      | Overlap            |
#'   |------|-----------------------------------------------|--------------------|
#'   | A    | Raw straddles the **bottom** edge             | rsAgeBot - ageTop  |
#'   | B    | Raw straddles the **top** edge                | ageBot - rsAgeTop  |
#'   | C    | Raw lies **entirely within** resampled interval | ageBot - ageTop  |
#'   | D    | Resampled lies **entirely within** raw sample | yrInterp           |
#'
#'   ## zoneDiv auto-correction
#'   If \code{zoneDiv[end]} exceeds the bottom age of the last raw sample,
#'   the terminal resampled intervals would have no overlapping raw data and
#'   \code{accI} would be \code{NA}.  These \code{NA}s propagate into
#'   \code{charBkg} and can hang the GMM in Phase 2.  The value is silently
#'   corrected to \code{lastAgeBotInData} and the user is notified
#'   (v2.0 behaviour, preserved here).
#'
#' @seealso [char_parameters()], [char_validate_params()], [CharAnalysis()]
char_pretreatment <- function(char_data, site, pretreatment,
                               results = NULL, plot_data = 1L) {

  zone_div  <- pretreatment$zoneDiv
  yr_interp <- pretreatment$yrInterp
  transform <- pretreatment$transform

  # =========================================================================
  # TRIM RECORD TO zoneDiv BOUNDS
  # =========================================================================
  if (length(zone_div) > 1L) {
    ybp_stop  <- zone_div[1L]
    ybp_start <- zone_div[length(zone_div)]
    keep      <- char_data[, 3L] >= ybp_stop & char_data[, 4L] <= ybp_start
    char_data <- char_data[keep, , drop = FALSE]
  }

  # =========================================================================
  # AUTO-CORRECT zoneDiv[end] IF IT EXCEEDS THE DATA BOUNDARY
  #
  #   If zoneDiv[end] extends beyond the bottom age of the last raw sample,
  #   the vectorised proportion matrix assigns NA to those terminal resampled
  #   intervals.  These NAs propagate into charBkg and can hang the GMM.
  #   Correct silently and notify the user (v2.0 behaviour).
  # =========================================================================
  last_age_bot <- char_data[nrow(char_data), 4L]
  n_zones      <- length(pretreatment$zoneDiv)

  if (pretreatment$zoneDiv[n_zones] > last_age_bot) {
    message("NOTE: zoneDiv[end] (", pretreatment$zoneDiv[n_zones],
            " yr BP) exceeds the bottom age of the last raw sample (",
            last_age_bot, " yr BP). zoneDiv[end] corrected to ",
            last_age_bot, " yr BP.")
    pretreatment$zoneDiv[n_zones] <- last_age_bot
    zone_div <- pretreatment$zoneDiv
  }

  # =========================================================================
  # SCREEN FOR MISSING VALUES (sample volume <= 0)
  #
  #   nGaps is defined unconditionally here so it is always safe to reference
  #   (v2.0 bug fix -- in v1.1 it was undefined when nMissingValues == 0).
  # =========================================================================
  missing_idx  <- which(char_data[, 5L] <= 0)
  n_missing    <- length(missing_idx)
  n_gaps       <- 0L
  gap_in       <- matrix(NA_integer_, nrow = 0L, ncol = 2L)

  if (n_missing > 0L) {

    # Identify contiguous runs: a new run starts whenever consecutive missing
    # indices are non-adjacent (diff > 1).  Prepending 99 makes the first
    # missing index always qualify as a run start (99 > 1 is TRUE).
    diff_idx <- c(99L, diff(missing_idx))
    start_in <- missing_idx[diff_idx > 1L]

    # A run ends at positions where the jump to the next missing index > 1,
    # plus always at the very last missing index.
    diff_end <- diff(missing_idx)
    end_in   <- c(missing_idx[diff_end > 1L], missing_idx[n_missing])

    gap_in <- cbind(start_in, end_in)
    n_gaps <- nrow(gap_in)

    cm_gaps <- sum(char_data[gap_in[, 2L] + 1L, 1L] -
                   char_data[gap_in[, 1L] - 1L, 1L])
    yr_gaps <- sum(char_data[gap_in[, 2L] + 1L, 3L] -
                   char_data[gap_in[, 1L] - 1L, 3L])

    message("NOTE: ", n_missing, " missing value(s); ", n_gaps,
            " gap(s) totalling ", cm_gaps, " cm and ",
            round(yr_gaps), " years.")
    message("      Values created via interpolation.")

    for (g in seq_len(n_gaps)) {
      i1 <- gap_in[g, 1L]
      i2 <- gap_in[g, 2L]

      # Anchor depths and volumes from the rows surrounding the gap
      x_anchor <- c(char_data[i1 - 1L, 1L], char_data[i2 + 1L, 1L])
      cm_step  <- diff(char_data[i1, 1:2])
      xi       <- seq(char_data[i1 - 1L, 1L], char_data[i2 + 1L, 1L],
                      by = cm_step)

      y_vol <- c(char_data[i1 - 1L, 5L], char_data[i2 + 1L, 5L])
      y_cnt <- c(char_data[i1 - 1L, 6L], char_data[i2 + 1L, 6L])

      yi_vol <- approx(x_anchor, y_vol, xout = xi)$y
      yi_cnt <- approx(x_anchor, y_cnt, xout = xi)$y

      # Trim if the interpolated sequence is longer than the gap
      # (can happen with variable sampling intervals around the gap;
      # mirrors the MATLAB trim logic exactly)
      gap_len <- length(i1:i2)
      inner_vol <- yi_vol[-c(1L, length(yi_vol))]
      inner_cnt <- yi_cnt[-c(1L, length(yi_cnt))]
      if (length(inner_vol) != gap_len) {
        inner_vol <- inner_vol[-2L]
        inner_cnt <- inner_cnt[-2L]
      }

      char_data[i1:i2, 5L] <- inner_vol
      char_data[i1:i2, 6L] <- inner_cnt
    }
  }

  # Non-contiguous sample warning (informational only; does not affect output)
  remaining_gap <- which((char_data[-1L, 1L] -
                           char_data[-nrow(char_data), 2L]) != 0)
  if (length(remaining_gap) > 0L) {
    rem_cm <- sum(char_data[remaining_gap + 1L, 1L] -
                  char_data[remaining_gap, 1L])
    rem_yr <- sum(char_data[remaining_gap + 1L, 3L] -
                  char_data[remaining_gap, 3L])
    message("NOTE: ", length(remaining_gap), " gap(s) in record totalling ",
            rem_cm, " cm and ", rem_yr, " yr.")
    message("      Resampling will occur across this gap; consider filling")
    message("      missing rows with -999 to interpolate over it instead.")
  }

  # =========================================================================
  # EXTRACT RAW SERIES
  # =========================================================================
  charcoal <- list(
    cm    = char_data[, 1L],
    count = char_data[, 6L],
    vol   = char_data[, 5L],
    con   = char_data[, 6L] / char_data[, 5L],
    ybp   = char_data[, 3L]
  )

  # =========================================================================
  # SEDIMENT ACCUMULATION RATE
  #
  #   Mirrors MATLAB: sedAcc = zeros(n,1); loop fills positions 1:(n-1).
  #   The last element therefore stays 0, which is reproduced here via
  #   c(diff/diff, 0).
  # =========================================================================
  sed_acc <- c(diff(charcoal$cm) / diff(charcoal$ybp), 0)

  # =========================================================================
  # AUTO yrInterp (if not user-specified)
  # =========================================================================
  if (pretreatment$yrInterp == 0) {
    yr_interp            <- round(median(diff(char_data[, 3L])))
    pretreatment$yrInterp <- yr_interp
    message("\nRecord resampled to median resolution of ",
            pretreatment$yrInterp, " years.")
  }

  # =========================================================================
  # BUILD RESAMPLED AGE VECTOR
  #
  #   Mirrors MATLAB colon operator:
  #     zoneDiv(1) : yrInterp : zoneDiv(end)
  # =========================================================================
  charcoal$ybpI <- seq(from = pretreatment$zoneDiv[1L],
                       to   = pretreatment$zoneDiv[length(pretreatment$zoneDiv)],
                       by   = pretreatment$yrInterp)

  age_top <- char_data[, 3L]
  age_bot <- char_data[, 4L]
  N_rs    <- length(charcoal$ybpI)
  N_raw   <- length(age_top)

  rs_age_top <- charcoal$ybpI
  rs_age_bot <- rs_age_top + pretreatment$yrInterp

  # =========================================================================
  # VECTORISED PROPORTION MATRIX  [N_rs x N_raw]
  #
  #   For each pair (resampled interval i, raw sample j), prop_matrix[i,j]
  #   is the fraction of raw sample j's age span that falls within the
  #   resampled interval i.
  #
  #   Broadcasting: both rs_* vectors ([N_rs]) and age_* vectors ([N_raw])
  #   are expanded to [N_rs x N_raw] matrices before comparison.
  #
  #     rsT[i,j] = rs_age_top[i]   (same value across each row)
  #     rsB[i,j] = rs_age_bot[i]
  #     aT[i,j]  = age_top[j]      (same value down each column)
  #     aB[i,j]  = age_bot[j]
  #
  #   Case A: raw sample straddles the BOTTOM edge of the resampled interval
  #           (ageTop inside, ageBot outside below)
  #           overlap = rsAgeBot[i] - ageTop[j]
  #
  #   Case B: raw sample straddles the TOP edge of the resampled interval
  #           (ageBot inside, ageTop outside above)
  #           overlap = ageBot[j] - rsAgeTop[i]
  #
  #   Case C: raw sample lies ENTIRELY WITHIN the resampled interval
  #           overlap = ageBot[j] - ageTop[j]
  #
  #   Case D: resampled interval lies ENTIRELY WITHIN the raw sample
  #           overlap = yrInterp
  #
  #   Note: Cases C and A share (aT >= rsT); Case C is the stricter subset
  #   (aB <= rsB), while Case A requires (aB > rsB).  The four cases are
  #   mutually exclusive and exhaustive for all overlapping pairs.
  # =========================================================================

  # Expand to [N_rs x N_raw] matrices
  rsT <- matrix(rs_age_top, nrow = N_rs, ncol = N_raw)            # rs_age_top[i]
  rsB <- matrix(rs_age_bot, nrow = N_rs, ncol = N_raw)            # rs_age_bot[i]
  aT  <- matrix(age_top,   nrow = N_rs, ncol = N_raw, byrow = TRUE) # age_top[j]
  aB  <- matrix(age_bot,   nrow = N_rs, ncol = N_raw, byrow = TRUE) # age_bot[j]

  caseA <- (aT >= rsT) & (aT  < rsB) & (aB >  rsB)   # straddles bottom edge
  caseB <- (aT  < rsT) & (aB <= rsB) & (aB >  rsT)   # straddles top edge
  caseC <- (aT >= rsT) & (aB <= rsB)                  # entirely within
  caseD <- (aT  < rsT) & (aB  > rsB)                  # raw contains interval

  prop_matrix <- caseA * (rsB - aT)              +
                 caseB * (aB  - rsT)             +
                 caseC * (aB  - aT)              +
                 caseD *  pretreatment$yrInterp

  prop_matrix <- prop_matrix / pretreatment$yrInterp   # normalise to fraction

  # =========================================================================
  # DERIVE RESAMPLED VALUES
  #
  #   For each resampled interval i, sum the contributions from all raw
  #   samples j where prop_matrix[i,j] > 0.  This mirrors the MATLAB loop:
  #     in = propMatrix(i,:) > 0;
  #     Charcoal.countI(i) = sum( count(in) .* propMatrix(i,in)' );
  # =========================================================================
  charcoal$countI <- rep(NA_real_, N_rs)
  charcoal$volI   <- rep(NA_real_, N_rs)
  charcoal$conI   <- rep(NA_real_, N_rs)
  sed_acc_I       <- rep(NA_real_, N_rs)

  for (i in seq_len(N_rs)) {
    in_idx <- prop_matrix[i, ] > 0
    if (any(in_idx)) {
      p <- prop_matrix[i, in_idx]
      charcoal$countI[i] <- sum(charcoal$count[in_idx] * p)
      charcoal$volI[i]   <- sum(charcoal$vol[in_idx]   * p)
      charcoal$conI[i]   <- sum(charcoal$con[in_idx]   * p)
      sed_acc_I[i]       <- sum(sed_acc[in_idx]          * p)
    }
  }

  # Resampled depths: linear interpolation of raw cm vs. age.
  # rule = 1: return NA outside the data range, matching MATLAB's interp1
  # default (which returns NaN for out-of-range queries).
  charcoal$cmI <- approx(charcoal$ybp, charcoal$cm,
                          xout = charcoal$ybpI, rule = 1L)$y

  # =========================================================================
  # CHARCOAL ACCUMULATION RATES
  # =========================================================================
  charcoal$acc  <- charcoal$con  * sed_acc     # raw CHAR (# cm-2 yr-1)
  charcoal$accI <- charcoal$conI * sed_acc_I   # resampled CHAR

  # =========================================================================
  # OPTIONAL LOG TRANSFORM
  # =========================================================================
  if (transform == 1L) {
    charcoal$accI <- log10(charcoal$accI + 1)
  }
  if (transform == 2L) {
    charcoal$accI <- log(charcoal$accI + 1)
  }

  # =========================================================================
  # RETURN
  # =========================================================================
  list(
    charcoal     = charcoal,
    pretreatment = pretreatment,
    gap_in       = gap_in
  )
}
