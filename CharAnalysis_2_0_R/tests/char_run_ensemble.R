# char_run_ensemble.R
# ---------------------------------------------------------------------
# Propagate chronological uncertainty through CharAnalysis by running
# the pipeline once per chronology in a 1000-member ensemble.
#
# Inputs
#   params_file  : path to *_charParams.csv
#   chron_file   : path to *_chronologies.csv
#                  col 1 = Sample_cm (depth, cm)
#                  cols 2:1001 = CalAge_1:CalAge_1000 (cal yr BP)
#
# Output: list `ensemble` with
#   peaks_matrix    n_steps x n_iter  -- peaksFinal (0/1) per iteration
#   charPeak_matrix n_steps x n_iter  -- C_peak per iteration
#   charBkg_matrix  n_steps x n_iter  -- C_background per iteration
#   prob_peak       n_steps vector    -- fraction of iterations with peak
#   age_grid        n_steps vector    -- common time axis (cal yr BP)
#   iter_age_grids  list (n_iter)     -- ybpI from each iteration (diagnostic)
#   n_iter          integer
#   params          params list used
#
# Notes
#   - yrInterp is fixed from a single reference run (median chronology) so
#     all iterations use an identical resampling interval.
#   - 15 charData depths have no exact match in the chronology matrix
#     (clustered at ~349 cm and ~446 cm); ages at those depths are linearly
#     interpolated from neighboring chron_matrix depths.  A warning is issued
#     at startup listing those depths.  The chronology matrix should ideally
#     cover all charData depths exactly.
#   - chron_matrix rows beyond the charData depth range are ignored.
# ---------------------------------------------------------------------

library(CharAnalysis)

setwd("C:/Users/philip.higuera/OneDrive - The University of Montana/1_phiguera/1_working/CharAnalysis")

# ---- paths ----------------------------------------------------------
params_file <- "CharAnalysis_2_0_R/inst/validation/CH10_charParams.csv"
chron_file  <- "CharAnalysis_2_0_R/tests/CH10_MCAgeDepth_1000_chronologies_2026_06_12.csv"


# ---- read inputs ----------------------------------------------------
message("Reading params and base charcoal data...")
params    <- char_parameters(params_file)
base_data <- params$char_data   # numeric matrix: cmTop, cmBot, ageTop, ageBot, charVol, charCount

message("Reading chronology matrix...")
chron_mat    <- read.csv(chron_file)   # Sample_cm | CalAge_1 ... CalAge_1000
chron_depths <- chron_mat[, 1]
n_iter       <- ncol(chron_mat) - 1L  # 1000


# ---- depth alignment ------------------------------------------------
# charData has two depth columns: cmTop (col 1) and cmBot (col 2).
# For each iteration we need an age at every cmTop and cmBot depth.
# We look those depths up in chron_mat by exact match on Sample_cm;
# the 15 depths that have no exact match are handled by linear interpolation.

cm_top           <- base_data[, 1L]
cm_bot           <- base_data[, 2L]
all_query_depths <- sort(unique(c(cm_top, cm_bot)))
missing_depths   <- all_query_depths[!all_query_depths %in% chron_depths]

if (length(missing_depths) > 0L) {
  warning(sprintf(
    paste0("char_run_ensemble: %d charData depth(s) have no exact match in ",
           "the chronology matrix and will be linearly interpolated.\n",
           "  Depths (cm): %s\n",
           "  The chronology matrix should ideally include ages at all ",
           "charData depths."),
    length(missing_depths),
    paste(round(missing_depths, 3L), collapse = ", ")
  ))
}


# ---- reference run: determine stable yrInterp -----------------------
# char_pretreatment() may auto-set yrInterp when it is 0 or missing in
# the params file.  Run once on the base data to capture the value, then
# fix it so every iteration uses the same resampling interval.
message("Reference run to determine yrInterp...")
pre_ref  <- suppressMessages(
  CharAnalysis:::char_pretreatment(base_data, params$site,
                                   params$pretreatment, params$results, plot_data = 0L)
)
yrInterp <- pre_ref$pretreatment$yrInterp
message(sprintf("  yrInterp = %g yr", yrInterp))

# Fix yrInterp so it is not re-derived from per-iteration chronologies
params$pretreatment$yrInterp <- yrInterp


# ---- common time grid -----------------------------------------------
# Define once from the median youngest and oldest ages across all
# 1000 chronologies, snapped to the nearest yrInterp boundary.
first_ages <- as.numeric(chron_mat[1L,         -1L])  # youngest samples
last_ages  <- as.numeric(chron_mat[nrow(chron_mat), -1L])  # oldest samples

start_age   <- round(median(first_ages) / yrInterp) * yrInterp
end_age     <- round(median(last_ages)  / yrInterp) * yrInterp
common_grid <- seq(start_age, end_age, by = yrInterp)
n_steps     <- length(common_grid)

message(sprintf(
  "Common grid: %.0f to %.0f cal yr BP  |  step = %g yr  |  %d steps",
  start_age, end_age, yrInterp, n_steps
))


# ---- pre-allocate output matrices -----------------------------------
peaks_matrix    <- matrix(NA_real_, nrow = n_steps, ncol = n_iter)
charPeak_matrix <- matrix(NA_real_, nrow = n_steps, ncol = n_iter)
charBkg_matrix  <- matrix(NA_real_, nrow = n_steps, ncol = n_iter)
iter_age_grids  <- vector("list", n_iter)


# ---- helper: snap iteration grid onto common_grid -------------------
# ybpI    : resampled age vector from one iteration
# values  : parallel values vector (same length as ybpI)
# Returns a length(common_grid) vector; positions outside the iteration
# range stay NA.
snap_to_grid <- function(ybpI, values, grid, step) {
  out  <- rep(NA_real_, length(grid))
  idx  <- round((ybpI - grid[1L]) / step) + 1L
  keep <- idx >= 1L & idx <= length(grid)
  out[idx[keep]] <- values[keep]
  out
}


# ---- main loop ------------------------------------------------------
message(sprintf("Running %d iterations...", n_iter))
t0 <- proc.time()[["elapsed"]]

for (i in seq_len(n_iter)) {

  # -- build per-iteration age vectors via depth lookup + interpolation
  # approx() interpolates chron ages at every cmTop/cmBot depth;
  # the 15 missing depths are handled here transparently.
  age_lookup <- approx(
    x      = chron_depths,
    y      = chron_mat[, i + 1L],
    xout   = all_query_depths,
    method = "linear",
    rule   = 2L          # extrapolate at extremes if needed (rare)
  )

  # Map interpolated ages back to data rows by position in all_query_depths
  idx_top <- match(cm_top, all_query_depths)
  idx_bot <- match(cm_bot, all_query_depths)

  data_i       <- base_data
  data_i[, 3L] <- age_lookup$y[idx_top]   # ageTop
  data_i[, 4L] <- age_lookup$y[idx_bot]   # ageBot

  # -- run pipeline (messages suppressed) ----------------------------
  pre <- suppressMessages(
    CharAnalysis:::char_pretreatment(data_i, params$site,
                                     params$pretreatment, params$results, plot_data = 0L)
  )

  charcoal <- suppressMessages(
    CharAnalysis:::char_smooth(pre$charcoal, pre$pretreatment,
                               params$smoothing, params$results, plot_data = 0L)
  )

  # C_peak: residuals (cPeak == 1) or ratios (cPeak == 2)
  charcoal$peak <- if (params$peak_analysis$cPeak == 1L) {
    charcoal$accI - charcoal$accIS
  } else {
    charcoal$accI / charcoal$accIS
  }

  char_thresh <- suppressMessages(
    if (params$peak_analysis$threshType == 1L) {
      CharAnalysis:::char_thresh_global(charcoal, pre$pretreatment,
                                        params$peak_analysis, params$site, params$results,
                                        plot_data = 0L, bkg_sens_in = 0L)
    } else {
      CharAnalysis:::char_thresh_local(charcoal, params$smoothing,
                                       params$peak_analysis, params$site, params$results,
                                       plot_data = 0L)
    }
  )

  peak_result <- suppressMessages(
    CharAnalysis:::char_peak_id(charcoal, pre$pretreatment,
                                params$peak_analysis, char_thresh)
  )
  charcoal    <- peak_result$charcoal
  char_thresh <- peak_result$char_thresh

  post_result <- suppressMessages(
    CharAnalysis:::char_post_process(charcoal, pre$pretreatment,
                                     params$peak_analysis, char_thresh, params$smoothing)
  )
  charcoal <- post_result$charcoal

  # -- extract outputs -----------------------------------------------
  ybpI        <- charcoal$ybpI
  peaks_final <- charcoal$charPeaks[, ncol(charcoal$charPeaks)]  # peaksFinal
  c_peak      <- charcoal$peak
  c_bkg       <- charcoal$accIS

  iter_age_grids[[i]] <- ybpI

  # -- align to common grid ------------------------------------------
  peaks_matrix[, i]    <- snap_to_grid(ybpI, peaks_final, common_grid, yrInterp)
  charPeak_matrix[, i] <- snap_to_grid(ybpI, c_peak,      common_grid, yrInterp)
  charBkg_matrix[, i]  <- snap_to_grid(ybpI, c_bkg,       common_grid, yrInterp)

  if (i %% 100L == 0L) {
    elapsed <- proc.time()[["elapsed"]] - t0
    message(sprintf("  Iteration %4d / %d  (%.0f s elapsed, ~%.0f s remaining)",
                    i, n_iter, elapsed, elapsed / i * (n_iter - i)))
  }
}

message(sprintf("Done. Total time: %.0f s", proc.time()[["elapsed"]] - t0))


# ---- assemble output ------------------------------------------------
ensemble <- list(
  peaks_matrix    = peaks_matrix,
  charPeak_matrix = charPeak_matrix,
  charBkg_matrix  = charBkg_matrix,
  prob_peak       = rowMeans(peaks_matrix, na.rm = TRUE),
  age_grid        = common_grid,
  iter_age_grids  = iter_age_grids,
  n_iter          = n_iter,
  params          = params
)
class(ensemble) <- c("CharEnsemble", "list")

message(sprintf("ensemble$peaks_matrix : %d steps x %d iterations",
                nrow(ensemble$peaks_matrix), ncol(ensemble$peaks_matrix)))
message(sprintf("ensemble$prob_peak    : range [%.3f, %.3f]",
                min(ensemble$prob_peak, na.rm = TRUE),
                max(ensemble$prob_peak, na.rm = TRUE)))


# ---- save results ---------------------------------------------------
out_path <- "CharAnalysis_2_0_R/tests/CH10_ensemble_results.rds"
saveRDS(ensemble, out_path)
message(sprintf("Saved ensemble to: %s", out_path))
