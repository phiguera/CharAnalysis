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
#   - charData depths with no exact match in the chronology matrix are
#     linearly interpolated from neighboring depths. A warning lists any
#     such depths at startup. The chronology matrix should ideally cover
#     all charData depths exactly.
#   - chron_matrix rows beyond the charData depth range are ignored.
#
# Parallel processing
#   Set n_cores > 1 to distribute iterations across multiple CPU cores.
#   This works on both Windows and Mac. See the startup message for guidance.
#   NOTE for developers: parallel workers load the *installed* CharAnalysis
#   package, not a devtools::load_all() session. Run devtools::install()
#   before testing the parallel path during development.
# ---------------------------------------------------------------------

library(CharAnalysis)

setwd("C:/Users/philip.higuera/OneDrive - The University of Montana/1_phiguera/1_working/CharAnalysis")

# =====================================================================
# USER SETTINGS
# =====================================================================

params_file <- "CharAnalysis_2_0_R/inst/validation/CH10_charParams.csv"
chron_file  <- "CharAnalysis_2_0_R/tests/CH10_MCAgeDepth_1000_chronologies_2026_06_12.csv"

# Number of CPU cores to use.
# Default: all available cores minus one (leaves one core free for the OS).
# Override by setting n_cores to a specific integer, e.g.: n_cores <- 4
# To run sequentially (no parallel), set: n_cores <- 1
n_cores <- max(1L, parallel::detectCores() - 1L)

# NOTE: parallel workers load the *installed* CharAnalysis package.
# If you have made recent changes to the package source, install first:
#   devtools::install("CharAnalysis_2_0_R", quiet = TRUE, upgrade = FALSE)

# Save results and figure to disk?
# 0 = prompt at end of run (default, interactive use)
# 1 = save automatically without prompting (use for batch/multi-site scripts)
# Both the data files (.rds, .csv) and the figure (.pdf) are controlled
# by this single switch.
save_results <- 0

# =====================================================================


# ---- startup message ------------------------------------------------
avail <- parallel::detectCores()
if (n_cores > avail) {
  warning(sprintf(
    "n_cores (%d) exceeds detected cores (%d); setting n_cores = %d.",
    n_cores, avail, avail
  ))
  n_cores <- avail
}

if (n_cores == 1L) {
  message(paste0(
    "\n-- Chronological uncertainty ensemble ",
    paste(rep("-", 38), collapse = ""), "\n",
    "  Cores      : 1 (sequential)\n",
    "\n",
    "  NOTE: 1000 ensemble runs can take >10-15 minutes on a single\n",
    "  core (benchmark: ~15 min on a standard laptop). Parallel\n",
    "  processing can reduce this substantially -- adding N cores\n",
    "  cuts runtime by approximately 1/N (e.g., 19 cores ~ 2 min).\n",
    "\n",
    "  Only 1 core was detected on this machine. To override, set\n",
    "  n_cores manually at the top of this script.\n",
    paste(rep("-", 58), collapse = ""), "\n"
  ))
} else {
  message(paste0(
    "\n-- Chronological uncertainty ensemble ",
    paste(rep("-", 38), collapse = ""), "\n",
    "  Cores      : ", n_cores, " of ", avail, " available",
    " (parallel, Mac and Windows compatible)\n",
    "  Runtime    : ~1/", n_cores, " of sequential time\n",
    "             : benchmark ~15 min (1 core) vs. ~2 min (19 cores)\n",
    "\n",
    "  NOTE: To adjust the number of cores used, set n_cores at the\n",
    "  top of this script. To run sequentially: n_cores <- 1\n",
    paste(rep("-", 58), collapse = ""), "\n"
  ))
}


# ---- read inputs ----------------------------------------------------
message("Reading params and base charcoal data...")
params    <- char_parameters(params_file)
base_data <- params$char_data

message("Reading chronology matrix...")
chron_mat    <- read.csv(chron_file)
chron_depths <- chron_mat[, 1]
n_iter       <- ncol(chron_mat) - 1L


# ---- depth alignment ------------------------------------------------
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
message("Reference run to determine yrInterp...")
pre_ref  <- suppressMessages(
  CharAnalysis:::char_pretreatment(base_data, params$site,
                                   params$pretreatment, params$results, plot_data = 0L)
)
yrInterp <- pre_ref$pretreatment$yrInterp
message(sprintf("  yrInterp = %g yr", yrInterp))
params$pretreatment$yrInterp <- yrInterp


# ---- common time grid -----------------------------------------------
first_ages <- as.numeric(chron_mat[1L,              -1L])
last_ages  <- as.numeric(chron_mat[nrow(chron_mat), -1L])

start_age   <- round(median(first_ages) / yrInterp) * yrInterp
end_age     <- round(median(last_ages)  / yrInterp) * yrInterp
common_grid <- seq(start_age, end_age, by = yrInterp)
n_steps     <- length(common_grid)

message(sprintf(
  "Common grid: %.0f to %.0f cal yr BP  |  step = %g yr  |  %d steps",
  start_age, end_age, yrInterp, n_steps
))


# ---- helper: snap iteration grid onto common_grid -------------------
snap_to_grid <- function(ybpI, values, grid, step) {
  out  <- rep(NA_real_, length(grid))
  idx  <- round((ybpI - grid[1L]) / step) + 1L
  keep <- idx >= 1L & idx <= length(grid)
  out[idx[keep]] <- values[keep]
  out
}


# ---- worker function ------------------------------------------------
# Encapsulates one iteration; called by both the sequential loop and
# parLapply. All inputs are read-only objects exported to each worker.
run_one_iteration <- function(i) {
  age_lookup <- approx(
    x      = chron_depths,
    y      = chron_mat[, i + 1L],
    xout   = all_query_depths,
    method = "linear",
    rule   = 2L
  )

  idx_top <- match(cm_top, all_query_depths)
  idx_bot <- match(cm_bot, all_query_depths)

  data_i       <- base_data
  data_i[, 3L] <- age_lookup$y[idx_top]
  data_i[, 4L] <- age_lookup$y[idx_bot]

  pre <- suppressMessages(
    CharAnalysis:::char_pretreatment(data_i, params$site,
                                     params$pretreatment, params$results, plot_data = 0L)
  )

  charcoal <- suppressMessages(
    CharAnalysis:::char_smooth(pre$charcoal, pre$pretreatment,
                               params$smoothing, params$results, plot_data = 0L)
  )

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

  list(
    peaks_final = snap_to_grid(charcoal$ybpI,
                               charcoal$charPeaks[, ncol(charcoal$charPeaks)],
                               common_grid, yrInterp),
    c_peak      = snap_to_grid(charcoal$ybpI, charcoal$peak,   common_grid, yrInterp),
    c_bkg       = snap_to_grid(charcoal$ybpI, charcoal$accIS,  common_grid, yrInterp),
    ybpI        = charcoal$ybpI
  )
}


# ---- main loop (sequential or parallel) -----------------------------
message(sprintf("Running %d iterations on %d core(s)...", n_iter, n_cores))
t0 <- proc.time()[["elapsed"]]

if (n_cores == 1L) {

  # Sequential path: progress reported every 100 iterations
  results <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    results[[i]] <- run_one_iteration(i)
    if (i %% 100L == 0L) {
      elapsed <- proc.time()[["elapsed"]] - t0
      message(sprintf("  Iteration %4d / %d  (%.1f min elapsed, ~%.1f min remaining)",
                      i, n_iter, elapsed / 60, elapsed / i * (n_iter - i) / 60))
    }
  }

} else {

  # Parallel path: socket cluster (works on Windows and Mac).
  # Iterations are run in 10 batches so progress can be reported
  # between batches, letting users know the run is active.
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)  # clean up even if error

  parallel::clusterEvalQ(cl, library(CharAnalysis))
  parallel::clusterExport(cl,
    varlist = c("chron_depths", "chron_mat", "all_query_depths",
                "cm_top", "cm_bot", "base_data", "params",
                "common_grid", "yrInterp", "snap_to_grid"),
    envir = environment()
  )

  # Split iteration indices into 10 equal batches
  batch_size <- ceiling(n_iter / 10L)
  batches    <- split(seq_len(n_iter),
                      ceiling(seq_len(n_iter) / batch_size))

  results <- vector("list", n_iter)
  for (b in seq_along(batches)) {
    batch_idx        <- batches[[b]]
    batch_results    <- parallel::parLapply(cl, batch_idx, run_one_iteration)
    results[batch_idx] <- batch_results

    elapsed   <- proc.time()[["elapsed"]] - t0
    completed <- max(batch_idx)
    remaining <- elapsed / completed * (n_iter - completed)
    message(sprintf("  Iteration %4d / %d  (%.1f min elapsed, ~%.1f min remaining)",
                    completed, n_iter, elapsed / 60, remaining / 60))
  }
}

message(sprintf("Done. Total time: %.1f min", (proc.time()[["elapsed"]] - t0) / 60))


# ---- assemble output matrices ---------------------------------------
peaks_matrix    <- matrix(NA_real_, nrow = n_steps, ncol = n_iter)
charPeak_matrix <- matrix(NA_real_, nrow = n_steps, ncol = n_iter)
charBkg_matrix  <- matrix(NA_real_, nrow = n_steps, ncol = n_iter)
iter_age_grids  <- vector("list", n_iter)

for (i in seq_len(n_iter)) {
  peaks_matrix[, i]    <- results[[i]]$peaks_final
  charPeak_matrix[, i] <- results[[i]]$c_peak
  charBkg_matrix[, i]  <- results[[i]]$c_bkg
  iter_age_grids[[i]]  <- results[[i]]$ybpI
}


# ---- per-depth chronological uncertainty ribbon ---------------------
# One row per depth in the chronology matrix.
# median_age  : median calibrated age across iterations (cal yr BP)
# ci95_lo/hi  : 2.5th and 97.5th percentile ages across iterations
# ci95_width  : full 95% CI width (yr)  -- used in the uncertainty ribbon panel
message("Computing per-depth chronological uncertainty (95% CI)...")
chron_mat_ages <- as.matrix(chron_mat[, -1L])
chron_ci <- data.frame(
  depth_cm   = chron_depths,
  median_age = apply(chron_mat_ages, 1L, median),
  ci95_lo    = apply(chron_mat_ages, 1L, quantile, probs = 0.025),
  ci95_hi    = apply(chron_mat_ages, 1L, quantile, probs = 0.975)
)
chron_ci$ci95_width <- chron_ci$ci95_hi - chron_ci$ci95_lo


# ---- assemble ensemble object ---------------------------------------
ensemble <- list(
  peaks_matrix    = peaks_matrix,
  charPeak_matrix = charPeak_matrix,
  charBkg_matrix  = charBkg_matrix,
  prob_peak       = rowMeans(peaks_matrix, na.rm = TRUE),
  age_grid        = common_grid,
  iter_age_grids  = iter_age_grids,
  n_iter          = n_iter,
  params          = params,
  chron_ci        = chron_ci      # per-depth CI ribbon (stored so plot works from RDS)
)
class(ensemble) <- c("CharEnsemble", "list")

message(sprintf("ensemble$peaks_matrix : %d steps x %d iterations",
                nrow(ensemble$peaks_matrix), ncol(ensemble$peaks_matrix)))
message(sprintf("ensemble$prob_peak    : range [%.3f, %.3f]",
                min(ensemble$prob_peak, na.rm = TRUE),
                max(ensemble$prob_peak, na.rm = TRUE)))


# ---- run peak matching ----------------------------------------------
source("CharAnalysis_2_0_R/tests/run_ensemble_analysis.R")


# ---- determine whether to save (data + figure PDF) ------------------
do_save <- if (save_results == 1L) {
  TRUE
} else {
  ans <- readline("Save results and figure PDF? (y/n): ")
  tolower(trimws(ans)) == "y"
}
save_pdf <- do_save   # passed through to plot_ensemble_figure.R


# ---- generate figure ------------------------------------------------
source("CharAnalysis_2_0_R/tests/plot_ensemble_figure.R")


# ---- save data files ------------------------------------------------
out_rds <- sprintf("CharAnalysis_2_0_R/tests/%s_ensemble_results.rds",
                   params$site)
out_csv <- sprintf("CharAnalysis_2_0_R/tests/%s_peakAgeUncertainty.csv",
                   out$site)

if (do_save) {
  saveRDS(ensemble, out_rds)
  write.csv(out_df, out_csv, row.names = FALSE)
  message(sprintf("Saved: %s", out_rds))
  message(sprintf("Saved: %s", out_csv))
} else {
  message("Results not saved.")
}
