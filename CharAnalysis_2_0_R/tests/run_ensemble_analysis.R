# run_ensemble_analysis.R
# -----------------------------------------------------------------------
# Reference run -> load ensemble -> nearest-neighbor peak matching.
# Saving is handled by the calling script (char_run_ensemble.R) or
# manually after sourcing standalone.
#
# Can be sourced standalone from a fresh session:
#   setwd("C:/Users/philip.higuera/OneDrive - The University of Montana/1_phiguera/1_working/CharAnalysis")
#   source("CharAnalysis_2_0_R/tests/run_ensemble_analysis.R")
#
# When sourced from char_run_ensemble.R, `out` and `ensemble` are
# already in memory and both load steps are skipped automatically.
# -----------------------------------------------------------------------

library(CharAnalysis)

setwd("C:/Users/philip.higuera/OneDrive - The University of Montana/1_phiguera/1_working/CharAnalysis")

# ---- 1. Reference run --------------------------------------------------
# Skipped if `out` already exists in the workspace (e.g., sourced from
# char_run_ensemble.R, or already run this session).
if (!exists("out")) {
  message("Loading package...")
  devtools::load_all("CharAnalysis_2_0_R")
  params_file <- "CharAnalysis_2_0_R/inst/validation/CH10_charParams.csv"
  message("Running reference CharAnalysis...")
  out <- CharAnalysis(params_file, plots = FALSE)
} else {
  message("Using existing 'out' from workspace.")
}


# ---- 2. Load ensemble --------------------------------------------------
# Skipped if `ensemble` already exists in the workspace.
if (!exists("ensemble")) {
  ensemble_rds <- sprintf("CharAnalysis_2_0_R/tests/%s_ensemble_results.rds",
                          out$site)
  if (!file.exists(ensemble_rds)) {
    stop("Ensemble RDS not found: ", ensemble_rds,
         "\nRun char_run_ensemble.R first to generate it.")
  }
  message("Loading ensemble from RDS...")
  ensemble <- readRDS(ensemble_rds)
} else {
  message("Using existing 'ensemble' from workspace.")
}

message(sprintf("  %d iterations  |  %d time steps  |  age grid %.0f-%.0f cal yr BP",
                ensemble$n_iter, length(ensemble$age_grid),
                min(ensemble$age_grid), max(ensemble$age_grid)))


# ---- 3. Reference peak ages --------------------------------------------
ccp           <- out$post$CharcoalCharPeaks
sel_col       <- ncol(ccp)
x             <- out$charcoal$ybpI
ref_peak_idx  <- which(ccp[, sel_col] > 0)
ref_peak_ages <- x[ref_peak_idx]
n_peaks       <- length(ref_peak_ages)

message(sprintf("\nReference run: %d peaks identified", n_peaks))
message(sprintf("  Age range: %.0f-%.0f cal yr BP",
                min(ref_peak_ages), max(ref_peak_ages)))


# ---- 4. Matching window ------------------------------------------------
yrInterp    <- out$pretreatment$yrInterp
fri_params  <- out$post$FRI_params_zone
valid_zones <- !is.na(fri_params[, 2L]) & fri_params[, 2L] > 0
if (any(valid_zones)) {
  mean_mFRI     <- mean(fri_params[valid_zones, 2L])
  match_halfwin <- max(mean_mFRI / 2, yrInterp)
} else {
  mean_mFRI     <- NA_real_
  match_halfwin <- 500
}
message(sprintf("  CharAnalysis mean mFRI = %.0f yr  |  matching half-window = %.0f yr  |  yrInterp = %.0f yr",
                mean_mFRI, match_halfwin, yrInterp))

raw_intervals <- abs(diff(sort(ref_peak_ages)))
pct_below_win <- mean(raw_intervals < match_halfwin) * 100
if (pct_below_win > 10) {
  message(sprintf(
    "  NOTE: %.0f%% of inter-peak intervals < matching window (%.0f yr).",
    pct_below_win, match_halfwin),
    "\n        Close peaks will compete for ensemble detections -- interpret probabilities with care.")
}


# ---- 5. Nearest-neighbor matching --------------------------------------
age_grid  <- ensemble$age_grid
peaks_mat <- ensemble$peaks_matrix
n_iter    <- ensemble$n_iter

assigned <- vector("list", n_peaks)
for (k in seq_len(n_peaks)) assigned[[k]] <- numeric(0)

for (k in seq_len(n_peaks)) {
  for (i in seq_len(n_iter)) {
    det_idx <- which(peaks_mat[, i] == 1L)
    if (length(det_idx) == 0L) next
    det_ages <- age_grid[det_idx]
    dists    <- abs(det_ages - ref_peak_ages[k])
    best     <- which.min(dists)
    if (dists[best] <= match_halfwin) {
      assigned[[k]] <- c(assigned[[k]], det_ages[best])
    }
  }
}


# ---- 6. Per-peak summaries ---------------------------------------------
peak_summaries <- data.frame(
  ref_age    = ref_peak_ages,
  n_iter     = n_iter,
  prob       = sapply(assigned, function(v) length(v) / n_iter),
  n_detected = sapply(assigned, length),
  mean_age   = sapply(assigned, function(v) if (length(v) > 0L) mean(v)            else NA_real_),
  median_age = sapply(assigned, function(v) if (length(v) > 0L) median(v)          else NA_real_),
  sd_age     = sapply(assigned, function(v) if (length(v) > 1L) sd(v)              else NA_real_),
  ci80_lo    = sapply(assigned, function(v) if (length(v) > 0L) quantile(v, 0.10)  else NA_real_),
  ci80_hi    = sapply(assigned, function(v) if (length(v) > 0L) quantile(v, 0.90)  else NA_real_),
  ci95_lo    = sapply(assigned, function(v) if (length(v) > 0L) quantile(v, 0.025) else NA_real_),
  ci95_hi    = sapply(assigned, function(v) if (length(v) > 0L) quantile(v, 0.975) else NA_real_)
)

any_over_1 <- sum(peak_summaries$prob > 1)
if (any_over_1 > 0L) {
  warning(sprintf("%d peaks have prob > 1.0 -- check matching window", any_over_1))
}

message("\nPeak summaries (first 10 rows):")
print(head(peak_summaries, 10L))
message(sprintf("  prob range: [%.3f, %.3f]",
                min(peak_summaries$prob), max(peak_summaries$prob)))


# ---- 7. Build output data frame (charResults-style column names) -------
# Kept in memory as `out_df`; saving is handled by the calling script
# (char_run_ensemble.R) or manually after sourcing standalone.
out_df <- data.frame(
  "age Top_i (yr BP)"       = peak_summaries$ref_age,
  "n ensemble (iterations)" = peak_summaries$n_iter,
  "detection freq (%)"      = round(peak_summaries$prob * 100, 1),
  "n detected (peaks)"      = peak_summaries$n_detected,
  "mean age (yr BP)"        = round(peak_summaries$mean_age,   0),
  "median age (yr BP)"      = round(peak_summaries$median_age, 0),
  "sd age (yr)"             = round(peak_summaries$sd_age,     1),
  "CI80 lo (yr BP)"         = round(peak_summaries$ci80_lo,    1),
  "CI80 hi (yr BP)"         = round(peak_summaries$ci80_hi,    1),
  "CI95 lo (yr BP)"         = round(peak_summaries$ci95_lo,    1),
  "CI95 hi (yr BP)"         = round(peak_summaries$ci95_hi,    1),
  check.names = FALSE
)

message("Objects in workspace: out, ensemble, peak_summaries, out_df")
