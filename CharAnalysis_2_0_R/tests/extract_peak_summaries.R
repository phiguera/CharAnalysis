# extract_peak_summaries.R
# -----------------------------------------------------------------------
# Nearest-neighbor peak matching: for each reference peak identified by
# the single-chronology CharAnalysis run, find all ensemble iterations
# that detected a peak within half the median FRI and record the detected
# age.  Outputs a CSV of per-peak summaries for figure production.
#
# Requires in workspace:
#   ensemble  -- CharEnsemble from char_run_ensemble.R / RDS
#   out       -- CharAnalysis output from reference run
#
# Output:
#   peak_summaries  data.frame in workspace
#   tests/CH10_peak_summaries.csv  on disk
# -----------------------------------------------------------------------

# ---- reference peak ages -----------------------------------------------
post         <- out$post
charcoal     <- out$charcoal
peak_analysis <- out$peak_analysis

# peaksFinal column from the reference run
ccp     <- post$CharcoalCharPeaks
sel_col <- ncol(ccp)
x       <- charcoal$ybpI
ref_peak_idx  <- which(ccp[, sel_col] > 0)
ref_peak_ages <- x[ref_peak_idx]

cat(sprintf("Reference run: %d peaks identified\n", length(ref_peak_ages)))
cat(sprintf("  Ages (cal yr BP): %s\n",
            paste(round(ref_peak_ages), collapse = ", ")))


# ---- matching window: half the median FRI ------------------------------
# Use the reference run's fire return intervals to set the window.
# If fewer than 2 peaks, fall back to 500 yr.
if (length(ref_peak_ages) >= 2L) {
  fri_median  <- median(abs(diff(sort(ref_peak_ages))))
  match_halfwin <- fri_median / 2
} else {
  match_halfwin <- 500
}
cat(sprintf("Median FRI = %.0f yr  |  matching half-window = %.0f yr\n",
            fri_median, match_halfwin))


# ---- nearest-neighbor matching -----------------------------------------
age_grid    <- ensemble$age_grid
peaks_mat   <- ensemble$peaks_matrix   # n_steps x n_iter
n_iter      <- ensemble$n_iter

# For each iteration, find detected peak ages and assign to nearest ref peak
assigned <- vector("list", length(ref_peak_ages))
for (k in seq_along(ref_peak_ages)) assigned[[k]] <- numeric(0)

for (i in seq_len(n_iter)) {
  det_idx  <- which(peaks_mat[, i] == 1L)
  if (length(det_idx) == 0L) next
  det_ages <- age_grid[det_idx]
  for (da in det_ages) {
    dists <- abs(da - ref_peak_ages)
    nearest <- which.min(dists)
    if (dists[nearest] <= match_halfwin) {
      assigned[[nearest]] <- c(assigned[[nearest]], da)
    }
  }
}


# ---- per-peak summaries ------------------------------------------------
peak_summaries <- data.frame(
  ref_age    = ref_peak_ages,
  prob       = sapply(assigned, function(v) length(v) / n_iter),
  n_detected = sapply(assigned, length),
  mean_age   = sapply(assigned, function(v) if (length(v) > 0) mean(v)   else NA_real_),
  median_age = sapply(assigned, function(v) if (length(v) > 0) median(v) else NA_real_),
  sd_age     = sapply(assigned, function(v) if (length(v) > 1) sd(v)     else NA_real_),
  ci80_lo    = sapply(assigned, function(v) if (length(v) > 0) quantile(v, 0.10) else NA_real_),
  ci80_hi    = sapply(assigned, function(v) if (length(v) > 0) quantile(v, 0.90) else NA_real_),
  ci95_lo    = sapply(assigned, function(v) if (length(v) > 0) quantile(v, 0.025) else NA_real_),
  ci95_hi    = sapply(assigned, function(v) if (length(v) > 0) quantile(v, 0.975) else NA_real_)
)

print(peak_summaries)


# ---- save --------------------------------------------------------------
out_csv <- "CharAnalysis_2_0_R/tests/CH10_peak_summaries.csv"
write.csv(peak_summaries, out_csv, row.names = FALSE)
cat(sprintf("Saved: %s\n", out_csv))
