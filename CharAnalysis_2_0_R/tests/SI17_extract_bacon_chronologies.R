# SI17_extract_bacon_chronologies.R
# ---------------------------------------------------------------------
# Extract 1000 chronologies from the SI17 Bacon MCMC posterior and
# save in the format expected by char_run_ensemble.R:
#   col 1  : Sample_cm  (depth, cm)
#   cols 2–1001: CalAge_1 … CalAge_1000  (cal yr BP)
#
# Run in the same R session immediately after Bacon("SI17", ...)
# so that the 'info' object is still in memory.
# ---------------------------------------------------------------------

library(rbacon)

setwd("C:/Users/philip.higuera/OneDrive - The University of Montana/1_phiguera/1_working/CharAnalysis")

# --- 1. Depths to extract --------------------------------------------
# Use cmTop from SI17_charData.csv (the 0.5-cm intervals, 0–431.5 cm)
chardata <- read.csv("CharAnalysis_2_0_R/tests/SI17_charData.csv",
                     check.names = FALSE)
depths <- chardata$cmTop
cat("Extracting ages at", length(depths), "depths",
    "(", min(depths), "–", max(depths), "cm)\n\n")

# --- 2. Sample 1000 MCMC indices -------------------------------------
# All Bacon.Age.d() calls return 8000-length vectors ordered by
# MCMC iteration, so using the same indices across depths gives
# properly correlated chronologies.
set.seed(42)
n_mcmc  <- 8000
n_iter  <- 1000
idx     <- sample(n_mcmc, n_iter, replace = FALSE)

# --- 3. Extract ages at every depth ----------------------------------
cat("Extracting ages (this may take a minute)...\n")
t0 <- proc.time()

age_matrix <- matrix(NA_real_, nrow = length(depths), ncol = n_iter)

for (i in seq_along(depths)) {
  all_ages       <- Bacon.Age.d(depths[i], set = get('info'), BCAD = FALSE)
  age_matrix[i,] <- all_ages[idx]
  if (i %% 100 == 0) {
    elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
    cat("  depth", i, "of", length(depths), "—", elapsed, "min elapsed\n")
  }
}

elapsed_total <- round((proc.time() - t0)["elapsed"] / 60, 1)
cat("Done. Total time:", elapsed_total, "min\n\n")

# --- 4. Assemble and write CSV ---------------------------------------
chron_df           <- as.data.frame(age_matrix)
names(chron_df)    <- paste0("CalAge_", seq_len(n_iter))
chron_df           <- cbind(Sample_cm = depths, chron_df)

out_file <- "CharAnalysis_2_0_R/tests/SI17_MCAgeDepth_1000_chronologies.csv"
write.csv(chron_df, out_file, row.names = FALSE)
cat("Saved:", out_file, "\n")
cat("Dimensions:", nrow(chron_df), "depths x", ncol(chron_df) - 1, "chronologies\n")

# --- 5. Quick sanity check -------------------------------------------
# Compare iteration medians to published ages
pub_ages    <- chardata$`ageTop (yr BP)`
bacon_meds  <- apply(age_matrix, 1, median)
diffs       <- bacon_meds - pub_ages
cat("\nSanity check (median of 1000 vs. published):\n")
cat("  Max |diff|:", round(max(abs(diffs)), 1), "yr\n")
cat("  Mean |diff|:", round(mean(abs(diffs)), 1), "yr\n")
