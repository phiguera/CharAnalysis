# SI17_bacon_compare.R
# ---------------------------------------------------------------------
# Compare new Bacon run (rbacon 3.5.2) median ages to the published
# age-depth model stored in SI17_charData.csv.
#
# Run this in the same R session immediately after Bacon("SI17", ...)
# completes, so that the 'info' object is still in memory.
# ---------------------------------------------------------------------

library(rbacon)

setwd("C:/Users/philip.higuera/OneDrive - The University of Montana/1_phiguera/1_working/CharAnalysis")

# --- 1. Load published ages -------------------------------------------
chardata <- read.csv("CharAnalysis_2_0_R/tests/SI17_charData.csv",
                     check.names = FALSE)
depths_char  <- chardata$cmTop
ages_pub     <- chardata$`ageTop (yr BP)`

# --- 2. Check Bacon.Age.d() return format ----------------------------
# This tells us whether we get a full MCMC vector or just summaries.
test <- Bacon.Age.d(50, set = get('info'), BCAD = FALSE)
cat("Bacon.Age.d(50) class:", class(test), "\n")
cat("Length:", length(test), "\n")
cat("First 5 values:", head(test, 5), "\n\n")

# --- 3. Get Bacon median ages at each charData depth -----------------
cat("Extracting Bacon ages at", length(depths_char), "depths...\n")
bacon_ages <- sapply(depths_char, function(d) {
  median(Bacon.Age.d(d, set = get('info'), BCAD = FALSE))
})

# --- 4. Compare -------------------------------------------------------
diffs <- bacon_ages - ages_pub

cat("\n--- Comparison: Bacon median vs. published ---\n")
cat("Max absolute difference:", round(max(abs(diffs)), 1), "yr\n")
cat("Mean absolute difference:", round(mean(abs(diffs)), 1), "yr\n")
cat("RMSE:", round(sqrt(mean(diffs^2)), 1), "yr\n")

# Flag any depth where difference exceeds 100 yr
large_diffs <- which(abs(diffs) > 100)
if (length(large_diffs) > 0) {
  cat("\nDepths with |diff| > 100 yr:\n")
  print(data.frame(
    cmTop      = depths_char[large_diffs],
    published  = ages_pub[large_diffs],
    bacon_med  = round(bacon_ages[large_diffs], 1),
    diff_yr    = round(diffs[large_diffs], 1)
  ))
} else {
  cat("No depths with |diff| > 100 yr.\n")
}

# --- 5. Plot ----------------------------------------------------------
plot(ages_pub, bacon_ages,
     xlab = "Published age (cal yr BP)",
     ylab = "Bacon median age (cal yr BP)",
     main = "SI17: Published vs. new Bacon run",
     pch  = 16, cex = 0.6, col = "steelblue")
abline(0, 1, col = "red", lwd = 1.5)
legend("topleft",
       legend = c("1:1 line"),
       col    = "red", lty = 1, lwd = 1.5, bty = "n")
