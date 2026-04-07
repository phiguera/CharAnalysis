# =============================================================================
# test_phase3.R
# Phase 3 validation: peak identification and minimum-count screening
#
# Datasets validated:
#   - CO  (Code Lake, local GMM threshold — primary test case)
#
# Columns validated here:
#   peaks 1, peaks 2, peaks 3, peaks Final
#   (peaks Insig. is assembled in CharPostProcess — validated in Phase 4)
#
# Run from the CharAnalysis_R/ directory:
#   source("tests/test_phase3.R")
#
# Or from any directory:
#   CHAR_R_ROOT <- "/path/to/CharAnalysis_R"
#   source(file.path(CHAR_R_ROOT, "tests/test_phase3.R"))
#
# ---------------------------------------------------------------------------
# VALIDATION APPROACH
# ---------------------------------------------------------------------------
# Phase 3 outputs (peak positions, counts) derive from threshold values that
# differ between R and MATLAB because EM implementations accumulate floating-
# point divergence across iterations (GaussianMixture
# MDL/first-last init).  This means peak positions cannot be expected to match
# exactly.  The test therefore has two sections:
#
# SECTION A — Structural invariants (strict PASS/FAIL)
#   Properties that must hold regardless of which threshold implementation is
#   used.  These test the correctness of the peak-flagging and
#   consecutive-removal logic itself:
#     1. No two adjacent samples in the same charPeaks column are both 1.
#     2. All minCountP values are in [0, 1] (or NaN for non-peaks).
#     3. peaksTotal equals colSums(charPeaks) for each column.
#
# SECTION B — Peak count comparison (informational)
#   R peak counts vs MATLAB peak counts per threshold column.  Differences
#   are expected due to GMM threshold differences.  The test reports the
#   counts and the agreement in peak positions (% of MATLAB peaks that also
#   appear in R output, and vice versa) so the practical impact can be
#   assessed.
#
# SECTION C — Sample-row table (visual inspection)
#   Side-by-side R vs MATLAB values at representative rows for charBkg,
#   charPeak, thresh FinalPos, peaks Final.
# =============================================================================

# ---- Paths ------------------------------------------------------------------
if (!exists("CHAR_R_ROOT")) {
  CHAR_R_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."),
                                mustWork = FALSE)
  if (!dir.exists(CHAR_R_ROOT)) CHAR_R_ROOT <- getwd()
}

source(file.path(CHAR_R_ROOT, "R", "utils_lowess.R"))
source(file.path(CHAR_R_ROOT, "R", "utils_gaussian_mixture.R"))
source(file.path(CHAR_R_ROOT, "R", "char_parameters.R"))
source(file.path(CHAR_R_ROOT, "R", "char_validate_params.R"))
source(file.path(CHAR_R_ROOT, "R", "char_pretreatment.R"))
source(file.path(CHAR_R_ROOT, "R", "char_smooth.R"))
source(file.path(CHAR_R_ROOT, "R", "char_thresh_global.R"))
source(file.path(CHAR_R_ROOT, "R", "char_thresh_local.R"))
source(file.path(CHAR_R_ROOT, "R", "char_peak_id.R"))
source(file.path(CHAR_R_ROOT, "R", "CharAnalysis.R"))

MATLAB_DIR <- file.path(CHAR_R_ROOT, "..", "CharAnalysis_2_0_MATLAB")

# =============================================================================
# Run CharAnalysis() once; reuse result throughout
# =============================================================================

cat("\nRunning CharAnalysis() for CO (Code Lake)...\n")
r_full <- suppressMessages(
  CharAnalysis(file.path(MATLAB_DIR, "CO_charParams.csv"))
)
cat("Done.\n\n")

# Load MATLAB reference CSV
m_raw   <- read.csv(file.path(MATLAB_DIR, "CO_charResults.csv"),
                    header = TRUE, check.names = FALSE)
m_names <- trimws(names(m_raw))
names(m_raw) <- m_names
mcol <- function(pattern) {
  idx <- grep(pattern, m_names, fixed = TRUE)
  if (!length(idx)) stop("MATLAB column not found: '", pattern, "'")
  as.numeric(m_raw[[idx[1L]]])
}

N <- length(r_full$charcoal$ybpI)
T_thresh <- ncol(r_full$charcoal$charPeaks)   # number of threshold columns

# =============================================================================
# SECTION A — Structural invariants  (strict PASS / FAIL)
# =============================================================================

cat(strrep("=", 72), "\n")
cat(" Phase 3 validation: Code Lake (CO) — local GMM threshold\n")
cat(strrep("=", 72), "\n\n")

cat("--- A. Structural invariants (strict PASS/FAIL) ---\n\n")

pass_all <- TRUE

## A1. No two consecutive flagged samples in any column
no_consec_ok <- TRUE
for (j in seq_len(T_thresh)) {
  pk <- r_full$charcoal$charPeaks[, j]
  if (any(pk[-N] > 0 & pk[-1L] > 0)) {
    cat(sprintf("  FAIL: consecutive peaks found in charPeaks column %d\n", j))
    no_consec_ok <- FALSE
    pass_all     <- FALSE
  }
}
if (no_consec_ok) {
  cat(sprintf("  PASS: no consecutive peaks in any of %d threshold columns\n",
              T_thresh))
}

## A2. minCountP values in [0, 1] or NA
mcp <- r_full$char_thresh$minCountP
mcp_vals <- mcp[!is.na(mcp)]
if (length(mcp_vals) == 0L) {
  cat("  INFO: minCountP is all NA (no peaks with >= 2 neighbours).\n")
} else if (all(mcp_vals >= 0 & mcp_vals <= 1)) {
  cat(sprintf("  PASS: all %d non-NA minCountP values in [0, 1]\n",
              length(mcp_vals)))
} else {
  out_of_range <- sum(mcp_vals < 0 | mcp_vals > 1)
  cat(sprintf("  FAIL: %d minCountP values outside [0, 1]\n", out_of_range))
  pass_all <- FALSE
}

## A3. peaksTotal == colSums(charPeaks)
col_sums <- colSums(r_full$charcoal$charPeaks)
totals_match <- all(r_full$charcoal$peaksTotal == col_sums)
if (totals_match) {
  cat(sprintf("  PASS: peaksTotal matches colSums(charPeaks) for all %d columns\n",
              T_thresh))
} else {
  cat("  FAIL: peaksTotal does not match colSums(charPeaks)\n")
  cat("    peaksTotal : ", paste(r_full$charcoal$peaksTotal, collapse = ", "), "\n")
  cat("    colSums    : ", paste(col_sums, collapse = ", "), "\n")
  pass_all <- FALSE
}

if (!pass_all) {
  stop("\n  One or more structural invariants FAILED.")
} else {
  cat("\n  All structural invariants: PASS\n\n")
}

# =============================================================================
# SECTION B — Peak count and position agreement  (informational)
# =============================================================================

cat("--- B. Peak count comparison vs MATLAB (informational) ---\n\n")

# MATLAB peak columns (binary 0/1): peaks 1, peaks 2, peaks 3, peaks Final
m_peaks_cols <- list(
  mcol("peaks 1"),
  mcol("peaks 2"),
  mcol("peaks 3"),
  mcol("peaks Final")
)
col_labels <- c("thresh 1 (95%)", "thresh 2 (99%)",
                "thresh 3 (99.9%)", "thresh Final (99%)")

cat(sprintf("%-20s  %8s  %8s  %8s  %9s  %9s\n",
            "Column", "R peaks", "M peaks", "Both", "R-only%", "M-only%"))
cat(strrep("-", 72), "\n")

for (j in seq_len(T_thresh)) {
  r_pk <- r_full$charcoal$charPeaks[, j] > 0
  m_pk <- m_peaks_cols[[j]] > 0

  n_r    <- sum(r_pk)
  n_m    <- sum(m_pk)
  n_both <- sum(r_pk & m_pk)
  r_only_pct <- if (n_r > 0) 100 * (n_r - n_both) / n_r else 0
  m_only_pct <- if (n_m > 0) 100 * (n_m - n_both) / n_m else 0

  cat(sprintf("%-20s  %8d  %8d  %8d  %8.1f%%  %8.1f%%\n",
              col_labels[j], n_r, n_m, n_both, r_only_pct, m_only_pct))
}

cat("\n  Columns: 'R peaks' / 'M peaks' = total peaks identified by R / MATLAB.\n")
cat("  'Both'   = positions where both R and MATLAB identified a peak.\n")
cat("  'R-only%'= % of R peaks not in MATLAB; 'M-only%'= % of MATLAB peaks\n")
cat("  not in R.  Differences arise because R and MATLAB's EM implementations\n")
cat("  accumulate floating-point divergence across iterations, settling at\n")
cat("  different local solutions despite identical initialisation and epsilon.\n")
cat("  The 'thresh Final' column drives FRIs and Weibull statistics.\n")
cat("  Discrepancy (~20%% in peak count) is documented and accepted.\n\n")

# =============================================================================
# SECTION C — Sample-row table for visual inspection
# =============================================================================

cat("--- C. Sample rows — visual inspection ---\n\n")

inspect_rows <- c(1L, 50L, 100L, 200L, 300L, 500L)
inspect_rows <- inspect_rows[inspect_rows <= N]

r_bkg        <- r_full$charcoal$accIS
r_peak_ser   <- r_full$charcoal$peak
r_thresh_fin <- r_full$char_thresh$pos[, T_thresh]
r_peaks_fin  <- r_full$charcoal$charPeaks[, T_thresh]

m_bkg        <- mcol("charBkg")
m_peak_ser   <- mcol("char Peak")
m_thresh_fin <- mcol("thresh FinalPos")
m_peaks_fin  <- mcol("peaks Final")

cat(sprintf("%-5s  %-10s  %-10s  %-10s  %-10s  %-8s  %-8s\n",
            "Row", "R_charBkg", "M_charBkg",
            "R_thrFin", "M_thrFin",
            "R_pkFin", "M_pkFin"))
cat(strrep("-", 68), "\n")

for (i in inspect_rows) {
  cat(sprintf("%-5d  %-10.5g  %-10.5g  %-10.5g  %-10.5g  %-8d  %-8d\n",
              i,
              r_bkg[i],        m_bkg[i],
              r_thresh_fin[i], m_thresh_fin[i],
              as.integer(r_peaks_fin[i]),
              as.integer(m_peaks_fin[i])))
}

# =============================================================================
# DONE
# =============================================================================

cat("\n")
cat(strrep("*", 72), "\n")
cat(" Phase 3 validation complete for Code Lake (CO).\n")
cat(strrep("*", 72), "\n\n")
