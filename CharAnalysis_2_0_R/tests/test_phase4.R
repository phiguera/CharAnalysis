# =============================================================================
# test_phase4.R
# Phase 4 validation: post-processing — fire-return intervals, fire frequency,
# Weibull statistics, and char_results output matrix
#
# Datasets validated:
#   - CO  (Code Lake, local GMM threshold — primary test case)
#
# Quantities validated:
#   peaks Insig.    (binary 0/1 at insignificant-peak positions)
#   peak Mag        (integrated C_peak above threshold, pieces cm-2 peak-1)
#   smPeak Frequ    (Lowess-smoothed fire frequency, peaks ka-1)
#   smFRIs          (Lowess-smoothed mean FRI, yr fire-1)
#   nFRI, mFRI      (per-zone FRI count and mean)
#   WBLb, WBLc      (Weibull scale and shape parameters)
#
# Run from the CharAnalysis_R/ directory:
#   source("tests/test_phase4.R")
#
# Or from any directory:
#   CHAR_R_ROOT <- "/path/to/CharAnalysis_R"
#   source(file.path(CHAR_R_ROOT, "tests/test_phase4.R"))
#
# ---------------------------------------------------------------------------
# VALIDATION APPROACH
# ---------------------------------------------------------------------------
# Phase 4 outputs depend heavily on Phase 3 peak positions, which differ
# between R and MATLAB due to irreducible floating-point divergence in the
# EM algorithm (documented in test_phase3.R).  Two test tiers apply:
#
# SECTION A — Structural invariants (strict PASS/FAIL)
#   Properties that must hold regardless of peak-position differences:
#     1. char_results has exactly 33 columns and N rows.
#     2. peakInsig is binary (0 or 1 only).
#     3. peakMagnitude >= 0 everywhere.
#     4. smoothedFireFrequ >= 0 everywhere.
#     5. smFRIs values (where non-NA) are positive.
#     6. Weibull parameters (WBLb, WBLc) are positive where fitted.
#
# SECTION B — Smoothed fire frequency comparison (informational)
#   R vs MATLAB smPeak Frequ.  This series is only indirectly affected by
#   peak-position differences, because the sliding window averages over a
#   wide (~1000 yr) band.  Moderate agreement is expected.
#
# SECTION C — Smoothed FRI comparison (informational)
#   R vs MATLAB smFRIs on the overlapping ybpI range.  The smoothed FRI
#   averages over many peaks so agreement is expected within ~20-30%.
#
# SECTION D — Per-zone Weibull / FRI statistics (informational)
#   R vs MATLAB nFRI, mFRI, WBLb, WBLc.  mFRI is computed from different
#   peak sets so exact agreement is not expected; the test reports both
#   values side-by-side.
#
# SECTION E — Sample-row table (visual inspection)
#   Side-by-side R vs MATLAB values at representative rows.
# =============================================================================

# ---- Paths ------------------------------------------------------------------
if (!exists("CHAR_R_ROOT")) {
  CHAR_R_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."),
                                mustWork = FALSE)
  if (!dir.exists(CHAR_R_ROOT)) CHAR_R_ROOT <- getwd()
}

source(file.path(CHAR_R_ROOT, "R", "charLowess.R"))
source(file.path(CHAR_R_ROOT, "R", "GaussianMixture.R"))
source(file.path(CHAR_R_ROOT, "R", "CharParameters.R"))
source(file.path(CHAR_R_ROOT, "R", "CharValidateParams.R"))
source(file.path(CHAR_R_ROOT, "R", "CharPretreatment.R"))
source(file.path(CHAR_R_ROOT, "R", "CharSmooth.R"))
source(file.path(CHAR_R_ROOT, "R", "CharThreshGlobal.R"))
source(file.path(CHAR_R_ROOT, "R", "CharThreshLocal.R"))
source(file.path(CHAR_R_ROOT, "R", "CharPeakID.R"))
source(file.path(CHAR_R_ROOT, "R", "CharPostProcess.R"))
source(file.path(CHAR_R_ROOT, "R", "CharWriteResults.R"))
source(file.path(CHAR_R_ROOT, "R", "CharAnalysis.R"))

MATLAB_DIR <- file.path(CHAR_R_ROOT, "..", "CharAnalysis_2_0_MATLAB")

# =============================================================================
# Run CharAnalysis() once; reuse throughout
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

N        <- length(r_full$charcoal$ybpI)
T_thresh <- ncol(r_full$charcoal$charPeaks)

# =============================================================================
# SECTION A — Structural invariants  (strict PASS / FAIL)
# =============================================================================

cat(strrep("=", 72), "\n")
cat(" Phase 4 validation: Code Lake (CO) — post-processing\n")
cat(strrep("=", 72), "\n\n")

cat("--- A. Structural invariants (strict PASS/FAIL) ---\n\n")

pass_all <- TRUE

## A1. char_results dimensions
cr <- r_full$char_results
if (is.null(cr)) {
  cat("  FAIL: char_results is NULL — char_post_process() may not have run\n")
  pass_all <- FALSE
} else if (!is.matrix(cr)) {
  cat(sprintf("  FAIL: char_results is %s, expected matrix\n", class(cr)))
  pass_all <- FALSE
} else if (nrow(cr) != N || ncol(cr) != 33L) {
  cat(sprintf("  FAIL: char_results is [%d x %d], expected [%d x 33]\n",
              nrow(cr), ncol(cr), N))
  pass_all <- FALSE
} else {
  cat(sprintf("  PASS: char_results is [%d x 33]\n", N))
}

## A2. peakInsig is binary (0 or 1)
pi_vals <- r_full$charcoal$peakInsig
if (is.null(pi_vals)) {
  cat("  FAIL: charcoal$peakInsig is NULL\n")
  pass_all <- FALSE
} else if (!all(pi_vals %in% c(0L, 1L))) {
  cat(sprintf("  FAIL: peakInsig has non-binary values: %s\n",
              paste(unique(pi_vals[!pi_vals %in% c(0L, 1L)]), collapse = ", ")))
  pass_all <- FALSE
} else {
  cat(sprintf("  PASS: peakInsig is binary (sum = %d)\n", sum(pi_vals)))
}

## A3. peakMagnitude >= 0
pm_vals <- r_full$charcoal$peakMagnitude
if (is.null(pm_vals)) {
  cat("  FAIL: charcoal$peakMagnitude is NULL\n")
  pass_all <- FALSE
} else if (any(!is.na(pm_vals) & pm_vals < 0)) {
  cat(sprintf("  FAIL: %d negative peakMagnitude values\n",
              sum(!is.na(pm_vals) & pm_vals < 0)))
  pass_all <- FALSE
} else {
  n_nonzero <- sum(pm_vals > 0, na.rm = TRUE)
  cat(sprintf("  PASS: peakMagnitude >= 0 everywhere (%d non-zero entries)\n",
              n_nonzero))
}

## A4. smoothedFireFrequ >= 0
sff <- r_full$charcoal$smoothedFireFrequ
if (is.null(sff)) {
  cat("  FAIL: charcoal$smoothedFireFrequ is NULL\n")
  pass_all <- FALSE
} else if (any(!is.na(sff) & sff < 0)) {
  cat(sprintf("  FAIL: %d negative smoothedFireFrequ values\n",
              sum(!is.na(sff) & sff < 0)))
  pass_all <- FALSE
} else {
  cat(sprintf("  PASS: smoothedFireFrequ >= 0 everywhere (range: %.3g – %.3g)\n",
              min(sff, na.rm = TRUE), max(sff, na.rm = TRUE)))
}

## A5. smFRIs (col 23) values positive where non-NA
smfris_r <- cr[, 23L]
smfris_ok <- smfris_r[!is.na(smfris_r)]
if (length(smfris_ok) == 0L) {
  cat("  INFO: smFRIs column 23 is all NA\n")
} else if (any(smfris_ok <= 0)) {
  cat(sprintf("  FAIL: %d non-positive smFRIs values\n",
              sum(smfris_ok <= 0)))
  pass_all <- FALSE
} else {
  cat(sprintf("  PASS: smFRIs positive where non-NA (%d values, range %.1f – %.1f yr)\n",
              length(smfris_ok), min(smfris_ok), max(smfris_ok)))
}

## A6. Weibull parameters positive where fitted
wbl_b <- cr[1L, 28L]   # WBLb, zone 1
wbl_c <- cr[1L, 31L]   # WBLc, zone 1
if (!is.na(wbl_b) && wbl_b <= 0) {
  cat(sprintf("  FAIL: WBLb zone 1 = %.4g (must be positive)\n", wbl_b))
  pass_all <- FALSE
} else if (!is.na(wbl_c) && wbl_c <= 0) {
  cat(sprintf("  FAIL: WBLc zone 1 = %.4g (must be positive)\n", wbl_c))
  pass_all <- FALSE
} else {
  wbl_b_str <- if (is.na(wbl_b)) "not fitted" else sprintf("%.4g", wbl_b)
  wbl_c_str <- if (is.na(wbl_c)) "not fitted" else sprintf("%.4g", wbl_c)
  cat(sprintf("  PASS: Weibull parameters non-negative (WBLb=%s, WBLc=%s for zone 1)\n",
              wbl_b_str, wbl_c_str))
}

if (!pass_all) {
  stop("\n  One or more structural invariants FAILED.")
} else {
  cat("\n  All structural invariants: PASS\n\n")
}

# =============================================================================
# SECTION B — Smoothed fire frequency comparison (informational)
# =============================================================================

cat("--- B. Smoothed fire frequency vs MATLAB (informational) ---\n\n")

m_sff  <- mcol("smPeak Frequ")
r_sff  <- r_full$charcoal$smoothedFireFrequ

both_ok <- !is.na(r_sff) & !is.na(m_sff)
if (sum(both_ok) > 0L) {
  abs_diff <- abs(r_sff[both_ok] - m_sff[both_ok])
  cat(sprintf("  Compared %d samples with non-NA values in both R and MATLAB\n",
              sum(both_ok)))
  cat(sprintf("  Mean |diff|: %.4f peaks ka-1\n", mean(abs_diff)))
  cat(sprintf("  Max  |diff|: %.4f peaks ka-1  (at sample %d)\n",
              max(abs_diff),
              which(both_ok)[which.max(abs_diff)]))
  cat(sprintf("  MATLAB mean: %.4f;  R mean: %.4f  peaks ka-1\n",
              mean(m_sff[both_ok]), mean(r_sff[both_ok])))
} else {
  cat("  No overlapping non-NA values to compare.\n")
}
cat("\n  Note: differences arise because R and MATLAB peak sets differ\n")
cat("  (~20%% discrepancy documented in Phase 3).\n\n")

# =============================================================================
# SECTION C — Smoothed FRI comparison (informational)
# =============================================================================

cat("--- C. Smoothed FRI (smFRIs) vs MATLAB (informational) ---\n\n")

m_smfri <- mcol("smFRIs")
r_smfri <- cr[, 23L]

both_fri <- !is.na(r_smfri) & !is.na(m_smfri)
if (sum(both_fri) > 0L) {
  abs_diff_fri <- abs(r_smfri[both_fri] - m_smfri[both_fri])
  cat(sprintf("  Compared %d samples with non-NA smFRIs in both R and MATLAB\n",
              sum(both_fri)))
  cat(sprintf("  Mean |diff|: %.2f yr fire-1\n", mean(abs_diff_fri)))
  cat(sprintf("  Max  |diff|: %.2f yr fire-1  (at sample %d)\n",
              max(abs_diff_fri),
              which(both_fri)[which.max(abs_diff_fri)]))
  cat(sprintf("  MATLAB mean: %.2f;  R mean: %.2f  yr fire-1\n",
              mean(m_smfri[both_fri]), mean(r_smfri[both_fri])))
} else {
  cat("  No overlapping non-NA smFRIs to compare.\n")
}
cat("\n")

# =============================================================================
# SECTION D — Per-zone FRI and Weibull statistics  (informational)
# =============================================================================

cat("--- D. Per-zone FRI / Weibull statistics vs MATLAB (informational) ---\n\n")

# MATLAB stores per-zone values in rows 1 and 2 of cols 24-33
# R stores them the same way in char_results
m_nfri  <- mcol("nFRIs")
m_mfri  <- mcol("mFRI")
m_wblb  <- mcol("WBLb")
m_wblc  <- mcol("WBLc")

# Find non-NA zones in MATLAB (non-empty rows)
m_zone_rows <- which(!is.na(m_nfri))

n_zones_r  <- nrow(r_full$post$FRI_params_zone)
n_zones_m  <- length(m_zone_rows)

cat(sprintf("  Zones fitted: R = %d, MATLAB = %d\n\n", n_zones_r, n_zones_m))

cat(sprintf("  %-8s  %-6s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n",
            "Zone", "Src", "nFRI", "mFRI", "WBLb", "WBLc",
            "mFRI_uCI", "mFRI_lCI"))
cat(strrep("-", 72), "\n")

for (z in seq_len(max(n_zones_r, n_zones_m))) {

  # R values (from char_results matrix, row z)
  if (z <= n_zones_r) {
    r_row  <- cr[z, 24:33]
    r_nfri <- r_row[1L];  r_mfri <- r_row[2L]
    r_muci <- r_row[3L];  r_mlci <- r_row[4L]
    r_wblb <- r_row[5L];  r_wblc <- r_row[8L]
    if (is.na(r_nfri)) {
      cat(sprintf("  %-8d  %-6s  (not fitted)\n", z, "R"))
    } else {
      r_str <- sprintf("%-8d  %-6s  %-8.1f  %-8.1f  %-8.4g  %-8.4g  %-8.1f  %-8.1f",
                       z, "R", r_nfri, r_mfri, r_wblb, r_wblc, r_muci, r_mlci)
      cat("  ", r_str, "\n")
    }
  }

  # MATLAB values
  if (z <= n_zones_m) {
    mz <- m_zone_rows[z]
    m_str <- sprintf("%-8d  %-6s  %-8.1f  %-8.1f  %-8.4g  %-8.4g  %-8.1f  %-8.1f",
                     z, "MATLAB",
                     m_nfri[mz], m_mfri[mz], m_wblb[mz], m_wblc[mz],
                     mcol("mFRI_uCI")[mz], mcol("mFRI_lCI")[mz])
    cat("  ", m_str, "\n")
  }
  if (z < max(n_zones_r, n_zones_m)) cat("\n")
}
cat("\n")
cat("  Note: nFRI and mFRI depend on peak positions, which differ ~20%% between\n")
cat("  R and MATLAB (documented EM floating-point divergence).  Weibull CIs are\n")
cat("  stochastic (bootstrapped), so small run-to-run variation is expected.\n\n")

# =============================================================================
# SECTION E — Sample-row table for visual inspection
# =============================================================================

cat("--- E. Sample rows — visual inspection ---\n\n")

inspect_rows <- c(1L, 10L, 25L, 50L, 100L, 200L, 300L)
inspect_rows <- inspect_rows[inspect_rows <= N]

m_insig  <- mcol("peaks Insig.")
m_pmag   <- mcol("peak Mag")
m_sff_v  <- mcol("smPeak Frequ")
m_smfri_v <- mcol("smFRIs")

cat(sprintf("%-5s  %-10s  %-10s  %-10s  %-10s  %-8s  %-8s\n",
            "Row", "R_smFreq", "M_smFreq",
            "R_smFRI", "M_smFRI",
            "R_Insig", "M_Insig"))
cat(strrep("-", 70), "\n")

for (i in inspect_rows) {
  cat(sprintf("%-5d  %-10.4g  %-10.4g  %-10.4g  %-10.4g  %-8d  %-8d\n",
              i,
              r_full$charcoal$smoothedFireFrequ[i],
              m_sff_v[i],
              ifelse(is.na(cr[i, 23L]), -999, cr[i, 23L]),
              ifelse(is.na(m_smfri_v[i]), -999, m_smfri_v[i]),
              as.integer(r_full$charcoal$peakInsig[i]),
              as.integer(m_insig[i])))
}

# =============================================================================
# SECTION F — CharWriteResults() round-trip test
# Write to a temp directory, read back, compare column headers and a
# selection of numeric values against the MATLAB reference.
# =============================================================================

cat("--- F. CharWriteResults() round-trip test ---\n\n")

tmp_dir  <- file.path(tempdir(), "charanalysis_test")
r_path   <- suppressMessages(
  CharWriteResults(r_full$char_results,
                   site    = "CO",
                   out_dir = tmp_dir)
)

if (!file.exists(r_path)) {
  cat("  FAIL: output file was not created at:", r_path, "\n")
} else {
  cat(sprintf("  File written: %s\n", r_path))

  r_csv <- read.csv(r_path, header = TRUE, check.names = FALSE)
  r_cnames <- trimws(names(r_csv))

  ## F1. Column count
  if (ncol(r_csv) != 33L) {
    cat(sprintf("  FAIL: CSV has %d columns, expected 33\n", ncol(r_csv)))
  } else {
    cat("  PASS: CSV has 33 columns\n")
  }

  ## F2. Header match vs MATLAB reference
  m_cnames <- trimws(names(m_raw))
  mismatched <- which(r_cnames != m_cnames)
  if (length(mismatched) > 0L) {
    cat(sprintf("  FAIL: %d column header(s) differ from MATLAB:\n",
                length(mismatched)))
    for (i in mismatched)
      cat(sprintf("    col %d: R='%s'  MATLAB='%s'\n",
                  i, r_cnames[i], m_cnames[i]))
  } else {
    cat("  PASS: All 33 column headers match MATLAB exactly\n")
  }

  ## F3. Row count
  if (nrow(r_csv) != N) {
    cat(sprintf("  FAIL: CSV has %d rows, expected %d\n", nrow(r_csv), N))
  } else {
    cat(sprintf("  PASS: CSV has %d rows\n", N))
  }

  ## F4. Spot-check numeric values at known non-NA rows
  #  Compare age Top_i (col 2), charBkg (col 7), thresh FinalPos (col 12)
  spot_rows <- c(1L, 10L, 50L, 100L, 200L)
  spot_cols <- c(2L, 7L, 12L, 22L)  # age, charBkg, thrFinalPos, smFireFreq
  spot_cnames <- c("age Top_i (yr BP)", "charBkg (# cm-2 yr-1)",
                   "thresh FinalPos (# cm-2 yr-1)", "smPeak Frequ (peaks 1ka-1)")
  tol <- 1e-4   # relative tolerance for round-trip

  cat("\n  Spot-check: R CSV vs in-memory char_results (round-trip precision)\n")
  cat(sprintf("  %-5s  %-30s  %-12s  %-12s  %-6s\n",
              "Row", "Column", "CSV_value", "Matrix_val", "Match?"))
  cat("  ", strrep("-", 68), "\n")

  all_spot_ok <- TRUE
  for (i in spot_rows) {
    for (k in seq_along(spot_cols)) {
      j   <- spot_cols[k]
      mat_val <- r_full$char_results[i, j]
      csv_val <- suppressWarnings(as.numeric(r_csv[i, j]))
      if (is.na(mat_val) && is.na(csv_val)) {
        match_str <- "NA=NA"
      } else if (is.na(mat_val) || is.na(csv_val)) {
        match_str <- "MISMATCH"
        all_spot_ok <- FALSE
      } else if (abs(mat_val) < 1e-15) {
        match_str <- if (abs(csv_val) < 1e-10) "OK" else "FAIL"
        if (match_str == "FAIL") all_spot_ok <- FALSE
      } else {
        rel_diff <- abs(csv_val - mat_val) / abs(mat_val)
        match_str <- if (rel_diff < tol) "OK" else sprintf("FAIL(%.1e)", rel_diff)
        if (rel_diff >= tol) all_spot_ok <- FALSE
      }
      cat(sprintf("  %-5d  %-30s  %-12.6g  %-12.6g  %-6s\n",
                  i, substr(spot_cnames[k], 1, 29),
                  ifelse(is.na(csv_val), NA_real_, csv_val),
                  ifelse(is.na(mat_val), NA_real_, mat_val),
                  match_str))
    }
  }

  if (all_spot_ok) {
    cat("\n  PASS: All spot-check values round-trip correctly\n")
  } else {
    cat("\n  WARN: Some spot-check values differ beyond tolerance\n")
  }
}

# =============================================================================
# DONE
# =============================================================================

cat("\n")
cat(strrep("*", 72), "\n")
cat(" Phase 4 validation complete for Code Lake (CO).\n")
cat(strrep("*", 72), "\n\n")
