# =============================================================================
# test_bkg_sensitivity.R
# Validation of char_bkg_sensitivity() against the MATLAB v2.0
# bkgCharSensitivity.m (MATLAB Figure 10; R Figure 5).
#
# Datasets:
#   CO local   (Code Lake, threshType = 2; 2x2 diagnostic path)
#   CO global  (threshType = 1; contour path)  -- built in tempdir()
#
# Quantities validated:
#   z      peak-count surface (global) / per-threshValue columns (local)
#   SNI_i  signal-to-noise index (scalar per window global; per-sample local)
#   GOF_i  goodness-of-fit p-value (sentinel global; per-sample local)
#   mFRI   mean / SD of fire-return intervals (local only)
#
# Run from the CharAnalysis_2_0_R/ directory:
#   source("tests/test_bkg_sensitivity.R")
# Or set the root explicitly:
#   CHAR_R_ROOT <- "/path/to/CharAnalysis_2_0_R"
#   source(file.path(CHAR_R_ROOT, "tests/test_bkg_sensitivity.R"))
#
# ---------------------------------------------------------------------------
# VALIDATION APPROACH (two tiers, mirroring the phase tests)
#
# SECTION A -- Structural invariants (strict PASS/FAIL).  Properties that must
#   hold for the function to be correct, independent of GMM floating-point
#   divergence between R and MATLAB:
#     dimensions of z / SNI_i / GOF_i / mFRI; window grid in 100-yr steps;
#     peak counts non-negative integers; GOF in [0,1]; figure object built.
#
# SECTION B -- Numeric comparison vs MATLAB references (informational).  Read
#   the CSVs written by CharAnalysis_2_0_MATLAB/z_dump_bkgSens_refs.m if present
#   and report agreement.  Because peak identities can differ slightly due to
#   irreducible EM/GMM divergence (documented in the phase tests and validation
#   reports), this tier reports differences rather than failing the run.  For
#   the global path the candidate-threshold grids are constructed differently
#   (R: baseline C_peak range; MATLAB: Figure 2 axis limits), so the z surfaces
#   may have different widths; only the directly comparable SNI_i is checked
#   there.
# ---------------------------------------------------------------------------

if (!exists("CHAR_R_ROOT")) CHAR_R_ROOT <- "."
if (!dir.exists(file.path(CHAR_R_ROOT, "R"))) {
  stop("Set CHAR_R_ROOT to the CharAnalysis_2_0_R directory (found no R/ ",
       "subfolder relative to '", CHAR_R_ROOT, "').")
}

# ---- source package sources in dependency order (snake_case) ----------------
.src <- function(f) source(file.path(CHAR_R_ROOT, "R", f))
.src("char_lowess.R")
.src("gaussian_mixture.R")
.src("char_parameters.R")
.src("char_validate_params.R")
.src("char_pretreatment.R")
.src("char_smooth.R")
.src("char_thresh_global.R")
.src("char_thresh_local.R")
.src("char_peak_id.R")
.src("char_post_process.R")
.src("char_write_results.R")
.src("char_plot_results.R")     # internal plot helpers used by the figure
.src("CharAnalysis.R")
.src("char_bkg_sensitivity.R")

# ---- tiny test harness ------------------------------------------------------
.n_pass <- 0L; .n_fail <- 0L
chk <- function(label, cond) {
  ok <- isTRUE(cond)
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (ok) .n_pass <<- .n_pass + 1L else .n_fail <<- .n_fail + 1L
}

co_params <- file.path(CHAR_R_ROOT, "tests", "CO_charParams.csv")
co_data    <- file.path(CHAR_R_ROOT, "tests", "CO_charData.csv")
stopifnot(file.exists(co_params), file.exists(co_data))

# =============================================================================
# Run the analyses
# =============================================================================
cat("\n=== Running CO local (threshType = 2) ===\n")
out_local  <- suppressMessages(CharAnalysis(co_params))
sens_local <- char_bkg_sensitivity(out_local, verbose = FALSE)

cat("\n=== Running CO global (threshType = 1, tempdir variant) ===\n")
tdir <- tempfile("co_global_"); dir.create(tdir)
file.copy(co_data, file.path(tdir, "COg_charData.csv"), overwrite = TRUE)
pl <- readLines(co_params)
pl <- sub(",threshType,2,", ",threshType,1,", pl, fixed = TRUE)
writeLines(pl, file.path(tdir, "COg_charParams.csv"))
out_global  <- suppressMessages(CharAnalysis(file.path(tdir, "COg_charParams.csv")))
sens_global <- char_bkg_sensitivity(out_global, verbose = FALSE)

have_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

# =============================================================================
# SECTION A -- Structural invariants
# =============================================================================
cat("\n--- A. Structural invariants ---\n")

## ---- Local ----
cat(" Local (threshType = 2):\n")
N_l   <- length(out_local$charcoal$peak)
nw_l  <- length(sens_local$bkSmooth)
n_tv  <- length(out_local$peak_analysis$threshValues)

chk("bkSmooth has >= 2 windows", nw_l >= 2L)
chk("bkSmooth in ascending 100-yr steps",
    all(diff(sens_local$bkSmooth) == 100) && all(sens_local$bkSmooth %% 100 == 0))
chk("z is [nWin x nThreshValues]",
    is.matrix(sens_local$z) && all(dim(sens_local$z) == c(nw_l, n_tv)))
chk("SNI_i is [nWin x N]",
    is.matrix(sens_local$SNI_i) && all(dim(sens_local$SNI_i) == c(nw_l, N_l)))
chk("GOF_i is [nWin x N]",
    is.matrix(sens_local$GOF_i) && all(dim(sens_local$GOF_i) == c(nw_l, N_l)))
chk("mFRI is [nWin x 2]",
    is.matrix(sens_local$mFRI) && all(dim(sens_local$mFRI) == c(nw_l, 2L)))
chk("z non-negative whole numbers",
    all(sens_local$z >= 0, na.rm = TRUE) &&
    all(sens_local$z == round(sens_local$z), na.rm = TRUE))
chk("GOF_i within [0, 1] (non-NA)",
    all(sens_local$GOF_i >= 0 & sens_local$GOF_i <= 1, na.rm = TRUE))
chk("SNI_i non-negative (non-NA)", all(sens_local$SNI_i >= 0, na.rm = TRUE))
if (have_ggplot)
  chk("figure object built", inherits(sens_local$fig, "ggplot") ||
                             inherits(sens_local$fig, "patchwork"))

## ---- Global ----
cat(" Global (threshType = 1):\n")
nw_g    <- length(sens_global$bkSmooth)
n_posb  <- sum(sens_global$thresh_bins > 0)

chk("bkSmooth has >= 2 windows", nw_g >= 2L)
chk("bkSmooth starts at 100 yr (global bkgMin)", sens_global$bkSmooth[1] == 100)
chk("thresh_bins has 251 candidate values",
    length(sens_global$thresh_bins) == 251L)
chk("z is [nWin x nPositiveBins]",
    is.matrix(sens_global$z) && all(dim(sens_global$z) == c(nw_g, n_posb)))
chk("SNI_i is length nWin", length(sens_global$SNI_i) == nw_g &&
                            is.null(dim(sens_global$SNI_i)))
chk("mFRI is NULL on global path", is.null(sens_global$mFRI))
chk("z non-negative whole numbers",
    all(sens_global$z >= 0, na.rm = TRUE) &&
    all(sens_global$z == round(sens_global$z), na.rm = TRUE))
if (have_ggplot)
  chk("figure object built", inherits(sens_global$fig, "ggplot") ||
                             inherits(sens_global$fig, "patchwork"))

# =============================================================================
# SECTION B -- Numeric comparison vs MATLAB references (informational)
# =============================================================================
cat("\n--- B. Comparison vs MATLAB references (informational) ---\n")

refs_dir <- if (!is.null(getOption("CHAR_MATLAB_REFS"))) {
  getOption("CHAR_MATLAB_REFS")
} else if (nzchar(Sys.getenv("CHAR_MATLAB_REFS"))) {
  Sys.getenv("CHAR_MATLAB_REFS")
} else {
  file.path(CHAR_R_ROOT, "..", "CharAnalysis_2_0_MATLAB", "bkgSens_refs")
}

read_ref <- function(name) {
  f <- file.path(refs_dir, name)
  if (!file.exists(f)) return(NULL)
  as.matrix(utils::read.csv(f, header = FALSE))
}

report_cmp <- function(label, r_mat, m_mat) {
  if (is.null(m_mat)) { cat(sprintf("  (skip) %s: no reference file\n", label)); return(invisible()) }
  r_mat <- as.matrix(r_mat); m_mat <- as.matrix(m_mat)
  if (!all(dim(r_mat) == dim(m_mat))) {
    cat(sprintf("  (info) %s: dim differ R[%s] vs MATLAB[%s] -- not elementwise comparable\n",
                label, paste(dim(r_mat), collapse="x"), paste(dim(m_mat), collapse="x")))
    return(invisible())
  }
  d   <- abs(r_mat - m_mat)
  rng <- range(c(r_mat, m_mat), na.rm = TRUE)
  cat(sprintf("  %s: max|diff| = %.4g, mean|diff| = %.4g (value range %.3g..%.3g)\n",
              label, max(d, na.rm = TRUE), mean(d, na.rm = TRUE), rng[1], rng[2]))
}

if (!dir.exists(refs_dir)) {
  cat(sprintf("  No reference directory at: %s\n", normalizePath(refs_dir, mustWork = FALSE)))
  cat("  Generate it with CharAnalysis_2_0_MATLAB/z_dump_bkgSens_refs.m, then\n")
  cat("  copy bkgSens_refs/ next to the MATLAB sources (or set CHAR_MATLAB_REFS).\n")
} else {
  cat(sprintf("  Reference directory: %s\n", refs_dir))
  cat(" Local:\n")
  report_cmp("z      (local)", sens_local$z,     read_ref("CO_local_z.csv"))
  report_cmp("SNI_i  (local)", sens_local$SNI_i, read_ref("CO_local_SNI_i.csv"))
  report_cmp("GOF_i  (local)", sens_local$GOF_i, read_ref("CO_local_GOF_i.csv"))
  cat(" Global (SNI_i directly comparable; z grids constructed differently):\n")
  report_cmp("SNI_i  (global)", matrix(sens_global$SNI_i, ncol = 1),
             read_ref("CO_global_SNI_i.csv"))
  report_cmp("z      (global)", sens_global$z, read_ref("CO_global_z.csv"))
}

# =============================================================================
# Summary
# =============================================================================
cat(sprintf("\n=== Section A: %d passed, %d failed ===\n", .n_pass, .n_fail))
if (.n_fail > 0L) stop("Structural invariants failed; see FAIL lines above.")
cat("Section A passed. Section B is informational (see above).\n")
