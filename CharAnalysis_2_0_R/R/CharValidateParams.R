#' Validate CharAnalysis input parameters
#'
#' Checks all user-supplied parameters for internal consistency and plausible
#' ranges.  Throws a descriptive error and halts if any check fails.  Emits
#' warnings for non-fatal but unusual settings.
#'
#' Mirrors \code{CharValidateParams.m} from the MATLAB v2.0 codebase.  All 15
#' checks and 2 non-fatal warnings are reproduced in the same order.
#'
#' @param char_data    Numeric matrix (n x 6+).
#' @param pretreatment List with \code{zoneDiv}, \code{yrInterp},
#'   \code{transform}.
#' @param smoothing    List with \code{method}, \code{yr}.
#' @param peak_analysis List with \code{cPeak}, \code{threshType},
#'   \code{threshMethod}, \code{threshValues}, \code{minCountP},
#'   \code{peakFrequ}.
#' @param results      List (not checked; included for API symmetry).
#'
#' @return \code{NULL} invisibly on success.  Stops with an error on failure.
#'
#' @seealso [char_parameters()], [CharAnalysis()]
char_validate_params <- function(char_data, pretreatment, smoothing,
                                  peak_analysis, results) {

  # 1. At least 6 columns ---------------------------------------------------
  if (ncol(char_data) < 6L) {
    stop("char_validate_params: char_data must have at least 6 columns ",
         "(cmTop, cmBot, ageTop, ageBot, charVol, charCount). Found ",
         ncol(char_data), ".")
  }

  # 2. No NA ages ------------------------------------------------------------
  if (anyNA(char_data[, 3L]) || anyNA(char_data[, 4L])) {
    stop("char_validate_params: char_data contains NA values in the age ",
         "columns (cols 3-4). All sample ages must be specified.")
  }

  # 3. Ages monotonically non-decreasing ------------------------------------
  age_tops <- char_data[, 3L]
  if (any(diff(age_tops) < 0)) {
    bad_row <- which(diff(age_tops) < 0)[1L]
    stop("char_validate_params: ageTop values (col 3) are not monotonically ",
         "non-decreasing. First violation at row ", bad_row, ".")
  }

  # 4. zoneDiv at least 2 values --------------------------------------------
  if (length(pretreatment$zoneDiv) < 2L) {
    stop("char_validate_params: zoneDiv must contain at least 2 values ",
         "(start and end of record).")
  }

  # 5. zoneDiv strictly ascending -------------------------------------------
  if (any(diff(pretreatment$zoneDiv) <= 0)) {
    stop("char_validate_params: zoneDiv must be strictly ascending ",
         "(youngest age first).")
  }

  # 6. yrInterp non-negative ------------------------------------------------
  if (pretreatment$yrInterp < 0) {
    stop("char_validate_params: yrInterp must be >= 0 ",
         "(0 = use median sample resolution).")
  }

  # 7. Smoothing window shorter than record ---------------------------------
  record_len <- max(pretreatment$zoneDiv) - min(pretreatment$zoneDiv)
  smooth_yr  <- smoothing$yr[length(smoothing$yr)]   # longest window specified
  if (smooth_yr >= record_len) {
    stop("char_validate_params: smoothing window (", smooth_yr,
         " yr) must be shorter than the record length (", record_len, " yr).")
  }

  # 8. Smoothing method in range 1-5 ----------------------------------------
  if (!smoothing$method %in% 1:5) {
    stop("char_validate_params: smoothing method must be 1-5. Got ",
         smoothing$method, ".")
  }

  # 9. threshType valid (1 or 2) --------------------------------------------
  if (!peak_analysis$threshType %in% c(1, 2)) {
    stop("char_validate_params: threshType must be 1 (global) or 2 (local). ",
         "Got ", peak_analysis$threshType, ".")
  }

  # 10. threshMethod valid (1, 2, or 3) -------------------------------------
  if (!peak_analysis$threshMethod %in% c(1, 2, 3)) {
    stop("char_validate_params: threshMethod must be 1, 2, or 3. Got ",
         peak_analysis$threshMethod, ".")
  }

  # 11. Local threshold cannot be user-defined ------------------------------
  if (peak_analysis$threshMethod == 1L && peak_analysis$threshType == 2L) {
    stop("char_validate_params: a locally defined threshold (threshType = 2) ",
         "cannot be user-defined (threshMethod = 1). Change input parameters.")
  }

  # 12. threshValues in (0, 1) when threshMethod > 1 ------------------------
  if (peak_analysis$threshMethod > 1L) {
    tv <- peak_analysis$threshValues[1:4]
    if (any(tv <= 0) || any(tv >= 1)) {
      stop("char_validate_params: all threshValues must be in the open ",
           "interval (0, 1) when threshMethod = ", peak_analysis$threshMethod,
           ". Got [", paste(tv, collapse = ", "), "].")
    }
  }

  # 13. cPeak valid (1 or 2) ------------------------------------------------
  if (!peak_analysis$cPeak %in% c(1, 2)) {
    stop("char_validate_params: cPeak must be 1 (residuals) or 2 (ratios). ",
         "Got ", peak_analysis$cPeak, ".")
  }

  # 14. minCountP in [0, 1] -------------------------------------------------
  if (peak_analysis$minCountP < 0 || peak_analysis$minCountP > 1) {
    stop("char_validate_params: minCountP must be in [0, 1]. Got ",
         peak_analysis$minCountP, ".")
  }

  # 15. peakFrequ positive --------------------------------------------------
  if (peak_analysis$peakFrequ <= 0) {
    stop("char_validate_params: peakFrequ must be positive. Got ",
         peak_analysis$peakFrequ, ".")
  }

  # -- Non-fatal warnings ---------------------------------------------------
  if (smooth_yr < 0.05 * record_len) {
    warning("char_validate_params: smoothing window (", smooth_yr,
            " yr) is less than 5% of the record length (", record_len,
            " yr). Consider a longer window.")
  }

  if (peak_analysis$threshType == 2L && pretreatment$yrInterp > 0) {
    samples_in_window <- smooth_yr / pretreatment$yrInterp
    if (samples_in_window < 30) {
      warning("char_validate_params: smoothing window (", smooth_yr,
              " yr) at yrInterp = ", pretreatment$yrInterp,
              " yr gives ~", floor(samples_in_window),
              " samples per local threshold window. ",
              "At least 30 are recommended.")
    }
  }

  message("      Parameter validation passed.")
  invisible(NULL)
}
