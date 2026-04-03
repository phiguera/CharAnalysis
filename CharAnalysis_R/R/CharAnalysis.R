#' Run the full CharAnalysis pipeline
#'
#' Top-level wrapper that calls each analytical stage in sequence and returns
#' all intermediate and final results as a single named list.
#'
#' Phases 1 and 2 (parameter reading, pretreatment, smoothing, C_peak
#' computation, and threshold determination) are implemented.  Phase 3
#' (peak identification, post-processing) and Phase 4 (output, figures)
#' will be added in subsequent phases.
#'
#' @param file_name Path to the \code{*_charParams.csv} (or \code{.xlsx})
#'   parameter file.
#'
#' @return Named list with the following elements (Phases 1–2):
#'   \describe{
#'     \item{charcoal}{List of raw and resampled series.  After Phase 2:
#'       also includes \code{accIS} (smoothed background) and \code{peak}
#'       (C_peak, either residuals or ratios).}
#'     \item{pretreatment}{Pretreatment parameter list (possibly updated by
#'       [char_pretreatment()] — e.g. \code{yrInterp} auto-set,
#'       \code{zoneDiv} end-value corrected).}
#'     \item{smoothing}{Smoothing parameter list.}
#'     \item{peak_analysis}{Peak-analysis parameter list.}
#'     \item{results}{Results / output parameter list.}
#'     \item{site}{Character string: site name.}
#'     \item{gap_in}{Integer matrix (nGaps x 2) of missing-value gap indices.}
#'     \item{char_thresh}{Threshold list returned by [char_thresh_global()] or
#'       [char_thresh_local()].  Contains \code{pos}, \code{neg}, \code{SNI},
#'       and \code{GOF}.}
#'   }
#'   Phases 3–4 will add \code{post}, etc.
#'
#' @seealso [char_parameters()], [char_validate_params()],
#'   [char_pretreatment()], [char_smooth()], [char_thresh_global()],
#'   [char_thresh_local()]
#'
#' @examples
#' \dontrun{
#'   out <- CharAnalysis("CO_charParams.csv")
#'   head(data.frame(ageTop_i = out$charcoal$ybpI,
#'                   charAcc_i = out$charcoal$accI,
#'                   charBkg_i = out$charcoal$accIS,
#'                   charPeak_i = out$charcoal$peak))
#' }
CharAnalysis <- function(file_name) {

  # (1) Read input file ---------------------------------------------------------
  message("(1) Reading input file...")
  params <- char_parameters(file_name)
  message("      ...done.")

  # (1b) Validate parameters ----------------------------------------------------
  message("(1b) Validating input parameters...")
  char_validate_params(params$char_data,
                       params$pretreatment,
                       params$smoothing,
                       params$peak_analysis,
                       params$results)

  # (2) Pretreatment ------------------------------------------------------------
  message("(2) Pretreating charcoal data...")
  pre <- char_pretreatment(params$char_data,
                            params$site,
                            params$pretreatment,
                            params$results,
                            plot_data = 0L)
  message("      ...done.")

  # (3) Smooth to estimate low-frequency C_background ---------------------------
  # Mirrors CharAnalysis.m step (3): CharSmooth()
  message("(3) Smoothing resampled CHAR to estimate low-frequency trends")
  message("    and calculating peak CHAR...")
  charcoal <- char_smooth(pre$charcoal,
                           pre$pretreatment,
                           params$smoothing,
                           params$results,
                           plot_data = 0L)

  # Guard: cannot compute ratio C_peak when background contains a zero.
  # Mirrors CharAnalysis.m lines 119-121.
  if (!is.null(charcoal$accIS) &&
      any(!is.na(charcoal$accIS)) &&
      min(charcoal$accIS, na.rm = TRUE) == 0 &&
      params$peak_analysis$cPeak == 2L) {
    stop("Cannot calculate C_peak (ratios) when C_background = 0; ",
         "change smoothing or cPeak parameters.")
  }

  # (3b) Compute peak CHAR (C_peak) --------------------------------------------
  # cPeak == 1 → residuals (accI - accIS)
  # cPeak == 2 → ratios    (accI / accIS)
  # Mirrors CharAnalysis.m lines 124-128.
  if (params$peak_analysis$cPeak == 1L) {
    charcoal$peak <- charcoal$accI - charcoal$accIS   # residuals
  } else {
    charcoal$peak <- charcoal$accI / charcoal$accIS   # ratios
  }
  message("      ...done.")

  # (4) Define thresholds -------------------------------------------------------
  # Mirrors CharAnalysis.m lines 131-141.
  message("(4) Defining possible thresholds for peak identification...")

  if (params$peak_analysis$threshType == 1L) {
    # Global threshold: one distribution fitted to the full C_peak record
    char_thresh <- char_thresh_global(charcoal,
                                       pre$pretreatment,
                                       params$peak_analysis,
                                       params$site,
                                       params$results,
                                       plot_data   = 0L,
                                       bkg_sens_in = 0L)
  } else {
    # Local threshold: per-sample sliding-window distribution
    char_thresh <- char_thresh_local(charcoal,
                                      params$smoothing,
                                      params$peak_analysis,
                                      params$site,
                                      params$results,
                                      plot_data = 0L)
  }
  message("      ...done.")

  # Assemble and return ---------------------------------------------------------
  list(
    charcoal      = charcoal,
    pretreatment  = pre$pretreatment,   # may differ from params$pretreatment
    smoothing     = params$smoothing,
    peak_analysis = params$peak_analysis,
    results       = params$results,
    site          = params$site,
    gap_in        = pre$gap_in,
    char_thresh   = char_thresh
  )
}
