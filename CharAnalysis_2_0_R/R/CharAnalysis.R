#' Run the full CharAnalysis pipeline
#'
#' Top-level wrapper that calls each analytical stage in sequence and returns
#' all intermediate and final results as a single named list.
#'
#' @param file_name Path to the \code{*_charParams.csv} (or \code{.xlsx})
#'   parameter file.  If omitted (or \code{NULL}), a system file-picker dialog
#'   is opened so you can navigate to the file interactively -- matching the
#'   MATLAB behaviour of typing \code{CharAnalysis} with no argument.  In
#'   non-interactive sessions (e.g. scripts, batch jobs) the argument is
#'   required and an error is thrown if it is missing.
#'
#' @return Named list with the following elements:
#'   \describe{
#'     \item{charcoal}{List of raw and resampled series.  After Phase 2:
#'       includes \code{accIS} (smoothed background) and \code{peak}
#'       (C_peak, either residuals or ratios).  After Phase 3: also includes
#'       \code{charPeaks} (\eqn{[N \times T]} binary peak matrix),
#'       \code{charPeaksThresh}, \code{peaksTotal}, and \code{threshFRI}.
#'       After Phase 4: also includes \code{peakInsig},
#'       \code{peakMagnitude}, \code{smoothedFireFrequ}, \code{peaksFrequ}.}
#'     \item{pretreatment}{Pretreatment parameter list (possibly updated by
#'       [char_pretreatment()] -- e.g. \code{yrInterp} auto-set,
#'       \code{zoneDiv} end-value corrected).}
#'     \item{smoothing}{Smoothing parameter list.}
#'     \item{peak_analysis}{Peak-analysis parameter list.}
#'     \item{results}{Results / output parameter list.}
#'     \item{site}{Character string: site name.}
#'     \item{gap_in}{Integer matrix (nGaps x 2) of missing-value gap indices.}
#'     \item{char_thresh}{Threshold list returned by [char_thresh_global()] or
#'       [char_thresh_local()].  Contains \code{pos}, \code{neg}, \code{SNI},
#'       \code{GOF}, and (after Phase 3) \code{minCountP} -- the
#'       \eqn{[N \times T]} matrix of Shuie-Bain p-values.}
#'     \item{post}{Post-processing list from [char_post_process()]: FRI
#'       series, smoothed FRI curve, per-zone Weibull statistics, and the
#'       assembled \code{char_results} output matrix (\eqn{N \times 33}).}
#'     \item{char_results}{Numeric matrix (\eqn{N \times 33}) matching the
#'       MATLAB \code{charResults} output exactly (alias of
#'       \code{post$char_results}).}
#'     \item{bkg_sensitivity}{Present only when \code{peak_analysis$bkgSens == 1}
#'       in the parameter file.  The list returned by [char_bkg_sensitivity()]
#'       (window grid, peak-count surface, SNI/GOF, FRI summaries, and the
#'       figure).}
#'   }
#'
#' @seealso [char_parameters()], [char_validate_params()],
#'   [char_pretreatment()], [char_smooth()], [char_thresh_global()],
#'   [char_thresh_local()], [char_peak_id()], [char_post_process()],
#'   [char_bkg_sensitivity()]
#'
#' @examples
#' \donttest{
#'   # Run the full pipeline on the bundled example dataset:
#'   params_file <- system.file("validation", "CO_charParams.csv",
#'                              package = "CharAnalysis")
#'   out <- CharAnalysis(params_file)
#'   # Phase 2 outputs
#'   head(data.frame(ageTop_i   = out$charcoal$ybpI,
#'                   charAcc_i  = out$charcoal$accI,
#'                   charBkg_i  = out$charcoal$accIS,
#'                   charPeak_i = out$charcoal$peak))
#'   # Phase 3 outputs
#'   sum(out$charcoal$charPeaks[, ncol(out$charcoal$charPeaks)])
#' }
#' @export
CharAnalysis <- function(file_name = NULL) {

  # If no file supplied, open an interactive picker (mirrors MATLAB behaviour
  # of prompting when CharAnalysis is called with no argument).
  if (is.null(file_name)) {
    if (!interactive()) {
      stop("file_name is required in non-interactive sessions (scripts, ",
           "batch jobs, RMarkdown, etc.).")
    }
    message("No file specified -- opening file picker...")
    file_name <- tryCatch(
      file.choose(),
      error = function(e) {
        # file.choose() not available (e.g. RStudio Server, some terminals);
        # fall back to a console prompt.
        readline("Enter path to *_charParams.csv or .xlsx file: ")
      }
    )
    if (!nzchar(file_name)) stop("No file selected. CharAnalysis cancelled.")
  }

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

  # Figure 1 (allFigures only): C_raw / C_resampled / C_background options.
  # Mirrors MATLAB CharPretreatment.m subplot 1 + CharSmooth.m subplot 2.
  if (isTRUE(params$results$allFigures == 1L)) {
    mini_out <- list(charcoal     = charcoal,
                     pretreatment = pre$pretreatment,
                     smoothing    = params$smoothing,
                     site         = params$site)
    char_plot_raw(mini_out)
  }

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
  # cPeak == 1 -> residuals (accI - accIS)
  # cPeak == 2 -> ratios    (accI / accIS)
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

  # SNI advisory: warn when any portion of the record has SNI < 3.0 ----------
  # Kelly et al. (2011) show that SNI > 3 consistently identifies records
  # appropriate for peak detection; values below this threshold indicate
  # portions of the record where peak detection may be unreliable.
  # Stored here and emitted at the end of step 7 so it appears after all
  # other output. Placed here (not inside char_thresh_local/global) to avoid
  # repeated firing during char_bkg_sensitivity().
  sni_advisory <- NULL
  if (!is.null(char_thresh$SNI)) {
    if (params$peak_analysis$threshType == 1L) {
      if (!is.na(char_thresh$SNI) && char_thresh$SNI < 3.0) {
        sni_advisory <- sprintf(
          "The record-wide SNI is %.2f, below the minimum value of 3.0 recommended by Kelly et al. (2011) to identify records suitable for peak detection. Carefully consider whether and how to interpret this record.",
          char_thresh$SNI
        )
      }
    } else {
      valid_sni <- char_thresh$SNI[!is.na(char_thresh$SNI)]
      if (length(valid_sni) > 0L) {
        pct_below <- mean(valid_sni < 3.0) * 100
        if (pct_below > 0) {
          sni_advisory <- sprintf(
            "%.1f%% of samples have an SNI < 3.0, the minimum value recommended by Kelly et al. (2011) to identify records suitable for peak detection. Carefully consider whether and how to interpret portions of this record with SNI < 3.0.",
            pct_below
          )
        }
      }
    }
  }

  # Figure 2 (allFigures only): threshold determination diagnostics.
  # Mirrors MATLAB CharThreshGlobal.m (single panel) or
  # CharThreshLocal.m (5x5 grid of local window distributions).
  if (isTRUE(params$results$allFigures == 1L)) {
    mini_out2 <- list(charcoal      = charcoal,
                      char_thresh   = char_thresh,
                      peak_analysis = params$peak_analysis,
                      pretreatment  = pre$pretreatment,
                      site          = params$site)
    char_plot_thresh_diag(mini_out2)
  }

  # (5) Identify peaks ----------------------------------------------------------
  # Mirrors CharAnalysis.m step (5): CharPeakID()
  message("(5) Identifying peaks and applying minimum-count screening...")
  peak_result <- char_peak_id(charcoal,
                               pre$pretreatment,
                               params$peak_analysis,
                               char_thresh)
  charcoal    <- peak_result$charcoal
  char_thresh <- peak_result$char_thresh
  message("      ...done.")

  # (6) Post-processing: FRIs, fire frequency, Weibull statistics --------------
  # Mirrors CharAnalysis.m step (6): CharPostProcess()
  message("(6) Post-processing: fire-return intervals, Weibull statistics...")
  post_result  <- char_post_process(charcoal,
                                     pre$pretreatment,
                                     params$peak_analysis,
                                     char_thresh,
                                     params$smoothing)
  charcoal     <- post_result$charcoal
  post         <- post_result$post
  char_results <- post_result$char_results
  message("      ...done.")

  # (7) Write results CSV -------------------------------------------------------
  # In the R package, CSV output is explicit: call char_write_results() directly
  # after CharAnalysis() returns.  The saveData flag from the parameter file is
  # stored in results$save and can be inspected by the caller, but no file is
  # written automatically here (prevents accidental overwrites of reference data).
  message("(7) Analysis complete.")
  message("    Save CSV:     char_write_results(out$char_results, out$site, out_dir = \"<your/path>\")")
  message("    All figures:  char_plot_all(out)  [Figs 1-2 only when allFigures = 1]")
  message("                  char_plot_all(out, save = TRUE, out_dir = \"<your/path>\")  # save PDFs")
  message("    One figure:   char_plot_raw(out)            # Fig 1: C_raw, C_interp, C_back options")
  message("                  char_plot_thresh_diag(out)    # Fig 2: threshold diagnostics")
  message("                  char_plot_peaks(out)          # Fig 3: peak analysis")
  message("                  char_plot_sni(out)            # Fig 4: threshold sensitivity and SNI")
  message("                  char_plot_cumulative(out)     # Fig 5: cumulative peaks")
  message("                  char_plot_fri(out)            # Fig 6: FRI distributions")
  message("                  char_plot_fire_history(out)   # Fig 7: continuous fire history")
  message("                  char_plot_zones(out)          # Fig 8: CHAR zone comparisons")
  message("    Sensitivity:  char_bkg_sensitivity(out)     # Fig 10: C_background window sensitivity")

  # Emit SNI advisory last, after all other output, so it is not buried.
  if (!is.null(sni_advisory)) {
    message("")
    cli::cli_text(cli::col_red("CAUTION"))
    cli::cli_alert_danger(sni_advisory)
    cli::cli_alert_info("Run {.code char_plot_sni(out)} to visualize the SNI series and identify which portion(s) of the record fall below the threshold.")
  }

  # Assemble and return ---------------------------------------------------------
  out <- list(
    charcoal      = charcoal,
    pretreatment  = pre$pretreatment,   # may differ from params$pretreatment
    smoothing     = params$smoothing,
    peak_analysis = params$peak_analysis,
    results       = params$results,
    site          = params$site,
    gap_in        = pre$gap_in,
    char_thresh   = char_thresh,
    post          = post,
    char_results  = char_results
  )
  class(out) <- c("CharAnalysis", "list")

  # (7) Background-window sensitivity analysis (optional) -----------------------
  # Gated by the params file: peak_analysis$bkgSens (CSV row 23, position 22).
  # Mirrors CharAnalysis.m step (7). Numeric results are always computed when
  # bkgSens == 1; the figure is shown only when allFigures == 1 (matching how
  # Figs 1-2 are handled), and is never auto-saved (saving is explicit, like
  # CSV output and the other figures).
  if (isTRUE(out$peak_analysis$bkgSens == 1L)) {
    message("(7) Running C_background sensitivity analysis (bkgSens = 1)...")
    out$bkg_sensitivity <- char_bkg_sensitivity(out, save = FALSE,
                                                verbose = TRUE)
    if (isTRUE(out$results$allFigures == 1L) &&
        !is.null(out$bkg_sensitivity$fig)) {
      print(out$bkg_sensitivity$fig)
    }
    message("      ...done. Save the figure with: ",
            "char_bkg_sensitivity(out, save = TRUE, out_dir = \"<your/path>\")")
  }

  out
}

# ── S3 methods ────────────────────────────────────────────────────────────────

#' @export
print.CharAnalysis <- function(x, ...) {
  n      <- nrow(x$char_results)
  n_peak <- sum(x$charcoal$charPeaks[, ncol(x$charcoal$charPeaks)],
                na.rm = TRUE)
  zones  <- x$pretreatment$zoneDiv
  n_zone <- max(1L, length(zones) - 1L)
  cat(sprintf(
    "CharAnalysis results \u2014 %s\n  %d interpolated samples  |  %d fire event%s  |  %d zone%s\n",
    x$site,
    n,
    n_peak, if (n_peak == 1) "" else "s",
    n_zone, if (n_zone == 1) "" else "s"
  ))
  invisible(x)
}

#' @export
summary.CharAnalysis <- function(object, ...) {
  n_peak  <- sum(object$charcoal$charPeaks[, ncol(object$charcoal$charPeaks)],
                 na.rm = TRUE)
  age_rng <- range(object$charcoal$ybpI, na.rm = TRUE)
  fpi     <- object$post$FRI_params_zone
  cat(sprintf("Site:            %s\n", object$site))
  cat(sprintf("Record length:   %d samples  (%.0f \u2013 %.0f yr BP)\n",
              nrow(object$char_results), age_rng[1], age_rng[2]))
  cat(sprintf("Smoothing:       method %d  (%d yr window)\n",
              object$smoothing$method, object$smoothing$yr))
  cat(sprintf("Peaks (final):   %d\n", n_peak))
  if (!is.null(fpi) && nrow(fpi) > 0) {
    cat("Zone statistics (nFRI / mFRI yr):\n")
    for (z in seq_len(nrow(fpi))) {
      cat(sprintf("  Zone %d:  n=%d  mFRI=%.1f yr\n",
                  z, fpi[z, 1], fpi[z, 2]))
    }
  }
  invisible(object)
}

#' @export
plot.CharAnalysis <- function(x, ...) {
  char_plot_all(x, ...)
}