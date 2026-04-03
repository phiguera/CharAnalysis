#' Read CharAnalysis parameter and data files
#'
#' Reads the \code{*_charParams.csv} (or \code{.xlsx}) parameter file and the
#' companion charcoal data file, then unpacks all analysis parameters into
#' named lists that mirror the MATLAB struct layout.
#'
#' @param file_name Path to the \code{*_charParams.csv} (or \code{.xlsx})
#'   parameter file.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{char_data}{Numeric matrix (n_samples x 6+): cmTop, cmBot, ageTop,
#'       ageBot, charVol, charCount.}
#'     \item{pretreatment}{List: \code{zoneDiv}, \code{yrInterp},
#'       \code{transform}.}
#'     \item{smoothing}{List: \code{method}, \code{yr}.}
#'     \item{peak_analysis}{List: \code{cPeak}, \code{threshType},
#'       \code{threshMethod}, \code{threshValues}, \code{minCountP},
#'       \code{peakFrequ}, \code{bkgSens}.}
#'     \item{results}{List: \code{saveFigures}, \code{save},
#'       \code{allFigures}.}
#'     \item{site}{Character string: site name derived from the filename stem.}
#'   }
#'
#' @details
#'   **CSV convention** (unchanged from MATLAB v1.1): the parameter file is
#'   named \code{<site>_charParams.csv} and the companion data file is
#'   \code{<site>_charData.csv} in the same directory.  The site name is the
#'   basename with the trailing \code{_charParams.csv} (15 characters) removed,
#'   mirroring MATLAB's \code{fileName(1:end-15)} idiom.
#'
#'   **Parameter vector layout**: column 3 ("Parameters") of the CSV, rows
#'   2–26 (25 data rows after the header), maps to positions 1–25 of the
#'   internal \code{charParams} vector exactly as in the MATLAB codebase.
#'   Unused \code{zoneDiv} slots are filled with \code{-9999} in the CSV
#'   (sentinel) or \code{NaN} in Excel; both are stripped before the list is
#'   returned.
#'
#' @seealso [char_validate_params()], [char_pretreatment()], [CharAnalysis()]
#'
#' @examples
#' \dontrun{
#'   p <- char_parameters("CO_charParams.csv")
#'   p$pretreatment$yrInterp   # interpolation interval (yr)
#'   p$smoothing$method        # smoothing method index (1–5)
#' }
char_parameters <- function(file_name) {

  ext <- tolower(tools::file_ext(file_name))

  # =========================================================================
  # EXCEL PATH (.xls / .xlsx)
  # =========================================================================
  if (ext %in% c("xls", "xlsx")) {

    # Charcoal data sheet
    char_data_raw <- readxl::read_excel(file_name, sheet = "charData",
                                        col_names = TRUE)
    char_data        <- as.matrix(char_data_raw[, seq_len(6L)])
    storage.mode(char_data) <- "double"

    # Parameter column: rows 2–26, column C (numeric values)
    params_raw  <- readxl::read_excel(file_name, sheet = "charParams",
                                      range = "C2:C26", col_names = FALSE)
    char_params <- as.numeric(unlist(params_raw))

    # Site name from cell G1 of the charData sheet
    site_raw <- readxl::read_excel(file_name, sheet = "charData",
                                   range = "G1:G1", col_names = FALSE)
    site <- if (!is.na(site_raw[[1L]][1L])) {
      as.character(site_raw[[1L]][1L])
    } else {
      warning("char_parameters: site name not found in cell G1 of charData ",
              "sheet. Using 'UnknownSite'.")
      "UnknownSite"
    }

  # =========================================================================
  # CSV PATH
  # =========================================================================
  } else {

    # Site name: strip the trailing '_charParams.csv' (15 chars) from the
    # base filename, matching MATLAB's  site = fileName(1:end-15)  idiom.
    base_name <- basename(file_name)
    site      <- substr(base_name, 1L, nchar(base_name) - 15L)
    dir_path  <- dirname(file_name)

    # Companion charcoal data file
    data_file <- file.path(dir_path, paste0(site, "_charData.csv"))
    if (!file.exists(data_file)) {
      stop("char_parameters: companion data file not found: ", data_file)
    }
    char_data        <- as.matrix(read.csv(data_file, header = TRUE,
                                           stringsAsFactors = FALSE))
    storage.mode(char_data) <- "double"

    # Parameter file: skip header, read 25 data rows.
    # Column 3 ("Parameters") contains the numeric values; all other columns
    # are discarded.  suppressWarnings because some rows may have empty strings
    # in the numeric column (handled cleanly by as.numeric -> NA).
    params_df   <- read.csv(file_name, header = TRUE,
                            stringsAsFactors = FALSE, nrows = 25L)
    char_params <- suppressWarnings(as.numeric(params_df[[3L]]))
  }

  # =========================================================================
  # UNPACK PARAMETER VECTOR INTO NAMED LISTS
  # Positions 1–25 match the charParams worksheet rows 2–26 in the .xlsx
  # template, identical to the MATLAB layout.
  # =========================================================================

  # -- Pretreatment (positions 1–10) ----------------------------------------
  zone_div <- char_params[1:8]
  zone_div <- zone_div[!is.na(zone_div) & zone_div != -9999]  # strip sentinels

  pretreatment <- list(
    zoneDiv   = zone_div,
    yrInterp  = char_params[9L],
    transform = char_params[10L]
  )

  # -- Smoothing (positions 11–12) ------------------------------------------
  smoothing <- list(
    method = char_params[11L],
    yr     = char_params[12L]
  )

  # -- PeakAnalysis (positions 13–22) ---------------------------------------
  peak_analysis <- list(
    cPeak        = char_params[13L],
    threshType   = char_params[14L],
    threshMethod = char_params[15L],
    threshValues = char_params[16:19],
    minCountP    = char_params[20L],
    peakFrequ    = char_params[21L],
    bkgSens      = char_params[22L]
  )

  # -- Results (positions 23–25) --------------------------------------------
  all_figs <- if (length(char_params) >= 25L && !is.na(char_params[25L])) {
    char_params[25L]
  } else {
    1  # default: show all diagnostic figures
  }

  results <- list(
    saveFigures = char_params[23L],
    save        = char_params[24L],
    allFigures  = all_figs
  )

  # =========================================================================
  # RETURN
  # =========================================================================
  list(
    char_data     = char_data,
    pretreatment  = pretreatment,
    smoothing     = smoothing,
    peak_analysis = peak_analysis,
    results       = results,
    site          = site
  )
}
