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
#'   **Parameter parsing (label-anchored).** Each parameter is located by the
#'   label in the "Variable" column (column 2), not by a fixed row position.
#'   A multi-row parameter (e.g. \code{zoneDiv}, \code{threshValues}) is written
#'   on its labeled row plus any number of continuation rows whose "Variable"
#'   cell is blank; the label is carried down across those blanks.  All values
#'   for a given label are taken from the "Parameters" column (column 3).  This
#'   means the number of \code{zoneDiv} rows is free to vary, and the values
#'   below it are no longer shifted by the size of the \code{zoneDiv} block,
#'   which was the failure mode of the earlier position-based reader.  Unused
#'   \code{zoneDiv} slots may be left blank or set to the sentinel \code{-9999}
#'   (CSV) or \code{NaN} (Excel); all are stripped before the list is returned.
#'
#'   Required scalar parameters must each resolve to exactly one numeric value.
#'   If any required label is missing, misspelled, or has a non-numeric
#'   "Parameters" cell, the function stops with a message naming the offending
#'   row(s) rather than letting an \code{NA} propagate to a downstream check.
#'
#' @seealso [char_validate_params()], [char_pretreatment()], [CharAnalysis()]
#'
#' @export
#'
#' @examples
#' # Read a bundled example parameter file:
#' params_file <- system.file("validation", "CO_charParams.csv",
#'                            package = "CharAnalysis")
#' p <- char_parameters(params_file)
#' p$pretreatment$yrInterp   # interpolation interval (yr)
#' p$smoothing$method        # smoothing method index (1-5)
char_parameters <- function(file_name) {

  ext <- tolower(tools::file_ext(file_name))

  # =========================================================================
  # READ RAW PARAMETER TABLE
  # Each branch produces two vectors of equal length:
  #   var_col : the "Variable" labels (column 2)
  #   val_col : the numeric "Parameters" values (column 3)
  # plus char_data and site.
  # =========================================================================
  if (ext %in% c("xls", "xlsx")) {

    # -- Charcoal data sheet --------------------------------------------------
    char_data_raw <- readxl::read_excel(file_name, sheet = "charData",
                                        col_names = TRUE)
    char_data        <- as.matrix(char_data_raw[, seq_len(6L)])
    storage.mode(char_data) <- "double"

    # -- Parameter sheet: columns B (Variable) and C (Parameters), rows 2-26 --
    params_raw <- readxl::read_excel(file_name, sheet = "charParams",
                                     range = "B2:C26", col_names = FALSE)
    var_col <- as.character(params_raw[[1L]])
    val_col <- suppressWarnings(as.numeric(params_raw[[2L]]))

    # -- Site name from cell G1 of the charData sheet -------------------------
    site_raw <- readxl::read_excel(file_name, sheet = "charData",
                                   range = "G1:G1", col_names = FALSE)
    site <- if (!is.na(site_raw[[1L]][1L])) {
      as.character(site_raw[[1L]][1L])
    } else {
      warning("char_parameters: site name not found in cell G1 of charData ",
              "sheet. Using 'UnknownSite'.")
      "UnknownSite"
    }

  } else {

    # -- Site name: strip the trailing '_charParams.csv' (15 chars) -----------
    base_name <- basename(file_name)
    site      <- substr(base_name, 1L, nchar(base_name) - 15L)
    dir_path  <- dirname(file_name)

    # -- Companion charcoal data file -----------------------------------------
    data_file <- file.path(dir_path, paste0(site, "_charData.csv"))
    if (!file.exists(data_file)) {
      stop("char_parameters: companion data file not found: ", data_file)
    }
    char_data        <- as.matrix(read.csv(data_file, header = TRUE,
                                           stringsAsFactors = FALSE))
    storage.mode(char_data) <- "double"

    # -- Parameter file: column 2 (Variable), column 3 (Parameters) -----------
    # suppressWarnings because blank/text value cells become NA via as.numeric.
    params_df <- read.csv(file_name, header = TRUE, stringsAsFactors = FALSE)
    var_col   <- as.character(params_df[[2L]])
    val_col   <- suppressWarnings(as.numeric(params_df[[3L]]))
  }

  # =========================================================================
  # LABEL-ANCHORED LOOKUP
  # Carry each label down over blank continuation rows, then fetch all values
  # belonging to a given label.
  # =========================================================================
  var_col <- trimws(var_col)
  var_col[is.na(var_col)] <- ""
  for (i in seq_along(var_col)) {
    if (i > 1L && var_col[i] == "") var_col[i] <- var_col[i - 1L]
  }

  getv <- function(label) val_col[var_col == label]

  # =========================================================================
  # FAIL FAST ON MISSING / UNREADABLE REQUIRED PARAMETERS
  # Each required label must resolve to exactly one numeric value.  A missing
  # label yields length 0; a blank or non-numeric cell yields NA.  Catch both
  # here, with the offending row(s) named, instead of letting NA reach a
  # downstream if() and trigger "missing value, need TRUE/FALSE".
  # =========================================================================
  required <- list(
    yrInterpolate           = getv("yrInterpolate"),
    transform               = getv("transform"),
    method                  = getv("method"),
    yr                      = getv("yr"),
    cPeak                   = getv("cPeak"),
    threshType              = getv("threshType"),
    threshMethod            = getv("threshMethod"),
    minCountP               = getv("minCountP"),
    peakFrequ               = getv("peakFrequ"),
    `Cbackground sensitivity` = getv("Cbackground sensitivity"),
    saveFigures             = getv("saveFigures"),
    saveData                = getv("saveData")
  )
  bad <- names(required)[vapply(required,
                                function(x) length(x) == 0L || is.na(x[1L]),
                                logical(1L))]

  thresh_values <- getv("threshValues")
  thresh_values <- thresh_values[!is.na(thresh_values)]
  if (length(thresh_values) == 0L) bad <- c(bad, "threshValues")

  if (length(bad) > 0L) {
    stop("char_parameters: could not read a numeric value for parameter ",
         "row(s): ", paste(bad, collapse = ", "), ".\n",
         "  Check that each label in the 'Variable' column matches the ",
         "template exactly\n  and that its 'Parameters' cell contains a ",
         "number.", call. = FALSE)
  }

  # =========================================================================
  # UNPACK INTO NAMED LISTS
  # =========================================================================

  # -- Pretreatment ---------------------------------------------------------
  zone_div <- getv("zoneDiv")
  zone_div <- zone_div[!is.na(zone_div) & zone_div != -9999]  # strip sentinels

  pretreatment <- list(
    zoneDiv   = zone_div,
    yrInterp  = getv("yrInterpolate")[1L],
    transform = getv("transform")[1L]
  )

  # -- Smoothing ------------------------------------------------------------
  smoothing <- list(
    method = getv("method")[1L],
    yr     = getv("yr")[1L]
  )

  # -- PeakAnalysis ---------------------------------------------------------
  # Keep up to four threshold values; the last is used for peak plotting and
  # analysis, matching the MATLAB convention.
  thresh_values <- thresh_values[seq_len(min(4L, length(thresh_values)))]

  peak_analysis <- list(
    cPeak        = getv("cPeak")[1L],
    threshType   = getv("threshType")[1L],
    threshMethod = getv("threshMethod")[1L],
    threshValues = thresh_values,
    minCountP    = getv("minCountP")[1L],
    peakFrequ    = getv("peakFrequ")[1L],
    bkgSens      = getv("Cbackground sensitivity")[1L]
  )

  # -- Results --------------------------------------------------------------
  # allFigures is optional; default to 1 (show all diagnostic figures).
  all_figs_raw <- getv("allFigures")
  all_figs <- if (length(all_figs_raw) >= 1L && !is.na(all_figs_raw[1L])) {
    all_figs_raw[1L]
  } else {
    1
  }

  results <- list(
    saveFigures = getv("saveFigures")[1L],
    save        = getv("saveData")[1L],
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
