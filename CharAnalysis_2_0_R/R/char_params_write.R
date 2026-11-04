# Writing and revising *_charParams.csv files
#
# Internal helpers for the parameter interface (ROADMAP_ParamsGUI.md, Phase 1).
# Both the R reader (char_parameters()) and the MATLAB reader
# (CharParameters.m) locate parameters by row position in column 3, so every
# file written here reproduces the template layout exactly: one header row,
# 25 parameter rows, and the template's Stage / Variable / Units / Description
# text in columns 1, 2, 4, and 5.
#
# Column layout:
#   1-5   Stage, Variable, Parameters, Units, Description and Options
#   6     spacer (empty)
#   7-16  Parameter Log: previous parameter sets, newest in column 7
#
# The MATLAB reader parses a fixed 16 fields per row (textscan in
# CharParameters.m), so a file may hold at most 10 log columns.
#
# AI-assisted development (Claude/Anthropic).

.CP_N_PARAMS <- 25L   # parameter rows after the header
.CP_LOG_COL  <- 7L    # first Parameter Log column
.CP_MAX_COL  <- 16L   # fields per row parsed by MATLAB CharParameters.m
.CP_MAX_LOG  <- .CP_MAX_COL - .CP_LOG_COL + 1L   # 10 log columns

#' Read a parameter file as a character table (all rows and columns)
#' @param file Path to a *_charParams.csv file.
#' @return Character matrix; row 1 is the header, rows 2-26 the parameters.
#' @noRd
char_params_read_table <- function(file) {
  if (!file.exists(file)) stop("Parameter file not found: ", file, call. = FALSE)
  n_fields <- utils::count.fields(file, sep = ",", quote = "\"",
                                  blank.lines.skip = TRUE)
  n_col <- max(5L, n_fields, na.rm = TRUE)
  tab <- utils::read.csv(file, header = FALSE, colClasses = "character",
                         col.names = paste0("V", seq_len(n_col)),
                         fill = TRUE, na.strings = character(0),
                         strip.white = FALSE, check.names = FALSE)
  tab <- as.matrix(tab)
  tab[] <- sub("\r$", "", tab)           # template mixes LF and CRLF endings
  dimnames(tab) <- NULL
  if (nrow(tab) < .CP_N_PARAMS + 1L) {
    stop("Parameter file has ", nrow(tab), " rows; expected a header plus ",
         .CP_N_PARAMS, " parameter rows: ", file, call. = FALSE)
  }
  tab
}

#' Path to the bundled parameter template
#' @noRd
char_params_template_path <- function() {
  f <- system.file("templates", "template_charParams.csv",
                   package = "CharAnalysis")
  if (!nzchar(f)) stop("Bundled template_charParams.csv not found.", call. = FALSE)
  f
}

#' Convert char_parameters() output to the 25-value parameter vector
#'
#' Inverse of the unpacking in char_parameters(): zoneDiv is padded to 8
#' values with the -9999 sentinel.
#' @param p List returned by char_parameters().
#' @return Numeric vector of length 25, in file row order.
#' @noRd
char_params_values <- function(p) {
  zd <- p$pretreatment$zoneDiv
  if (length(zd) > 8L) stop("zoneDiv has more than 8 values.", call. = FALSE)
  v <- c(zd, rep(-9999, 8L - length(zd)),
         p$pretreatment$yrInterp, p$pretreatment$transform,
         p$smoothing$method, p$smoothing$yr,
         p$peak_analysis$cPeak, p$peak_analysis$threshType,
         p$peak_analysis$threshMethod, p$peak_analysis$threshValues,
         p$peak_analysis$minCountP, p$peak_analysis$peakFrequ,
         p$peak_analysis$bkgSens,
         p$results$saveFigures, p$results$save, p$results$allFigures)
  if (length(v) != .CP_N_PARAMS) {
    stop("Expected ", .CP_N_PARAMS, " parameter values; got ", length(v), ".",
         call. = FALSE)
  }
  as.numeric(v)
}

# -- internal utilities --------------------------------------------------------

.cp_check_values <- function(values) {
  if (!is.numeric(values) || length(values) != .CP_N_PARAMS) {
    stop("`values` must be a numeric vector of length ", .CP_N_PARAMS, ".",
         call. = FALSE)
  }
  if (any(!is.finite(values))) {
    stop("`values` contains NA or non-finite entries at position(s): ",
         paste(which(!is.finite(values)), collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

.cp_check_name <- function(file) {
  if (!grepl("_charParams\\.csv$", basename(file))) {
    stop("Parameter file must be named <site>_charParams.csv: ",
         basename(file), call. = FALSE)
  }
  invisible(TRUE)
}

# Numbers written plainly (no scientific notation, no trailing zeros), so the
# MATLAB %n field and as.numeric() both parse them.
.cp_fmt_num <- function(x) {
  out <- formatC(x, digits = 15L, format = "fg", flag = "")
  trimws(out)
}

# Minimal CSV quoting: quote only fields containing a comma, quote, or line
# break. Column 3 values are numeric strings and therefore never quoted.
.cp_csv_field <- function(x) {
  needs <- grepl("[\",\r\n]", x)
  x[needs] <- paste0("\"", gsub("\"", "\"\"", x[needs]), "\"")
  x
}

# Write via a temporary file in the same folder, then swap it in, so an
# interrupted save cannot leave a partial parameter file.
.cp_write_atomic <- function(tab, file) {
  if (anyNA(iconv(tab, from = "", to = "ASCII"))) {
    stop("Parameter file content must be plain ASCII (required by the MATLAB ",
         "reader).", call. = FALSE)
  }
  lines <- apply(tab, 1L, function(r) paste(.cp_csv_field(r), collapse = ","))
  tmp <- tempfile(pattern = ".charParams_", tmpdir = dirname(file),
                  fileext = ".tmp")
  con <- file(tmp, open = "wb")
  writeLines(lines, con, sep = "\n")
  close(con)
  if (!identical(readLines(tmp, warn = FALSE), unname(lines))) {
    unlink(tmp)
    stop("Verification of the temporary parameter file failed; nothing was ",
         "written.", call. = FALSE)
  }
  if (!file.rename(tmp, file)) {
    ok <- file.copy(tmp, file, overwrite = TRUE)
    unlink(tmp)
    if (!ok) stop("Could not write parameter file: ", file, call. = FALSE)
  }
  invisible(file)
}

# -- exported-to-package internals ---------------------------------------------

#' Write a new parameter file from the template
#'
#' Refuses to overwrite an existing file; use char_params_revise() instead.
#' @param values Numeric vector of 25 parameter values in file row order.
#' @param file Destination path, named <site>_charParams.csv.
#' @param template Template path (default: bundled template).
#' @return The file path, invisibly.
#' @noRd
char_params_new <- function(values, file, template = NULL) {
  .cp_check_name(file)
  .cp_check_values(values)
  if (file.exists(file)) {
    stop("A parameter file already exists and will not be overwritten: ",
         file, "\nRevise it instead (char_params_revise()).", call. = FALSE)
  }
  if (is.null(template)) template <- char_params_template_path()
  tmpl <- char_params_read_table(template)
  rows <- seq_len(.CP_N_PARAMS + 1L)
  tab  <- matrix("", nrow = length(rows), ncol = .CP_LOG_COL)
  tab[, 1:5] <- tmpl[rows, 1:5]
  # Keep the template's "Parameter Log" label in the header of column 7;
  # log columns themselves start empty.
  if (ncol(tmpl) >= .CP_LOG_COL) tab[1L, .CP_LOG_COL] <- tmpl[1L, .CP_LOG_COL]
  tab[-1L, 3L] <- .cp_fmt_num(values)
  .cp_write_atomic(tab, file)
}

#' Revise an existing parameter file, logging the previous values
#'
#' If any value differs, the previous column-3 values move to a new log column
#' (column 7, headed with a time stamp and package version) and the new values
#' are written to column 3. Existing log columns, including hand-entered ones,
#' shift right unchanged. Log columns with no content in the parameter rows
#' (e.g., the template's empty labelled column) are dropped. Rows below the
#' parameter block, if any, are kept.
#' @param file Existing <site>_charParams.csv.
#' @param values Numeric vector of 25 new parameter values.
#' @param time Time stamp for the log header.
#' @param version Package version recorded in the log header.
#' @return TRUE (invisibly) if the file was revised; FALSE if nothing changed.
#' @noRd
char_params_revise <- function(file, values, time = Sys.time(),
                               version = NULL) {
  .cp_check_name(file)
  .cp_check_values(values)
  tab <- char_params_read_table(file)
  prow <- 2:(.CP_N_PARAMS + 1L)

  old_str <- tab[prow, 3L]
  old     <- suppressWarnings(as.numeric(old_str))
  if (!anyNA(old) && all(old == values)) {
    message("No parameter values changed; ", basename(file), " was not rewritten.")
    return(invisible(FALSE))
  }

  if (ncol(tab) < 6L) tab <- cbind(tab, matrix("", nrow(tab), 6L - ncol(tab)))
  log_idx <- if (ncol(tab) >= .CP_LOG_COL) .CP_LOG_COL:ncol(tab) else integer(0)
  keep    <- log_idx[vapply(log_idx, function(j) any(nzchar(trimws(tab[prow, j]))),
                            logical(1))]
  if (length(keep) >= .CP_MAX_LOG) {
    stop(basename(file), " already holds ", length(keep), " logged parameter ",
         "sets, the maximum (", .CP_MAX_LOG, ") that keeps the file readable ",
         "by the MATLAB version of CharAnalysis.\nSave a copy of the file to ",
         "keep its history, then delete older log columns (rightmost) before ",
         "revising.", call. = FALSE)
  }

  if (is.null(version)) {
    version <- tryCatch(as.character(utils::packageVersion("CharAnalysis")),
                        error = function(e) "unknown")
  }
  new_col <- rep("", nrow(tab))
  new_col[1L]   <- sprintf("Replaced %s, CharAnalysis %s",
                           format(time, "%Y-%m-%d %H:%M"), version)
  new_col[prow] <- old_str

  out <- cbind(tab[, 1:6, drop = FALSE], new_col, tab[, keep, drop = FALSE])
  dimnames(out) <- NULL
  out[prow, 3L] <- .cp_fmt_num(values)
  .cp_write_atomic(out, file)
  invisible(TRUE)
}
