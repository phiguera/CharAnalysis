#' Write the CharAnalysis results matrix to a CSV file
#'
#' Writes the \eqn{N \times 33} \code{char_results} matrix (assembled by
#' [char_post_process()]) to a CSV file whose column headers and numeric
#' format match the MATLAB CharAnalysis v2.0 \code{*_charResults.csv} output
#' exactly.
#'
#' @param char_results Numeric matrix (\eqn{N \times 33}) returned in
#'   \code{out$char_results} by [CharAnalysis()].
#' @param site         Character string: site name, used to build the
#'   output filename (\code{<site>_charResults.csv}).
#' @param out_dir      Directory in which to create the file.  Defaults to
#'   the current working directory.  Created if it does not exist.
#' @param digits       Number of significant digits for numeric output.
#'   Default 7, matching MATLAB's \code{fprintf} default precision.
#'   Use \code{NULL} for R's full double precision (15 digits).
#'
#' @return Invisibly returns the full path to the file written.
#'
#' @details
#'   ## Column order and headers
#'   Columns follow the MATLAB \code{charResults} layout:
#'   \enumerate{
#'     \item cm Top_i (cm)
#'     \item age Top_i (yr BP)
#'     \item char Count_i (#)
#'     \item char Vol_i (cm3)
#'     \item char Con_i (# cm-3)
#'     \item char Acc_i (# cm-2 yr-1)
#'     \item charBkg (# cm-2 yr-1)
#'     \item char Peak (# cm-2 yr-1)
#'     \item thresh 1 (# cm-2 yr-1)
#'     \item thresh 2 (# cm-2 yr-1)
#'     \item thresh 3 (# cm-2 yr-1)
#'     \item thresh FinalPos (# cm-2 yr-1)
#'     \item thresh FinalNeg (# cm-2 yr-1)
#'     \item SNI (index)
#'     \item thresh GOF (p-val)
#'     \item peaks 1
#'     \item peaks 2
#'     \item peaks 3
#'     \item peaks Final
#'     \item peaks Insig.
#'     \item peak Mag (# cm-2 peak-1)
#'     \item smPeak Frequ (peaks 1ka-1)
#'     \item smFRIs (yr fire-1)
#'     \item nFRIs (#)
#'     \item mFRI (yr fire-1)
#'     \item mFRI_uCI (yr fire-1)
#'     \item mFRI_lCI (yr fire-1)
#'     \item WBLb (yr)
#'     \item WBLb_uCI (yr)
#'     \item WBLb_lCI (yr)
#'     \item WBLc (unitless)
#'     \item WBLc_uCI (unitless)
#'     \item WBLc_lCI (unitless)
#'   }
#'
#'   ## NA / empty handling
#'   \code{NA} values are written as empty fields (no quotes), matching
#'   MATLAB's blank-cell convention.  This applies to:
#'   \itemize{
#'     \item \code{smFRIs} rows beyond the smoothed-FRI coverage window;
#'     \item zone-statistics columns (24-33) for rows beyond the last zone;
#'     \item any column not computed for a given run configuration.
#'   }
#'
#'   ## CI column convention
#'   MATLAB stores bootstrap CIs as \code{[quantile(2.5\%), quantile(97.5\%)]}
#'   in the columns labelled \code{uCI} / \code{lCI} respectively (i.e.
#'   \code{uCI} = lower bound, \code{lCI} = upper bound -- MATLAB's own
#'   labelling is inverted).  The R output follows the same convention so that
#'   column indices are identical to the MATLAB reference file.
#'
#' @seealso [char_post_process()], [CharAnalysis()]
#'
#' @examples
#' \dontrun{
#'   out <- CharAnalysis("CO_charParams.csv")
#'   char_write_results(out$char_results, out$site,
#'                      out_dir = "Results")
#' }
CharWriteResults <- function(CharResults,
                                site,
                                out_dir = ".",
                                digits  = 7L) {

  # ---- Validate inputs -------------------------------------------------------
  if (!is.matrix(char_results) || ncol(char_results) != 33L)
    stop("char_results must be a numeric matrix with exactly 33 columns.")
  if (!is.character(site) || length(site) != 1L || nchar(site) == 0L)
    stop("site must be a non-empty character string.")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    message("char_write_results: created output directory '", out_dir, "'")
  }

  # ---- Column headers (matching MATLAB charResults exactly) ------------------
  col_headers <- c(
    "cm Top_i (cm)",
    "age Top_i (yr BP)",
    "char Count_i (#)",
    "char Vol_i (cm3)",
    "char Con_i (# cm-3)",
    "char Acc_i (# cm-2 yr-1)",
    "charBkg (# cm-2 yr-1)",
    "char Peak (# cm-2 yr-1)",
    "thresh 1 (# cm-2 yr-1)",
    "thresh 2 (# cm-2 yr-1)",
    "thresh 3 (# cm-2 yr-1)",
    "thresh FinalPos (# cm-2 yr-1)",
    "thresh FinalNeg (# cm-2 yr-1)",
    "SNI (index)",
    "thresh GOF (p-val)",
    "peaks 1",
    "peaks 2",
    "peaks 3",
    "peaks Final",
    "peaks Insig.",
    "peak Mag (# cm-2 peak-1)",
    "smPeak Frequ (peaks 1ka-1)",
    "smFRIs (yr fire-1)",
    "nFRIs (#)",
    "mFRI (yr fire-1)",
    "mFRI_uCI (yr fire-1)",
    "mFRI_lCI (yr fire-1)",
    "WBLb (yr)",
    "WBLb_uCI (yr)",
    "WBLb_lCI (yr)",
    "WBLc (unitless)",
    "WBLc_uCI (unitless)",
    "WBLc_lCI (unitless)"
  )

  # ---- Format numeric matrix -------------------------------------------------
  # Convert each column to character, respecting the digits argument.
  # NA  -> "" (blank cell, matching MATLAB blank-cell convention).
  # Integer-valued columns (peaks, peakInsig) are written without decimals.

  integer_cols <- 16:20   # peaks 1-Final, peakInsig -- always 0/1 integers

  fmt_val <- function(x, is_int = FALSE) {
    if (is.na(x))   return("")
    if (is_int)     return(as.character(as.integer(x)))
    if (is.null(digits)) return(as.character(x))
    # Match MATLAB's %g-style: use formatC with "g" format
    formatC(x, digits = digits, format = "g", flag = "-")
  }

  N   <- nrow(char_results)
  out <- matrix("", nrow = N, ncol = 33L)

  for (j in seq_len(33L)) {
    is_int <- j %in% integer_cols
    out[, j] <- vapply(char_results[, j], fmt_val,
                       FUN.VALUE = character(1L),
                       is_int    = is_int)
  }

  # ---- Assemble data frame and write -----------------------------------------
  df <- as.data.frame(out, stringsAsFactors = FALSE)
  names(df) <- col_headers

  out_path <- file.path(out_dir,
                         paste0(site, "_charResults.csv"))

  # write.csv adds row names and quotes by default; use write.table for control
  utils::write.table(df,
                     file      = out_path,
                     sep       = ",",
                     col.names = TRUE,
                     row.names = FALSE,
                     quote     = FALSE,
                     na        = "")

  message("char_write_results: wrote '", out_path, "'")
  invisible(out_path)
}
