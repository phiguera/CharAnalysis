# char_extract_bacon_chronologies.R
# ---------------------------------------------------------------------
# Development script for a function destined for the CharAnalysis
# package (R/ directory). Extracts n_iter chronologies from a
# completed rbacon MCMC posterior at user-supplied depths.
#
# Usage (in the same R session as Bacon()):
#   source("CharAnalysis_2_0_R/tests/char_extract_bacon_chronologies.R")
#   chron <- char_extract_bacon_chronologies(
#               depths   = chardata$cmTop,
#               out_file = "CharAnalysis_2_0_R/tests/SI17_MCAgeDepth_1000_chronologies.csv"
#            )
# ---------------------------------------------------------------------


#' Extract a chronology ensemble from a completed rbacon run
#'
#' Samples \code{n_iter} age-depth chronologies from the MCMC posterior
#' of a \pkg{rbacon} model and returns them as a data frame suitable for
#' use with \code{char_run_ensemble()}.
#'
#' @param depths Numeric vector of depths (cm) at which to extract ages.
#'   Typically \code{charData$cmTop} for the target site.
#' @param n_iter Integer. Number of chronologies to extract. Default 1000.
#' @param set The rbacon \code{info} object. Defaults to \code{get('info')},
#'   which is available in the session immediately after \code{Bacon()} runs.
#' @param seed Integer random seed for reproducible sampling. Default 42.
#' @param out_file Optional character path for a CSV output file. If
#'   \code{NULL} (default), no file is written.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A data frame with columns \code{Sample_cm}, \code{CalAge_1}, …,
#'   \code{CalAge_n_iter} (cal yr BP). Returned invisibly when
#'   \code{out_file} is supplied, visibly otherwise.
#'
#' @details
#' Each call to \code{Bacon.Age.d()} returns a numeric vector of length
#' equal to the number of stored MCMC iterations, ordered by iteration.
#' This function samples \code{n_iter} row indices once and applies them
#' consistently across all depths, ensuring that ages at different depths
#' within a single chronology come from the same MCMC iteration and are
#' therefore properly correlated.
#'
#' @seealso \code{\link{char_run_ensemble}}
#'
#' @examples
#' \dontrun{
#'   library(rbacon)
#'   Bacon("MySite", coredir = "~/Bacon_runs/", ask = FALSE)
#'
#'   chardata <- read.csv("MySite_charData.csv", check.names = FALSE)
#'   chron <- char_extract_bacon_chronologies(
#'     depths   = chardata$cmTop,
#'     out_file = "MySite_MCAgeDepth_1000_chronologies.csv"
#'   )
#' }

char_extract_bacon_chronologies <- function(depths,
                                            n_iter   = 1000L,
                                            set      = get('info'),
                                            seed     = 42L,
                                            out_file = NULL,
                                            verbose  = TRUE) {

  # --- input checks ---------------------------------------------------
  if (!is.numeric(depths) || length(depths) < 2L)
    stop("'depths' must be a numeric vector of length >= 2.")

  if (!is.numeric(n_iter) || length(n_iter) != 1L || n_iter < 1L)
    stop("'n_iter' must be a positive integer.")
  n_iter <- as.integer(n_iter)

  if (!requireNamespace("rbacon", quietly = TRUE))
    stop("Package 'rbacon' is required. Install with install.packages('rbacon').")

  # --- determine available MCMC iterations ----------------------------
  n_mcmc <- length(rbacon::Bacon.Age.d(depths[1L], set = set, BCAD = FALSE))

  if (n_iter > n_mcmc)
    stop("Requested n_iter (", n_iter, ") exceeds available MCMC iterations (",
         n_mcmc, "). Re-run Bacon with a larger 'ssize'.")

  if (verbose)
    cat("Available MCMC iterations:", n_mcmc, "\n",
        "Sampling:", n_iter, "chronologies at",
        length(depths), "depths (", min(depths), "-", max(depths), "cm)\n\n")

  # --- sample indices once --------------------------------------------
  set.seed(seed)
  idx <- sample(n_mcmc, n_iter, replace = FALSE)

  # --- extract --------------------------------------------------------
  t0 <- proc.time()
  age_matrix <- matrix(NA_real_, nrow = length(depths), ncol = n_iter)

  for (i in seq_along(depths)) {
    all_ages        <- rbacon::Bacon.Age.d(depths[i], set = set, BCAD = FALSE)
    age_matrix[i, ] <- all_ages[idx]

    if (verbose && i %% 100L == 0L) {
      elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
      cat("  depth", i, "of", length(depths), "—", elapsed, "min\n")
    }
  }

  if (verbose) {
    elapsed_total <- round((proc.time() - t0)["elapsed"] / 60, 1)
    cat("Done. Total time:", elapsed_total, "min\n\n")
  }

  # --- assemble output data frame -------------------------------------
  chron_df        <- as.data.frame(age_matrix)
  names(chron_df) <- paste0("CalAge_", seq_len(n_iter))
  chron_df        <- cbind(Sample_cm = depths, chron_df)

  # --- write CSV if requested -----------------------------------------
  if (!is.null(out_file)) {
    write.csv(chron_df, out_file, row.names = FALSE)
    if (verbose)
      cat("Saved:", out_file, "\n",
          "Dimensions:", nrow(chron_df), "depths x", n_iter, "chronologies\n")
    return(invisible(chron_df))
  }

  chron_df
}
