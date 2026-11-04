# Tests for the parameter-file writer (ROADMAP_ParamsGUI.md, Phase 1)
#
# 1. Round trip: for each validation site, read -> write a new file -> re-read
#    gives identical parameters and identical CharAnalysis() results.
# 2. A new file is never written over an existing one.
# 3. Revision: previous values move to a time-stamped log column; unchanged
#    values do not rewrite the file; hand-entered log columns are preserved;
#    the 10-log limit is enforced; revised files still read correctly.
# 4. Files not named <site>_charParams.csv are rejected.
#
# Runs against the installed package. Set CHARANALYSIS_VAL_DIR and
# CHARANALYSIS_TEMPLATE to use source files instead.
# AI-assisted development (Claude/Anthropic).

suppressPackageStartupMessages(library(CharAnalysis))

fn <- function(n) {
  if (exists(n, envir = globalenv(), mode = "function")) get(n, envir = globalenv())
  else utils::getFromNamespace(n, "CharAnalysis")
}
char_params_values     <- fn("char_params_values")
char_params_new        <- fn("char_params_new")
char_params_revise     <- fn("char_params_revise")
char_params_read_table <- fn("char_params_read_table")
if (exists("char_parameters", envir = globalenv())) {
  char_parameters <- get("char_parameters", envir = globalenv())
  CharAnalysis    <- get("CharAnalysis", envir = globalenv())
}

val_dir  <- Sys.getenv("CHARANALYSIS_VAL_DIR",
                       system.file("validation", package = "CharAnalysis"))
template <- Sys.getenv("CHARANALYSIS_TEMPLATE", "")
template <- if (nzchar(template)) template else NULL

expect_error <- function(expr, pattern) {
  msg <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  if (is.null(msg)) stop("Expected an error matching '", pattern, "'; none raised.")
  if (!grepl(pattern, msg)) stop("Error message did not match '", pattern, "':\n", msg)
  invisible(TRUE)
}
results <- character(0)
pass <- function(label) { results <<- c(results, label); cat("PASS:", label, "\n") }

work <- file.path(tempdir(), "params_writer_test")
unlink(work, recursive = TRUE); dir.create(work)

# -- 1. Round trip -------------------------------------------------------------
sites <- c("CO", "CO_compensated", "CH10", "SI17", "RA07")
for (s in sites) {
  orig <- file.path(val_dir, paste0(s, "_charParams.csv"))
  p1   <- char_parameters(orig)
  v    <- char_params_values(p1)
  raw  <- as.numeric(char_params_read_table(orig)[2:26, 3])
  stopifnot(identical(v, raw))

  d <- file.path(work, s); dir.create(d)
  file.copy(file.path(val_dir, paste0(s, "_charData.csv")), d)
  new <- file.path(d, paste0(s, "_charParams.csv"))
  char_params_new(v, new, template = template)
  p2 <- char_parameters(new)
  stopifnot(identical(p1, p2))

  set.seed(1); o1 <- suppressMessages(CharAnalysis(orig, plots = FALSE, save = FALSE))
  set.seed(1); o2 <- suppressMessages(CharAnalysis(new,  plots = FALSE, save = FALSE))
  stopifnot(identical(o1$char_results, o2$char_results))
  pass(paste0("round trip, identical parameters and results: ", s))
}

# -- 2. No overwrite -----------------------------------------------------------
co_new <- file.path(work, "CO", "CO_charParams.csv")
before <- readLines(co_new)
expect_error(char_params_new(char_params_values(char_parameters(co_new)), co_new,
                             template = template), "will not be overwritten")
stopifnot(identical(before, readLines(co_new)))
pass("new file refuses to overwrite an existing file")

# -- 3. Revision ---------------------------------------------------------------
v0 <- char_params_values(char_parameters(co_new))

# 3a. Unchanged values: no rewrite
stopifnot(isFALSE(suppressMessages(char_params_revise(co_new, v0, version = "test"))))
stopifnot(identical(before, readLines(co_new)))
pass("unchanged values do not rewrite the file")

# 3b. One change: old values logged in column 7, new values in column 3
v1 <- v0; v1[9] <- 20                          # yrInterpolate 15 -> 20
t1 <- as.POSIXct("2026-11-03 14:32", tz = "UTC")
stopifnot(isTRUE(char_params_revise(co_new, v1, time = t1, version = "test")))
tab <- char_params_read_table(co_new)
stopifnot(ncol(tab) == 7L,
          tab[1, 7] == "Replaced 2026-11-03 14:32, CharAnalysis test",
          identical(as.numeric(tab[2:26, 7]), v0),
          identical(as.numeric(tab[2:26, 3]), v1),
          char_parameters(co_new)$pretreatment$yrInterp == 20)
pass("revision logs previous values with time stamp; new values read correctly")

# 3c. Hand-entered log column is preserved and shifts right
tab <- char_params_read_table(co_new)
hand <- c("my notes, run A", rep("", 25)); hand[c(2, 10)] <- c("-51", "15")
writeLines(apply(cbind(tab, hand), 1, function(r) {
  q <- grepl("[\",]", r); r[q] <- paste0("\"", gsub("\"", "\"\"", r[q]), "\""); paste(r, collapse = ",")
}), co_new)
v2 <- v1; v2[12] <- 600                         # smoothing window 500 -> 600
char_params_revise(co_new, v2, time = t1 + 60, version = "test")
tab <- char_params_read_table(co_new)
stopifnot(ncol(tab) == 9L,
          identical(as.numeric(tab[2:26, 7]), v1),
          identical(as.numeric(tab[2:26, 8]), v0),
          tab[1, 9] == "my notes, run A", tab[2, 9] == "-51", tab[10, 9] == "15",
          identical(char_params_values(char_parameters(co_new)), v2))
pass("hand-entered log column preserved; logs ordered newest first")

# 3d. Limit of 10 log columns (16 columns in total)
v <- v2
while (ncol(char_params_read_table(co_new)) < 16L) {
  v[12] <- v[12] + 10
  char_params_revise(co_new, v, version = "test")
}
stopifnot(ncol(char_params_read_table(co_new)) == 16L,
          max(utils::count.fields(co_new, sep = ",", quote = "\"")) == 16L)
before <- readLines(co_new)
v[12] <- v[12] + 10
expect_error(char_params_revise(co_new, v, version = "test"), "maximum \\(10\\)")
stopifnot(identical(before, readLines(co_new)))
set.seed(1); o3 <- suppressMessages(CharAnalysis(co_new, plots = FALSE, save = FALSE))
stopifnot(nrow(o3$char_results) > 0)
pass("10-log limit enforced; a full 16-column file still runs")

# -- 4. Naming -----------------------------------------------------------------
expect_error(char_params_new(v0, file.path(work, "CO_params.csv"), template = template),
             "_charParams.csv")
pass("files not named <site>_charParams.csv are rejected")

cat("\nAll", length(results), "parameter-writer tests passed.\n")
