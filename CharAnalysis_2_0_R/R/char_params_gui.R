# Parameter interface (ROADMAP_ParamsGUI.md, Phases 2-3; prototype V3)
#
# A Shiny gadget that guides the user through every parameter in a
# <site>_charParams.csv file, then writes the file with char_params_new() or
# char_params_revise() (char_params_write.R). Launched from CharAnalysis()
# when no parameter file is given, or directly with char_params_gui().
#
# Layout: tabs across the top (Data, Pretreatment, Smoothing, Peaks, Results,
# Review); each tab is a table with one row per parameter (name and plain
# definition | setting | description), an info icon with more detail, notes
# shown in the row they apply to, and Cancel / Save buttons at the bottom.
#
# shiny and miniUI are suggested packages; the user is asked before they are
# installed.
#
# AI-assisted development (Claude/Anthropic).

#' Create or revise a parameter file with a graphical interface
#'
#' Opens a form (in the RStudio Viewer or a web browser) containing the
#' parameters in a \code{<site>_charParams.csv} file, with a summary of the
#' charcoal record, plain-language definitions, and advisory notes. On saving,
#' the parameter file is written next to the data file. When revising an
#' existing file, the previous values are kept in the file's Parameter Log
#' columns.
#'
#' @param file Path to either a data file named \code{<site>_charData.csv}
#'   (create a new parameter file) or an existing parameter file named
#'   \code{<site>_charParams.csv} (revise it). If \code{NULL}, you are asked
#'   which to do and a file picker opens.
#'
#' @return The path to the saved parameter file, invisibly, or \code{NULL} if
#'   the interface was cancelled. Run the analysis with
#'   \code{CharAnalysis(<path>)}.
#'
#' @details Requires the suggested packages \pkg{shiny} and \pkg{miniUI}.
#'   Data files must follow the \code{<site>_charData.csv} naming convention,
#'   because \code{CharAnalysis()} finds the data file from the parameter file
#'   name. An existing parameter file is never changed without its previous
#'   values being logged; a file can hold up to 10 logged parameter sets.
#'   The \code{allFigures} setting is not shown: it is used only by the MATLAB
#'   version, and its existing value (or 1 for new files) is kept.
#'   \code{CharAnalysis()} called with no file offers the same interface and
#'   runs the analysis after saving.
#'
#' @seealso [CharAnalysis()], [char_parameters()]
#' @export
char_params_gui <- function(file = NULL) {
  .cp_gui(file = file, mode = NULL, button = "Save")
}

# Shared entry point for char_params_gui() and CharAnalysis().
# mode: NULL (ask), "new" (choose a data file), or "revise" (choose a params file).
.cp_gui <- function(file = NULL, mode = NULL, button = "Save") {
  if (!interactive()) {
    stop("The parameter interface requires an interactive R session.", call. = FALSE)
  }
  .cp_require_gui_pkgs()

  if (is.null(file)) {
    if (is.null(mode)) {
      choice <- utils::menu(
        c("Create a new parameter file from a data file (<site>_charData.csv)",
          "Revise an existing parameter file (<site>_charParams.csv)"),
        title = "Parameter interface: what would you like to do?")
      if (choice == 0L) return(invisible(NULL))
      mode <- c("new", "revise")[choice]
    }
    message(if (mode == "new") "Select the data file (<site>_charData.csv)..."
            else "Select the parameter file (<site>_charParams.csv)...")
    file <- .cp_choose_file()
  }

  file <- normalizePath(file, winslash = "/", mustWork = TRUE)
  base <- basename(file)

  if (grepl("_charData\\.csv$", base)) {
    mode        <- "new"
    site        <- sub("_charData\\.csv$", "", base)
    data_file   <- file
    params_file <- file.path(dirname(file), paste0(site, "_charParams.csv"))
    if (file.exists(params_file)) {
      ans <- utils::menu(c("Yes, revise the existing file", "No, cancel"),
                         title = paste0("A parameter file already exists for ",
                                        "this site:\n  ", params_file,
                                        "\nRevise it instead?"))
      if (ans != 1L) return(invisible(NULL))
      mode <- "revise"
    }
  } else if (grepl("_charParams\\.csv$", base)) {
    mode        <- "revise"
    site        <- sub("_charParams\\.csv$", "", base)
    params_file <- file
    data_file   <- file.path(dirname(file), paste0(site, "_charData.csv"))
    if (!file.exists(data_file)) {
      stop("The companion data file was not found:\n  ", data_file,
           "\nCharAnalysis requires the data file and parameter file to sit ",
           "in the same folder, named <site>_charData.csv and ",
           "<site>_charParams.csv.", call. = FALSE)
    }
  } else {
    stop(.cp_naming_message(base), call. = FALSE)
  }

  char_data <- as.matrix(utils::read.csv(data_file, header = TRUE,
                                         stringsAsFactors = FALSE))
  storage.mode(char_data) <- "double"
  summ <- .cp_data_summary(char_data)

  if (mode == "new") {
    init <- as.numeric(char_params_read_table(char_params_template_path())[2:26, 3])
    init[1:8] <- c(summ$age_top, summ$age_bot, rep(-9999, 6))  # start/end only
    init[9]   <- summ$res_suggest
    old       <- NULL
    n_logs    <- 0L
  } else {
    tab    <- char_params_read_table(params_file)
    old    <- suppressWarnings(as.numeric(tab[2:26, 3]))
    init   <- old
    n_logs <- if (ncol(tab) >= .CP_LOG_COL) {
      sum(vapply(.CP_LOG_COL:ncol(tab),
                 function(j) any(nzchar(trimws(tab[2:26, j]))), logical(1)))
    } else 0L
  }

  path <- .cp_run_gadget(mode, site, data_file, params_file, char_data, summ,
                         init, old, n_logs, button)
  if (is.null(path)) {
    message("Parameter interface cancelled; no file written.")
    return(invisible(NULL))
  }
  message("Parameter file: ", path)
  message("To run (or re-run) this analysis from the console:\n  CharAnalysis(\"",
          path, "\")")
  invisible(path)
}

# -- helpers -------------------------------------------------------------------

.cp_require_gui_pkgs <- function() {
  pkgs <- c("shiny", "miniUI")
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need) == 0L) return(invisible(TRUE))
  ok <- isTRUE(utils::askYesNo(paste0(
    "The parameter interface needs the package(s): ", paste(need, collapse = ", "),
    ". Install now?")))
  if (!ok) {
    stop("Install ", paste(need, collapse = " and "), " to use the parameter ",
         "interface: install.packages(c(", paste0("\"", need, "\"", collapse = ", "),
         "))", call. = FALSE)
  }
  utils::install.packages(need)
  invisible(TRUE)
}

.cp_choose_file <- function() {
  f <- tryCatch(file.choose(), error = function(e) readline("Enter file path: "))
  if (!nzchar(f)) stop("No file selected.", call. = FALSE)
  f
}

.cp_naming_message <- function(base) {
  paste0("'", base, "' does not follow the CharAnalysis naming convention.\n",
         "Data files must be named <site>_charData.csv (e.g., CO_charData.csv), ",
         "and parameter files <site>_charParams.csv. CharAnalysis finds the data ",
         "file from the parameter file name, so this convention is required. ",
         "Rename the file and try again.")
}

# Summary of the charcoal record used for display, suggestions, and notes.
.cp_data_summary <- function(char_data) {
  n   <- nrow(char_data)
  res <- diff(char_data[, 3L])
  list(
    n           = n,
    cm_top      = unname(char_data[1L, 1L]), cm_bot  = unname(char_data[n, 2L]),
    age_top     = unname(char_data[1L, 3L]), age_bot = unname(char_data[n, 4L]),
    res_median  = stats::median(res),
    res_min     = min(res), res_max = max(res),
    res_suggest = round(stats::median(res)),     # as CharAnalysis does for 0
    n_missing   = sum(char_data[, 5L] <= 0, na.rm = TRUE),
    n_gaps      = sum((char_data[-1L, 1L] - char_data[-n, 2L]) != 0, na.rm = TRUE)
  )
}

# Row key for each of the 25 parameter positions (used to place messages).
.cp_pos_key <- c(rep("zoneDiv", 8), "yrInterp", "transform", "method", "smoYr",
                 "cPeak", "threshType", "threshMethod", rep("threshValues", 4),
                 "minCountP", "peakFrequ", "bkgSens", "saveFigures", "saveData",
                 "allFigures")

.cp_param_labels <- c(
  "Zone start", paste("Zone value", 2:8),
  "yrInterpolate", "transform", "Smoothing method", "Smoothing window (yr)",
  "cPeak", "threshType", "threshMethod",
  "Peaks 1 threshold", "Peaks 2 threshold", "Peaks 3 threshold",
  "Peaks Final threshold", "minCountP", "peakFrequ", "C_background sensitivity",
  "saveFigures", "saveData", "allFigures")

.cp_key_tab <- c(data = "Data", zoneDiv = "Pretreatment", yrInterp = "Pretreatment",
                 transform = "Pretreatment", method = "Smoothing", smoYr = "Smoothing",
                 cPeak = "Peaks", threshType = "Peaks", threshMethod = "Peaks",
                 threshValues = "Peaks", minCountP = "Peaks", peakFrequ = "Results",
                 bkgSens = "Results", saveFigures = "Results", saveData = "Results",
                 allFigures = "Results")

# Same unpacking as char_parameters(); used to run char_validate_params().
.cp_unpack <- function(v) {
  zd <- v[1:8]; zd <- zd[!is.na(zd) & zd != -9999]
  list(pretreatment  = list(zoneDiv = zd, yrInterp = v[9], transform = v[10]),
       smoothing     = list(method = v[11], yr = v[12]),
       peak_analysis = list(cPeak = v[13], threshType = v[14],
                            threshMethod = v[15], threshValues = v[16:19],
                            minCountP = v[20], peakFrequ = v[21], bkgSens = v[22]),
       results       = list(saveFigures = v[23], save = v[24], allFigures = v[25]))
}

# Which parameter row a char_validate_params() message belongs to.
.cp_key_of <- function(msg) {
  rules <- c("smoothing window" = "smoYr", "smoothing method" = "method",
             "zoneDiv" = "zoneDiv", "yrInterp" = "yrInterp",
             "locally defined threshold" = "threshMethod",
             "threshValues" = "threshValues", "threshType" = "threshType",
             "threshMethod" = "threshMethod", "cPeak" = "cPeak",
             "minCountP" = "minCountP", "peakFrequ" = "peakFrequ")
  for (p in names(rules)) if (grepl(p, msg, ignore.case = TRUE)) return(rules[[p]])
  "data"
}

# Comma-separated zone boundaries -> numeric vector (or an error message).
.cp_parse_bounds <- function(txt) {
  txt <- trimws(if (is.null(txt)) "" else txt)
  if (!nzchar(txt)) return(list(values = numeric(0), error = NULL))
  parts <- trimws(strsplit(txt, "[,;]")[[1]])
  parts <- parts[nzchar(parts)]
  num   <- suppressWarnings(as.numeric(parts))
  if (anyNA(num)) {
    return(list(values = numeric(0), error = paste0(
      "Zone boundaries must be numbers separated by commas; could not read: ",
      paste(parts[is.na(num)], collapse = ", "), ".")))
  }
  if (length(num) > 6L) {
    return(list(values = numeric(0), error = paste0(
      "At most six zone boundaries can be entered (found ", length(num), ").")))
  }
  list(values = num, error = NULL)
}

# Errors block saving; notes advise. Each is a data frame (key, text).
.cp_check <- function(v, summ, char_data, bounds_error = NULL) {
  err  <- data.frame(key = character(0), text = character(0))
  warn <- err
  add  <- function(df, k, t) rbind(df, data.frame(key = k, text = t))

  if (summ$n_missing > 0L) {
    warn <- add(warn, "data", paste0(summ$n_missing, " sample(s) have missing ",
      "values (volume <= 0); CharAnalysis will interpolate across them."))
  }
  if (summ$n_gaps > 0L) {
    warn <- add(warn, "data", paste0(summ$n_gaps, " break(s) in sample depths ",
      "(cmTop not equal to the previous cmBot); resampling will occur across ",
      "them. Consider filling missing rows with -999 to interpolate instead."))
  }
  if (!is.null(bounds_error)) err <- add(err, "zoneDiv", bounds_error)

  if (length(v) != 25L || anyNA(v)) {
    for (i in which(is.na(v))) {
      err <- add(err, .cp_pos_key[i], paste0(.cp_param_labels[i], " is empty."))
    }
    return(list(errors = err, warnings = warn))
  }

  p  <- .cp_unpack(v)
  zd <- p$pretreatment$zoneDiv
  if (length(zd) > 2L) {
    inner <- zd[2:(length(zd) - 1L)]
    if (any(inner <= zd[1] | inner >= zd[length(zd)])) {
      err <- add(err, "zoneDiv", "Zone boundaries must fall between Start and End.")
    }
  }

  w_msgs <- character(0)
  e_msg <- tryCatch(withCallingHandlers(
    suppressMessages(char_validate_params(char_data, p$pretreatment, p$smoothing,
                                          p$peak_analysis, p$results)),
    warning = function(w) {
      w_msgs <<- c(w_msgs, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
  if (is.character(e_msg)) {
    m   <- sub("^char_validate_params: ", "", e_msg)
    err <- add(err, .cp_key_of(m), m)
  }
  for (m in sub("^char_validate_params: ", "", w_msgs)) {
    warn <- add(warn, .cp_key_of(m), m)
  }

  if (length(zd) >= 2L) {
    end <- zd[length(zd)]
    if (zd[1] < summ$age_top) warn <- add(warn, "zoneDiv", paste0("Start (",
      zd[1], ") is younger than the top of the record (", summ$age_top, ")."))
    if (zd[1] > summ$age_top) warn <- add(warn, "zoneDiv", paste0("Start (",
      zd[1], ") is older than the top of the record (", summ$age_top,
      "); samples above it will be excluded."))
    if (end > summ$age_bot) warn <- add(warn, "zoneDiv", paste0("End (", end,
      ") is older than the bottom of the record (", summ$age_bot,
      "); CharAnalysis will reset it to ", summ$age_bot, "."))
    if (end < summ$age_bot) warn <- add(warn, "zoneDiv", paste0("End (", end,
      ") is younger than the bottom of the record (", summ$age_bot,
      "); samples below it will be excluded."))
  }
  if (v[9] > 0 && v[9] < summ$res_median) {
    warn <- add(warn, "yrInterp", paste0("The interpolation interval (", v[9],
      " yr) is finer than the median sample resolution (",
      round(summ$res_median, 1), " yr), so interpolated values will not be ",
      "independent measurements."))
  }
  list(errors = err, warnings = warn)
}

# -- display helpers -----------------------------------------------------------

# Escape text and subscript variable names: C_peak -> C<sub>peak</sub>.
.cp_sub <- function(x) {
  x <- htmltools::htmlEscape(paste(x, collapse = ""))
  shiny::HTML(gsub("\\bC_([A-Za-z]+)", "C<sub>\\1</sub>", x, perl = TRUE))
}

.cp_box <- function(text, type) {
  shiny::div(class = if (type == "error") "cp-note cp-err" else "cp-note",
             shiny::strong(if (type == "error") "Must fix: " else "Note: "),
             .cp_sub(text))
}

.cp_msgs_ui <- function(chk, key = NULL, prefix = FALSE) {
  e <- chk$errors; w <- chk$warnings
  if (!is.null(key)) { e <- e[e$key == key, ]; w <- w[w$key == key, ] }
  if (nrow(e) + nrow(w) == 0L) return(NULL)
  lab <- function(df) if (prefix) paste0(.cp_key_tab[df$key], ": ", df$text) else df$text
  shiny::tagList(lapply(lab(e), .cp_box, type = "error"),
                 lapply(lab(w), .cp_box, type = "note"))
}

# Info icon with a hover/focus popup. `tip` is a character vector; elements of
# the form "Label: text" are shown with the label in bold.
.cp_tip <- function(tip) {
  if (is.null(tip)) return(NULL)
  lines <- lapply(tip, function(t) {
    if (grepl("^[^:]{1,40}: ", t)) {
      shiny::div(shiny::strong(sub(":.*$", ":", t)), .cp_sub(sub("^[^:]*: ", "", t)))
    } else shiny::div(.cp_sub(t))
  })
  shiny::span(class = "cp-tip", tabindex = "0", `aria-label` = "More information",
              shiny::icon("circle-info"), shiny::span(class = "cp-bub", lines))
}

# One parameter row: name + definition | input | description + messages.
.cp_row <- function(name, def, input, desc, key, tip = NULL) {
  shiny::tags$tr(
    shiny::tags$td(shiny::tags$code(.cp_sub(name)), .cp_tip(tip),
                   shiny::span(class = "cp-def", .cp_sub(def))),
    shiny::tags$td(input),
    shiny::tags$td(class = "cp-desc", .cp_sub(desc),
                   shiny::uiOutput(paste0("msg_", key))))
}

.cp_table <- function(...) {
  shiny::tags$table(class = "cp-tbl",
    shiny::tags$colgroup(shiny::tags$col(style = "width:27%"),
                         shiny::tags$col(style = "width:31%"),
                         shiny::tags$col(style = "width:42%")),
    shiny::tags$tr(shiny::tags$th("Parameter"), shiny::tags$th("Setting"),
                   shiny::tags$th("Description")),
    ...)
}

.cp_html_table <- function(df, header = TRUE) {
  shiny::tags$table(class = "table table-condensed",
    if (header) shiny::tags$thead(shiny::tags$tr(lapply(names(df), shiny::tags$th))),
    shiny::tags$tbody(lapply(seq_len(nrow(df)), function(i)
      shiny::tags$tr(lapply(df[i, ], function(x) shiny::tags$td(.cp_sub(as.character(x))))))))
}

.cp_css <- "
.cp-title{padding:8px 14px;font-weight:600;border-bottom:1px solid #ddd;background:#f7f7f7}
.cp-tbl{width:100%;table-layout:fixed;border-collapse:collapse;margin-top:6px}
.cp-tbl th{text-align:left;font-weight:600;font-size:12px;color:#555;padding:6px;border-bottom:1px solid #bbb}
.cp-tbl td{vertical-align:top;padding:8px 6px;border-bottom:1px solid #e5e5e5;line-height:1.45}
.cp-tbl td .form-group{margin-bottom:0}
.cp-tbl td .form-group > label.control-label{margin:0;padding:0;height:0;display:block}
.cp-stack{margin-bottom:10px}
.cp-stack:last-child{margin-bottom:0}
.cp-tbl code{color:#222;background:transparent;padding:0;font-size:12.5px;font-weight:600}
.cp-def{display:block;color:#444;margin-top:2px}
.cp-desc{color:#555;font-size:12.5px}
.cp-lab{font-size:12px;color:#555;margin:0 0 2px 0;line-height:1.3}
.cp-note{font-size:12px;padding:5px 8px;margin-top:6px;border-left:4px solid #b26b00;background:#fff4e0;color:#5c3b00}
.cp-err{border-left-color:#8b0000;background:#fdecea;color:#8b0000}
.cp-tip{position:relative;display:inline-block;margin-left:5px;color:#555;cursor:help;outline:none}
.cp-tip .cp-bub{visibility:hidden;position:absolute;left:0;top:20px;z-index:50;width:300px;padding:8px 10px;background:#fff;border:1px solid #bbb;border-radius:4px;color:#222;font-size:12px;font-weight:normal;line-height:1.5;box-shadow:0 2px 6px rgba(0,0,0,.15)}
.cp-tip:hover .cp-bub,.cp-tip:focus .cp-bub{visibility:visible}
.cp-tip .cp-bub div{margin-bottom:4px}
.nav-tabs{margin-bottom:4px}
"

.cp_run_gadget <- function(mode, site, data_file, params_file, char_data, summ,
                           init, old, n_logs, button) {
  app <- .cp_gadget_app(mode, site, data_file, params_file, char_data, summ,
                        init, old, n_logs, button)
  shiny::runGadget(app$ui, app$server, stopOnCancel = FALSE,
                   viewer = shiny::dialogViewer(app$title, width = 1000,
                                                height = 850))
}

# Builds the gadget's UI and server (separate from running it, for testing).
.cp_gadget_app <- function(mode, site, data_file, params_file, char_data, summ,
                           init, old, n_logs, button) {
  sr  <- function(x) shiny::span(class = "sr-only", x)   # label for screen readers
  num <- function(id, label, value, ...) {
    shiny::numericInput(id, sr(label), value, width = "100%", ...)
  }
  sel <- function(id, label, choices, value) {
    shiny::selectInput(id, sr(label), choices = choices, width = "100%",
                       selected = as.character(value))
  }
  stacked <- function(lab, input) {
    shiny::div(class = "cp-stack", shiny::p(class = "cp-lab", lab), input)
  }
  yn  <- c("No" = "0", "Yes" = "1")

  zinit   <- init[1:8]; zinit <- zinit[!is.na(zinit) & zinit != -9999]
  z_start <- if (length(zinit) >= 1L) zinit[1] else NA
  z_end   <- if (length(zinit) >= 2L) zinit[length(zinit)] else NA
  z_int   <- if (length(zinit) > 2L) zinit[2:(length(zinit) - 1L)] else numeric(0)
  all_figs <- if (is.na(init[25])) 1 else init[25]   # hidden; kept as is

  summary_table <- data.frame(
    Measure = c("Data file", "Samples", "Depth range (cm)",
                "Age range (cal. yr BP)", "Record length (yr)",
                "Sample resolution (yr per sample)",
                "Samples with missing values", "Breaks in sample depths"),
    Value = c(basename(data_file), summ$n,
              paste(summ$cm_top, "to", summ$cm_bot),
              paste(summ$age_top, "to", summ$age_bot),
              summ$age_bot - summ$age_top,
              sprintf("median %.1f (range %.1f to %.1f)", summ$res_median,
                      summ$res_min, summ$res_max),
              summ$n_missing, summ$n_gaps))

  title <- paste0("CharAnalysis parameters: ", site,
                  if (mode == "new") " (new file)" else " (revising)")

  tab_data <- shiny::tabPanel("Data",
    .cp_html_table(summary_table, header = FALSE),
    if (mode == "revise") shiny::helpText(paste0("Revising ", basename(params_file),
      ". The file holds ", n_logs, " of a maximum of 10 logged parameter sets.")),
    shiny::helpText("Parameters will be saved to ", basename(params_file),
                    ", in the same folder as the data file."),
    shiny::uiOutput("msg_data"))

  tab_pre <- shiny::tabPanel("Pretreatment", .cp_table(
    .cp_row("zoneDiv", "Analysis period and zone boundaries (cal. yr BP).",
      shiny::tagList(
        stacked("Start", num("z_start", "Start of record", z_start)),
        stacked("Zone boundaries (optional; separate multiple values with commas)",
                shiny::textInput("z_bounds", sr("Zone boundaries"),
                                 value = paste(z_int, collapse = ", "),
                                 placeholder = "e.g., 1000, 3500, 5500",
                                 width = "100%")),
        stacked("End", num("z_end", "End of record", z_end))),
      paste(if (mode == "new") "Start and end are pre-filled from the data (top age of the youngest sample, bottom age of the oldest)."
            else "Values from the existing parameter file.",
            "Enter one or more zone boundaries, youngest first, separated by commas (e.g., 1000, 3500, 5500). Leave blank for a single zone."),
      "zoneDiv",
      tip = c("Zones divide the record for plotting and for comparing fire return intervals between periods.",
              "Up to six boundaries can be entered; they must fall between the start and end.")),
    .cp_row("yrInterpolate", "Time step the record is resampled to (yr).",
      num("yrInterp", "yrInterpolate", init[9], min = 0),
      sprintf(paste("Counts, volumes, and depths are interpolated before CHAR is",
                    "calculated. Median sample resolution is %.1f yr (suggested: %d).",
                    "Enter 0 to use the median automatically."),
              summ$res_median, as.integer(summ$res_suggest)),
      "yrInterp",
      tip = c("Resampling to equal time steps puts every sample on the same time scale, which the smoothing and threshold steps require.",
              "An interval finer than the sample resolution creates values that are not independent measurements.")),
    .cp_row("transform", "Whether CHAR is transformed before analysis.",
      sel("transform", "transform", c("0: No transformation" = "0",
          "1: Base-10 log" = "1", "2: Natural log" = "2"), init[10]),
      "For log transformation, 1 is added to all CHAR values first.",
      "transform",
      tip = "Log transformation compresses large CHAR values relative to small ones.")))

  tab_smo <- shiny::tabPanel("Smoothing", .cp_table(
    .cp_row("method", "How low-frequency CHAR (C_background) is estimated.",
      sel("method", "Smoothing method", c("1: Lowess" = "1",
          "2: Lowess, robust to outliers" = "2", "3: Moving average" = "3",
          "4: Moving median" = "4", "5: Moving mode" = "5"), init[11]),
      paste("The smoother used to estimate the slowly varying background",
            "component of the record. Hover over the icon for what each option does."),
      "method",
      tip = c("Lowess: locally weighted regression within the window.",
              "Lowess, robust to outliers: lowess that iteratively down-weights outliers, reducing the influence of peaks on C_background.",
              "Moving average: mean of the values in the window.",
              "Moving median: running median of the window, then smoothed with lowess.",
              "Moving mode: running mode of the window (from a 100-bin histogram), then smoothed with lowess.")),
    .cp_row("yr", "Window over which C_background is estimated (yr).",
      num("smoYr", "Smoothing window", init[12], min = 1),
      paste("Years to smooth the record over. C_background sensitivity (Results",
            "tab) shows how results change with this window."),
      "smoYr",
      tip = "For local thresholds, this window also sets the length of the moving window in which thresholds are defined.")))

  tab_peak <- shiny::tabPanel("Peaks", .cp_table(
    .cp_row("cPeak", "How high-frequency CHAR (C_peak) is calculated.",
      sel("cPeak", "cPeak", c("1: Residuals" = "1", "2: Ratios" = "2"), init[13]),
      paste("Residuals: C_peak = C_interpolated - C_background.",
            "Ratios: C_peak = C_interpolated / C_background."),
      "cPeak",
      tip = c("Residuals: C_peak is measured in CHAR units above C_background.",
              "Ratios: C_peak is measured relative to C_background.")),
    .cp_row("threshType", "Whether one threshold or locally varying thresholds separate peaks from noise.",
      sel("threshType", "threshType", c("1: Global" = "1", "2: Local" = "2"), init[14]),
      paste("Global: one threshold for the whole record. Local: thresholds are",
            "defined within a moving window (equal to the length selected in",
            "Smoothing) and thus vary through time."),
      "threshType",
      tip = c("Global: one noise distribution and threshold for the whole record.",
              "Local: a noise distribution and threshold are estimated within each window, so thresholds vary through time.")),
    .cp_row("threshMethod", "How threshold values are determined.",
      sel("threshMethod", "threshMethod", c("1: User-defined" = "1",
          "2: Gaussian noise" = "2", "3: Gaussian mixture" = "3"), init[15]),
      paste("For methods 2 and 3, threshold values are percentiles of the noise",
            "distribution (e.g., 0.95). For method 1, values are in C_peak units."),
      "threshMethod",
      tip = c("User-defined: you supply thresholds in C_peak units (global thresholds only).",
              "Gaussian noise: noise is modeled as a Gaussian distribution with mean 0 (residuals) or 1 (ratios).",
              "Gaussian mixture: C_peak values are modeled as a mixture of two Gaussian populations, noise and signal; thresholds are percentiles of the noise population.")),
    .cp_row("threshValues", "Thresholds evaluated; Peaks Final is used for analysis.",
      shiny::tagList(
        stacked("Peaks 1", num("t1", "Peaks 1 threshold", init[16])),
        stacked("Peaks 2", num("t2", "Peaks 2 threshold", init[17])),
        stacked("Peaks 3", num("t3", "Peaks 3 threshold", init[18])),
        stacked("Peaks Final", num("t4", "Peaks Final threshold", init[19]))),
      paste("Peaks 1-3 are for comparison. Peaks Final is used for peak",
            "identification and all calculations in the Results tab."),
      "threshValues",
      tip = "Names match the peaks 1, peaks 2, peaks 3, and peaks Final columns in the _charResults.csv output."),
    .cp_row("minCountP", "Cut-off for screening peaks that could arise from counting noise.",
      num("minCountP", "minCountP", init[20], min = 0, max = 1),
      paste("A peak is flagged and excluded if the minimum count in the 75 yr",
            "before it has more than this probability of coming from the same",
            "Poisson distribution as the peak's maximum count. Set to 0.99 to turn off."),
      "minCountP",
      tip = "Screens out peaks that could result from sampling variability when charcoal counts are low. Flagged peaks are still plotted.")))

  tab_res <- shiny::tabPanel("Results", .cp_table(
    .cp_row("peakFrequ", "Window for smoothing fire frequency and fire return intervals (yr).",
      num("peakFrequ", "peakFrequ", init[21], min = 1),
      "Years over which fire frequency and fire return intervals, from Peaks Final, are smoothed.",
      "peakFrequ"),
    .cp_row("C_background sensitivity", "Whether the analysis is repeated across smoothing windows.",
      sel("bkgSens", "C_background sensitivity", yn, init[22]),
      "Shows how sensitive results are to the C_background window. Recommended for the first analysis of a record.",
      "bkgSens",
      tip = "Repeats the analysis across a range of smoothing windows and plots how the results change."),
    .cp_row("saveData", "Save the results table.",
      sel("saveData", "saveData", yn, init[24]),
      paste0("Writes ", site, "_charResults.csv to the folder holding the parameter ",
             "file, replacing any earlier results file. Results can always be ",
             "recreated from the saved parameters."),
      "saveData"),
    .cp_row("saveFigures", "Save the figures.",
      sel("saveFigures", "saveFigures", yn, init[23]),
      "Writes all figures as PDF files to the same folder, replacing earlier copies.",
      "saveFigures")))

  tab_rev <- shiny::tabPanel("Review", shiny::uiOutput("review"))

  ui <- miniUI::miniPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML(.cp_css))),
    shiny::div(class = "cp-title", title),
    miniUI::miniContentPanel(
      shiny::tabsetPanel(id = "tabs", type = "tabs", selected = "Pretreatment",
                         tab_data, tab_pre, tab_smo, tab_peak, tab_res, tab_rev)),
    miniUI::miniButtonBlock(
      shiny::actionButton("cancel", "Cancel"),
      shiny::actionButton("done", button, class = "btn-primary"),
      border = "top"))

  server <- function(input, output, session) {
    nz <- function(x) if (is.null(x) || length(x) == 0L) NA_real_ else suppressWarnings(as.numeric(x))

    bounds <- shiny::reactive(.cp_parse_bounds(input$z_bounds))
    vals <- shiny::reactive({
      zd <- c(nz(input$z_start), bounds()$values, nz(input$z_end))
      zd <- c(zd, rep(-9999, max(0L, 8L - length(zd))))[1:8]
      unname(c(zd,
        vapply(c("yrInterp", "transform", "method", "smoYr", "cPeak",
                 "threshType", "threshMethod", "t1", "t2", "t3", "t4",
                 "minCountP", "peakFrequ", "bkgSens", "saveFigures", "saveData"),
               function(id) nz(input[[id]]), numeric(1)),
        all_figs))
    })
    chk <- shiny::reactive(.cp_check(vals(), summ, char_data, bounds()$error))
    changes <- shiny::reactive({
      if (is.null(old)) return(NULL)
      v <- vals()
      d <- which(is.na(old) != is.na(v) | (!is.na(old) & !is.na(v) & old != v))
      if (length(d) == 0L) return(data.frame())
      data.frame(Parameter = .cp_param_labels[d],
                 Previous  = .cp_fmt_num(old[d]),
                 New       = ifelse(is.na(v[d]), "(empty)", .cp_fmt_num(v[d])))
    })

    for (k in unique(c("data", .cp_pos_key))) local({
      key <- k
      output[[paste0("msg_", key)]] <- shiny::renderUI(.cp_msgs_ui(chk(), key))
    })

    output$review <- shiny::renderUI({
      c0 <- chk(); ch <- changes()
      shiny::tagList(
        shiny::p(if (nrow(c0$errors) > 0L)
          shiny::strong(nrow(c0$errors), " problem(s) must be fixed before saving.")
          else "No problems prevent saving."),
        .cp_msgs_ui(c0, prefix = TRUE),
        if (!is.null(ch)) shiny::tagList(
          shiny::h4("Changes from the existing file"),
          if (nrow(ch) == 0L) shiny::p("No values changed; the file will not be rewritten.")
          else .cp_html_table(ch)))
    })

    do_save <- function() {
      v   <- vals()
      res <- tryCatch({
        if (mode == "new") char_params_new(v, params_file)
        else char_params_revise(params_file, v)
        TRUE
      }, error = function(e) conditionMessage(e))
      if (isTRUE(res)) {
        shiny::stopApp(normalizePath(params_file, winslash = "/"))
      } else {
        shiny::showModal(shiny::modalDialog(title = "Could not save", res,
                                            easyClose = TRUE))
      }
    }

    shiny::observeEvent(input$done, {
      c0 <- chk()
      if (nrow(c0$errors) > 0L) {
        shiny::showModal(shiny::modalDialog(
          title = "Fix these problems before saving",
          shiny::tags$ul(lapply(paste0(.cp_key_tab[c0$errors$key], ": ", c0$errors$text),
                                function(t) shiny::tags$li(.cp_sub(t)))),
          easyClose = TRUE))
        return()
      }
      ch <- changes()
      if (mode == "revise" && nrow(ch) == 0L) {
        shiny::stopApp(normalizePath(params_file, winslash = "/"))  # unchanged
      } else if (mode == "revise") {
        shiny::showModal(shiny::modalDialog(
          title = "Confirm changes",
          shiny::p("The previous values will be kept in a new, time-stamped",
                   " Parameter Log column."),
          .cp_html_table(ch),
          footer = shiny::tagList(shiny::modalButton("Go back"),
                                  shiny::actionButton("confirm", button,
                                                      class = "btn-primary"))))
      } else {
        do_save()
      }
    })
    shiny::observeEvent(input$confirm, { shiny::removeModal(); do_save() })
    shiny::observeEvent(input$cancel, shiny::stopApp(NULL))
  }

  list(ui = ui, server = server, title = title)
}
