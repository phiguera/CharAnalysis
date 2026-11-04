# *CharAnalysis* Parameter Interface: Development Plan
*Last updated: November 2026*

---

## Purpose

Running *CharAnalysis* requires a parameter file (`<site>_charParams.csv`) with
fixed row positions, index-coded options, and conventions (e.g., `-9999` padding
in `zoneDiv`) that are easy to get wrong, particularly for users new to R. This
plan describes a graphical interface, run from within R or RStudio, that guides
users through every parameter and writes a valid parameter file. The analysis
still runs in R; the interface replaces hand-editing the parameter file, not the
R package. The goals are to lower the barrier to use (install the package and
type `CharAnalysis()`), keep the analysis transparent (summarize the data and
explain each choice), and serve as an on-ramp to R (produce a standard parameter
file and print the command that reproduces the analysis).

---

## User workflow

1. The user types `CharAnalysis()` with no argument.
2. A console prompt offers three options: run an existing parameter file (the
   current file picker), revise an existing parameter file, or create a new one.
3. To create a new file, the user selects a data file, which must be named
   `<site>_charData.csv`; the site name is taken from the file name. If the file
   does not follow this convention, the interface stops and explains the naming
   requirement. If `<site>_charParams.csv` already exists, the interface offers
   to revise that file instead.
4. The interface opens (RStudio Viewer or default browser), showing a summary of
   the data record and a form containing every parameter, grouped by stage.
   When revising, the form is pre-filled with the file's current values.
5. On saving, the interface writes `<site>_charParams.csv` to the folder holding
   the data file, prints the equivalent command
   (`CharAnalysis("<path>/<site>_charParams.csv")`), and runs the analysis.

**Revising a file.** If any value has changed, the interface lists the changes
and asks for confirmation. The previous values are then moved into a new
Parameter Log column, headed with a time stamp and the package version (e.g.,
`Replaced 2026-11-03 14:32, CharAnalysis 2.1.0`), and the new values are written
to the Parameters column. The newest log column sits next to the current values;
older log columns, including any entered by hand, are kept unchanged to the
right. If nothing has changed, the file is not rewritten. A file can hold up to
10 logged parameter sets: the MATLAB reader parses a fixed 16 columns per row,
and 10 log columns is the most that keeps files readable by both versions. When
a file reaches this limit, the interface stops and asks the user to save a copy
and remove older log columns before revising.

The interface will also be available directly as an exported function,
`char_params_gui()`. Behaviour in non-interactive sessions is unchanged.

---

## Design principles

- **The parameter file remains the single record of an analysis.** The interface
  writes the same structure as `template_charParams.csv` (all rows, labels,
  descriptions, and `-9999` padding), so files remain readable by the R and
  MATLAB versions, both of which read parameters by row position, and editable
  by hand. Nothing is stored elsewhere. In new
  files, the Parameter Log columns are left empty.
- **Previous parameters are never lost.** An existing file is changed only
  through the revision path above, after confirmation, with prior values moved
  to the Parameter Log. Files are written to a temporary copy and then swapped in,
  so an interrupted save cannot leave a partial file.
- **One set of validation rules.** The interface uses `char_validate_params()`,
  so interface-built and hand-built files are checked identically.
- **No new required dependencies.** `shiny` and `miniUI` are suggested packages,
  installed on first use if missing. The interface runs on Windows and macOS,
  from the R console or RStudio.

---

## Interface content

| Section | Parameters | Controls and guidance |
|---|---|---|
| Data summary | (read-only) | Number of samples; depth and age range; sample resolution (median and range); gaps and missing values. |
| Pretreatment | `zoneDiv`, `yrInterpolate`, `transform` | Only the start and end of the record are pre-filled (top age of the youngest sample and bottom age of the oldest sample); both are editable. Internal zone boundaries (up to six) are entered by the user and are never suggested. `-9999` padding handled automatically. Interpolation interval suggested from the median sample resolution. Warnings if the interval is finer than the median resolution or a boundary falls outside the record. |
| Smoothing | `method`, `yr` | Smoothing options listed by name rather than index. Warning for implausible windows (e.g., longer than the record or a zone). |
| Peak analysis | `cPeak`, `threshType`, `threshMethod`, `threshValues`, `minCountP` | Threshold labels and units change with the threshold method (percentiles vs. C_peak units); the form states that the last value is used for peak analysis. |
| Results | `peakFrequ`, background sensitivity, `saveFigures`, `saveData`, `allFigures` | Background sensitivity recommended for the first analysis of a record; note that `allFigures` has no effect in R. |

Each field carries inline help adapted from the template descriptions, with
references where relevant (Higuera et al. 2009, 2010; Kelly et al. 2011).
Warnings advise but do not block; only values that `char_validate_params()`
rejects prevent saving. The existing safeguards after an analysis remain
unchanged: the signal-to-noise index (SNI) advisory and diagnostic figures 1 to 5
before the analytical figures.

---

## Development phases

**Status (November 2026):** Phases 1-3 are implemented in R and MATLAB on the
`dev` branch and were tested interactively by the author (RStudio; MATLAB
R2022a). Documentation in the vignette and User's Guide (Section 2.3) is
written. Phase 4 testing with new users and the release remain.

**Phase 1: Parameter writer.** An internal function that writes a
template-structured parameter file from a list of values. Round-trip tests for
the CO, CH10, SI17, and RA07 validation files confirm that read, write, and
re-read produce identical parameters and identical analysis results. Further
tests confirm that revisions add log columns correctly, preserve hand-entered
log content, and leave files readable by both the R and MATLAB readers.

**Phase 2: Interface.** The parameter form, the data summary, the
existing-or-new prompt in `CharAnalysis()`, and the exported `char_params_gui()`.

**Phase 3: Guidance.** Data-informed suggestions, warnings, and inline help.

**Phase 4: Testing and release.** Testing with graduate students and researchers
new to R; revision; documentation in the vignette and User's Guide; release in a
minor version.

---

## MATLAB version

The MATLAB version (`CharAnalysis_2_0_MATLAB`) follows the same workflow:
`CharAnalysis` with no input offers the same three options, and
`CharParamsGUI` (a `uifigure` window; MATLAB R2018b or later) creates or
revises the parameter file. The writer functions (`CharParamsNew`,
`CharParamsRevise`) produce the same file layout and Parameter Log as the R
writer, so files are interchangeable. Differences from R: `allFigures` is
shown, because in MATLAB it controls Figures 1, 2, and 9; per-parameter notes
are listed at the bottom of each tab, with the affected inputs highlighted,
rather than inside each row; and saving results and figures follows the
existing MATLAB behaviour (PDF and TIFF figures), now written to the folder
holding the parameter file.

---

## Open questions

- **Excel input.** The reader accepts `.xlsx` files; the interface is planned to
  write CSV only.
---

*This plan reflects current intentions and is subject to change. Feedback and
collaboration are welcome; please use the Issues tab at
https://github.com/phiguera/CharAnalysis/issues.*

*This plan was developed with the assistance of Claude, an AI assistant by
Anthropic, and has been reviewed by the author.*
