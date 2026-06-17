# *CharAnalysis* 2.0.4 (development version)

## Development scripts (ROADMAP Section 3 — Chronological Uncertainty)

Three exploratory scripts added to `tests/` (not yet exported functions):

- **`char_run_ensemble.R`**: runs the CharAnalysis pipeline across an N-member
  age-depth ensemble. Saves output as `CH10_ensemble_results.rds`.
- **`run_ensemble_analysis.R`**: nearest-neighbor peak matching that computes
  P(peak) and timing confidence intervals for each reference peak.
  Matching half-window = max(mean Weibull mFRI / 2, `yrInterp`); a
  reference-peak-first loop guarantees P(peak) ≤ 1.0.
- **`plot_ensemble_figure.R`**: four-panel ensemble figure. Panel (c) shows
  a beeswarm of circles sized by P(peak) on a universal fixed 5-bin scale
  (0–20, 21–40, 41–60, 61–80, 81–100%). Panel (d) shows jittered 95% CI bars on peak timing
  with a dot at the reference peak age. See script header for full methods.

Applied to CH10 (Chickaree Lake): 59 reference peaks, P(peak) 0.58–1.00,
typical 95% CI widths 40–80 yr.

## Breaking changes

Output figures are reorganized into **diagnostic** (Figs 1-5) and
**analytical** (Figs 6-9) categories. Figure numbers and PDF filenames have
changed from the previous release. MATLAB figure numbering is unaffected; see
the User's Guide (Section 7.2) for a note on the R/MATLAB divergence.

| Figure | Function | Previous number |
|--------|----------|-----------------|
| 1 | `char_plot_raw()` | 1 (unchanged) |
| 2 | `char_plot_thresh_diag()` | 2 (unchanged) |
| 3 | `char_plot_peaks()` | 3 (unchanged) |
| 4 | `char_plot_sni()` | 4 (unchanged) |
| **5** | `char_bkg_sensitivity()` | **10** |
| **6** | `char_plot_zones()` | **8** |
| **7** | `char_plot_cumulative()` | **5** |
| **8** | `char_plot_fri()` | **6** |
| **9** | `char_plot_fire_history()` | **7** |

PDF filenames saved via `char_plot_all()` or the new wrappers reflect the new
numbers (e.g. `site_05_sensitivity_to_C_background.pdf`,
`site_06_zone_comparisons.pdf`).

The `allFigures` parameter (parameter-file row 25) is still read without error
but has no effect on figure output in the R package. All diagnostic figures
(Figs 1-5) are now produced whenever `char_plot_diagnostic()` is called.
Use `plots = FALSE` in `CharAnalysis()` to suppress automatic figure output.

## New features

- Two new wrapper functions expose the diagnostic/analytical split directly:
  `char_plot_diagnostic(out)` produces Figs 1-5 (parameter evaluation);
  `char_plot_analysis(out)` produces Figs 6-9 (fire-history interpretation).
  `char_plot_all()` is retained as a convenience wrapper that calls both.

- `CharAnalysis()` gains a `plots` argument (default `TRUE`). When `TRUE`,
  diagnostic figures (Figs 1-5) are produced automatically at the end of the
  pipeline via `char_plot_diagnostic()`. Set `plots = FALSE` to suppress all
  automatic figure output in scripts or batch jobs. Analytical figures are
  always produced by an explicit call to `char_plot_analysis(out)`.

- `CharAnalysis()` now emits a styled advisory (via the `cli` package) at the
  end of all console output when the signal-to-noise index (SNI) falls below
  3.0, the minimum value recommended by Kelly et al. (2011) to identify records
  suitable for peak detection. For the local-threshold path (`threshType != 1`)
  the message reports the percentage of samples below the threshold; for the
  global-threshold path (`threshType = 1`) it reports the single record-wide
  SNI value. In both cases the user is directed to `char_plot_sni(out)` for a
  visual assessment. The advisory is rendered with a red bordered block and
  colour-coded bullets (`cli_alert_danger` / `cli_alert_info`), with graceful
  plain-text fallback in environments that do not support ANSI colour.

- New function `char_bkg_sensitivity()` ports the background-window
  sensitivity analysis from the MATLAB v2.0 `bkgCharSensitivity.m`
  (MATLAB Figure 10; R Figure 5). It re-runs the smoothing -> C_peak -> threshold ->
  peak-identification pipeline across a range of background smoothing
  windows and summarises the result. On the global-threshold path it
  produces a filled-contour surface of peaks identified versus threshold
  and window width; on the local-threshold path it produces a four-panel
  diagnostic (noise goodness-of-fit, signal-to-noise index, summed median
  SNI + GOF, and mean fire-return interval with peak count). Figures use
  `ggplot2`/`patchwork`, consistent with the other `char_plot_*` figures,
  and save as a single PDF.
- `CharAnalysis()` now honours the `bkgSens` parameter (parameter-file
  row 23). When `bkgSens == 1` the sensitivity analysis runs automatically
  and its results are attached to the returned object as
  `out$bkg_sensitivity`; the figure (Fig 5) is shown automatically when
  `plots = TRUE` (default) in `CharAnalysis()`, and is never auto-saved.
- `char_thresh_global()` gains an optional `thresh_bins` argument
  (default `NULL`, preserving previous behaviour). It lets the threshold
  grid be held fixed across smoothing windows, the explicit R realisation
  of MATLAB's `bkgSensIn` mechanism. No change to existing results.

## Chronological uncertainty — exploratory scripts (ROADMAP section 2)

Two scripts in `tests/` begin the implementation of ROADMAP section 2
(propagating chronological uncertainty into fire-event characterization).
These are not yet exported package functions; they demonstrate the design
and serve as the basis for future development.

- `tests/char_run_ensemble.R`: runs the full CharAnalysis pipeline across a
  user-supplied ensemble of age-depth chronologies. For each iteration, one
  chronology is substituted into the input data, depths are matched to the
  chronology matrix by the `Sample_cm` column (with linear interpolation at
  the 16 charData depths lacking an exact match), and peak detection is
  re-run. The binary `peaksFinal` result is snapped onto a common time grid
  (defined from the median start and end ages across the ensemble, at the
  fixed `yrInterp` from a reference run). Output is a `CharEnsemble` list
  containing `peaks_matrix` (time steps × iterations), per-iteration
  `charPeak_matrix` and `charBkg_matrix`, and `prob_peak` (fraction of
  iterations identifying a peak at each time step).

- `tests/plot_prob_peak.R`: four-panel `ggplot2` + `patchwork` figure
  summarising the ensemble. Panel (a) replicates the
  C_interpolated / C_background / peak-ID panel from `char_plot_sni()`;
  panel (b) is the signal-to-noise index (SNI); panel (c) shows reference
  peaks as proportional circles scaled by P(peak); panel (d) shows P(peak)
  as grey bars with the same peak markers used in panel (a).

An example 1000-member chronology ensemble for Chickaree Lake (CH10,
Dunnette et al. 2014) is included as
`tests/CH10_MCAgeDepth_1000_chronologies_2026_06_12.csv`. The ensemble was
generated by `MC_AgeDepth.m`, a MATLAB program written by P. Higuera for
Monte Carlo sampling of age-depth models. Future development will add
equivalent input from `rbacon` and other R-based Bayesian age-depth programs.

## Bug fixes

- `char_post_process()` now handles the global-threshold path
  (`threshType = 1`) correctly. Previously a global run errored in
  post-processing with `subscript out of bounds`, because the number of
  reported threshold columns was taken from `ncol(charPeaks)` (the full
  positive candidate-threshold grid, ~hundreds of columns) instead of the
  number of `threshValues`, and the single global negative-threshold
  column was indexed as though it had one column per `threshValue`. Both
  the column count and the negative-threshold index are now derived from
  the number of `threshValues`. Local-threshold results are unchanged.
  This path was not previously exercised end to end.

# *CharAnalysis* 2.0.3

Patch release with substantive vignette improvements and several
factual corrections to the worked example. No changes to analytical
behaviour or to any function signature.

- The vignette now live-renders all worked-example figures from the
  bundled Code Lake dataset; the previous build showed only static
  code blocks. Figure chunks are gated on the availability of
  `ggplot2`, `patchwork`, and `ggtext` (all in `Suggests`) and skip
  cleanly if any is missing.
- Two diagnostic figures, `char_plot_raw()` (Fig. 1) and
  `char_plot_thresh_diag()` (Fig. 2), are now included in the
  vignette alongside the other figures.
- A new bundled parameter file `CO_compensated_charParams.csv`
  (with companion `CO_compensated_charData.csv`) is shipped in
  `inst/validation/`. It is identical to the standard `CO_charParams.csv`
  except that the working threshold percentile (`threshValues[4]`) is
  lowered from 0.99 to 0.95 to compensate for language-induced drift in
  the Gaussian mixture model (GMM) noise estimation. The vignette uses
  this file for the worked example. With this compensation, R identifies
  50 peaks for Code Lake (close to the published MATLAB v2.0 result of
  48 peaks) and reproduces the significant decrease in fire-return
  intervals from Zone 1 to Zone 2 reported in Higuera et al. (2009).
  The strict 1-to-1 reference configuration remains available as
  `CO_charParams.csv`.
- The vignette's *Comparison with MATLAB v2.0* section gained a
  `Threshold` column in the validation table and a new "CO (compensated)"
  row; the narrative below the table was rewritten in terms of
  compensation for GMM drift rather than ad-hoc tuning.
- Vignette factual corrections: Code Lake is in Alaska (not Colorado);
  the smoothing-method label for Code Lake in the validation table is
  Method 4 (moving median), not Method 1 (lowess); the discussion of
  smoothing-related differences was adjusted to identify Method 2
  (robust lowess) as the only method that diverges between R and MATLAB.
- Fixed a stale `system.file()` path in the vignette
  (`extdata` → `validation`) that would have returned an empty string
  if the chunk had been evaluated.

# *CharAnalysis* 2.0.2

Patch release addressing the third round of CRAN reviewer feedback,
plus a related rendering fix discovered during local testing. No
changes to analytical behaviour.

- `Title` field shortened to "Peak Detection and Fire History from
  Sediment-Charcoal Records" (62 characters) to satisfy the CRAN
  convention of titles under 65 characters.
- `char_write_results()` and `char_plot_all()` no longer default
  `out_dir` to the working directory. `out_dir` is now required for
  `char_write_results()` and required when `save = TRUE` for
  `char_plot_all()`. This brings the package into compliance with the
  CRAN policy against writing to the user's home filespace by default.
  Users should pass an explicit path; `tempdir()` is acceptable for a
  transient location.
- Vignette updated to write its example output to `tempdir()` instead of
  a `Results/` directory in the user's working directory, with a note
  that real users would substitute their own path.
- Fixed axis-label rendering in `char_plot_peaks()`,
  `char_plot_fire_history()` (peak-magnitude panel), and
  `char_plot_zones()` (CDF and box-plot panels). These labels were
  still using `expression(paste(...))` (plotmath) syntax and rendered
  as raw text under the `ggtext::element_markdown()` axis-title theme
  introduced in 2.0.1. They now use the same conditional HTML-tag /
  plain-text pattern as the FRI and fire-frequency labels, so super-
  and subscripts render correctly when `ggtext` is available.

# *CharAnalysis* 2.0.1

Patch release addressing CRAN reviewer feedback on the initial submission.
No changes to analytical behaviour.

- DESCRIPTION: condensed to a single paragraph and removed the paragraph
  separators that were rendering as double periods in CRAN metadata.
- `char_parameters()` is now exported. Its help page previously contained
  an example for an unexported function.
- Replaced `\dontrun{}` wrappers in all examples with `\donttest{}` (or
  unwrapped entirely, where the example runs in under 5 seconds). All
  examples now use the bundled validation dataset via `system.file()` and
  write any output to `tempdir()`.

# *CharAnalysis* 2.0.0

First R implementation of *CharAnalysis*, a direct translation of
*CharAnalysis* v2.0 (MATLAB). Analytical outputs are validated against
four benchmark datasets; user testing is ongoing. The package is in the
[experimental](https://lifecycle.r-lib.org/articles/stages.html#experimental)
lifecycle stage: the API may change as user feedback is incorporated.

Please report issues at <https://github.com/phiguera/CharAnalysis/issues>.

## New features relative to MATLAB v2.0

- Full R package with `devtools::install_github()` installation.
- `CharAnalysis()` returns a named list with S3 class `"CharAnalysis"`;
  `print()`, `summary()`, and `plot()` methods are provided.
- Eight publication-quality ggplot2 figures (`char_plot_peaks()`,
  `char_plot_fire_history()`, `char_plot_cumulative()`, `char_plot_fri()`,
  `char_plot_zones()`, `char_plot_raw()`, `char_plot_thresh_diag()`,
  `char_plot_sni()`); save all to PDF with `char_plot_all()`.
- `char_write_results()` writes the 33-column output matrix to CSV with
  column headers and numeric format matching the MATLAB output exactly.

## Known differences from MATLAB v2.0

**Robust lowess background (smoothing method 2).** MATLAB's Curve Fitting
Toolbox `smooth(..., 'rlowess')` and the R `char_lowess()` port handle
NaN gaps differently inside the bisquare robustness iteration, producing
slightly different C_background series for records with missing values.
For gap-free records the difference is negligible (< 0.001); for records
with NaN gaps the maximum absolute difference across validated datasets
is 0.267 (Chickaree Lake). Smoothing method 1 (plain lowess) is
unaffected and agrees to within floating-point noise.

**GMM peak counts.** The Gaussian mixture model EM algorithm accumulates
floating-point arithmetic differently in R and MATLAB, causing slightly
different threshold values in some local windows. Peak counts differ by
10–20% across validated datasets, with the direction varying by dataset.

**MATLAB figures not implemented in R.** The threshold-sensitivity detail
plot (MATLAB Fig. 9) and multi-site comparison plot (MATLAB Fig. 10) are
not implemented in this R package. The background-sensitivity analysis
(MATLAB Fig. 10) is implemented as `char_bkg_sensitivity()` (R Fig. 5)
but produces a different layout than the MATLAB version.

**Smoothed FRI column (col 23).** The R package computes smoothed
fire-return intervals in this column; MATLAB v2.0 stores NaN.

## Validation summary

| Dataset | Site | charBkg max\|diff\| | Peaks R | Peaks MATLAB |
|---------|------|-------------------|---------|-------------|
| CO | Code Lake, AK | ~5 × 10^-6^ | 39 | 48 |
| CH10 | Chickaree Lake, CO | 0.267 | 59 | 50 |
| SI17 | Silver Lake, CO | 0.130 | 25 | 31 |
| RA07 | Raven Lake, AK | < 0.001 | 15 | 17 |

Full validation details are in `inst/z_Validation_report_R_vs_MATLAB_V_2.0.md`.

## Citation

If you use *CharAnalysis* in published research, please cite:

> Higuera, P.E., L.B. Brubaker, P.M. Anderson, F.S. Hu, and T.A. Brown.
> 2009. Vegetation mediated the impacts of postglacial climate change on
> fire regimes in the south-central Brooks Range, Alaska. *Ecological
> Monographs* 79:201–219. <https://doi.org/10.1890/07-2019.1>

If you used Version 2.0 specifically, please also cite the software:

> Higuera, P.E. 2026. *CharAnalysis*: Diagnostic and analytical tools for
> peak analysis in sediment-charcoal records (Version 2.0). Zenodo.
> <https://doi.org/10.5281/zenodo.19304064>
