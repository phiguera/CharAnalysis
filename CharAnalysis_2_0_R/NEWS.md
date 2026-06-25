# *CharAnalysis* 2.0.4 (development version)

## Development scripts (ROADMAP Section 3 — Chronological Uncertainty, continued)

### CO validation site; ensemble bug fixes; age-depth figure (2026-06-24)

**`char_run_ensemble.R`** — three bug fixes:

- **NA bridging** in `run_one_iteration`: after `char_pretreatment`, any boundary
  NAs in `charcoal$accI` are now repaired with `approx(..., rule = 2L)` before
  `char_smooth` is called. Root cause: for some MC chronologies the resampled age
  grid extends slightly beyond the trimmed data, leaving the first or last
  resampled interval with no overlapping raw sample. `char_smooth` bridges these
  internally in `acc_clean`, but `charcoal$peak = accI - accIS` uses the raw
  `accI`, so the NA propagated to `char_post_process` and caused `if (NA)` errors.
- **Stale workspace fix**: `out <- CharAnalysis(params_file, plots = FALSE)` is
  now called explicitly in `char_run_ensemble.R` before sourcing
  `run_ensemble_analysis.R`. Previously, `run_ensemble_analysis.R` guarded its
  `CharAnalysis()` call with `if (!exists("out"))`, which silently reused a stale
  object from a prior site if one remained in the workspace.
- **Input validation**: column class check added after `read.csv(chron_file)`;
  catches the common mistake of selecting a `*_peakAgeUncertainty.csv` output
  file instead of `*_chronologies.csv`, with an informative error message.

**CO (Code Lake, Alaska) — third validation site** added to the chronological
uncertainty ensemble workflow:

- Results: 50 reference peaks (arithmetic mean FRI 148 yr), median detection
  frequency 98.4%, median peak timing SD 21 yr — comparable to CH10 in
  chronological precision. Matching windows narrow (median ±64 yr), dominated
  by the mFRI/2 floor. 12 ensemble-only peaks (6 near-reference, 6 independent).
- Files added: `tests/CO_MCAgeDepth_1000_chronologies.csv`,
  `tests/CO_ensemble_results.rds`, `tests/CO_peakAgeUncertainty.csv`,
  `tests/CO_ensemble_figure.png`, `inst/validation/CO_peakAgeUncertainty.csv`.

**`tests/plot_agedepth_vignette.R`** — updated from two to three sites:

- CO added in green (`#1b7837`); CH10 blue, SI17 black.
- **Axis offset bug fixed**: replaced `limits = c(x_hi, x_lo)` + `expand` in
  `scale_x_reverse()` with `expand = expansion(0)` + `coord_cartesian(xlim = ...)`.
  The `limits`/`expand` combination on reversed ggplot2 scales mis-applies padding
  in data-coordinate space, shifting all lines away from the surface axis edge by
  ~2% of the depth range. `coord_cartesian` sets visual bounds post-reversal and
  is unambiguous. General rule: always use `coord_cartesian()` to set display
  limits on reversed axes, not `limits` + `expand`.

**Vignette** `vignettes/chronological_uncertainty_draft.md`:

- Section 3.3 (CO) added as a third worked example.
- Section 3 framing paragraph and summary table updated to include CO.
- Figure 2 (age-depth comparison) updated: CO ribbon added; caption corrected
  from "CO (red)" to "CO (green)".

## Bug fix (2026-06-24)

- `run_ensemble_analysis.R`: fixed `object 'mFRI_floor' not found` error in the
  Section 7b merging pass; `mFRI_floor` (removed in Round 8 per-zone refactor)
  replaced with `mFRI_floor_for_age(median(grp_ages))`. Ensemble-only peak count
  for CH10 updated from 10 to 9 as a result.

## Development scripts (ROADMAP Section 3 — Chronological Uncertainty, continued)

### Vignette finalization: methods figure, subscripts, title, x-axis crop (2026-06-23, continued)

**Vignette** `vignettes/chronological_uncertainty_draft.md` — final round of edits:

- Title changed to "Integrating Chronological Uncertainty into *CharAnalysis*".
- *CharAnalysis* italicized in all prose instances (skipping code blocks and paths).
- Step 3 gained a third bullet for ensemble-only peaks, including the 10%
  detection frequency threshold and note that it is user-adjustable.
- Subscripts switched from HTML `<sub>` tags to Pandoc `~subscript~` notation
  (`~k~`, `~z~`, `~95,k~`) for correct rendering in RStudio Visual editor and
  Quarto/HTML output.
- Methods illustration figure (Figure 1) integrated into Section 2 after Step 4,
  with full panel-by-panel caption; downstream figures renumbered 2–4.
- Figure 1 caption updated to note x-axis shows only the most recent 3,000 years
  of the CH10 record (full record ~6,200 cal yr BP).
- rbacon cited as Blaauw & Christen (2011); MCAgeDepth cited as Higuera et al.
  (2009) in Step 1 prose.

**`plot_chronUncertainty_methods.R`** — x-axis crop:

- Added `x_crop` user setting (oldest age to display); defaults to `NA` (full
  record). CH10 override set to 3,000 cal yr BP via `switch(out$site, ...)`.
- `ref_peak_ages_vis` filters benchmark peak vlines to the visible window so
  out-of-range peaks do not add clutter.
- Added `ggsave()` saving to `{site}_chronUncertainty_methods.png`
  (6 × 10 in, 300 dpi).

AI-assisted development (Claude/Anthropic).

---

### Per-zone mFRI, zone indicators, vignette restructure, age-depth figure (2026-06-23)

**`run_ensemble_analysis.R`** — per-zone mFRI floor:

- Replaced global `mean_mFRI` / `mFRI_floor` with a per-zone implementation.
  `mFRI_floor_by_zone[]` computes the mFRI/2 floor for each zone from the
  benchmark run's Weibull mean FRI. `mFRI_floor_for_age(age)` looks up the
  correct floor for any peak age, falling back to the nearest zone boundary
  for ages outside the record. `mFRI_floor_k` and `match_halfwin_k` are now
  vectors aligned to `ref_peak_ages`, with per-zone values.
- `zone_div` and `n_zones` are now defined once in Section 4 (removed
  duplicate definition that previously appeared in Section 9).
- Single-zone records (e.g. CH10) are unaffected; multi-zone records
  (e.g. SI17) now use zone-specific mFRI floors.

**`plot_ensemble_figure.R`** — zone boundary indicators:

- Added dashed grey vertical lines (`inner_zones`, `zone_vlines`) at interior
  zone boundaries in panels (b), (c), and (d). Centered zone labels added to
  panel (b). Single-zone records produce `NULL` for `zone_vlines`; ggplot
  ignores it silently.

**`plot_chronUncertainty_methods.R`** (renamed from `plot_methods_figure.R`):

- **Ensemble-only peaks** added to the bottom panel (Option B): grey open
  circles on the same y-axis, plotted via `cluster_unmatched()` with
  `ens_only_min_freq` threshold (default 40%); horizontal 95% CI bars shown
  where ≥2 iterations contribute.
- **Per-iteration deduplication** in `cluster_unmatched()`: detection
  frequency now counts unique iterations, not total detections, preventing
  double-counting when one iteration contributes two nearby unmatched peaks.
- **Zone boundary indicators** on all panels (same `zone_vlines` pattern as
  `plot_ensemble_figure.R`).
- **CH10-specific CHAR y-axis ceiling** via `switch(out$site, "CH10" = 50, NA)`
  in the user settings block; prevents extremely high outlier values from
  compressing the visible CHAR range.
- Bug fixes: `inherit.aes = FALSE` on `geom_errorbarh`/`geom_point` for
  ensemble-only peaks; `if/else` parsing fixed with explicit braces;
  `print(wrap_plots(...))` added so the figure displays when sourced from
  RStudio.

**Vignette** `vignettes/chronological_uncertainty_draft.md`:

- **Step 3** gained a full adaptive matching window description: formula
  (window_k = max(mFRI_z/2, CI_95,k/2)), ecological rationale for the
  mFRI_z/2 floor ("once a peak shifts more than half the mean FRI, it is
  likely a different event"), and a note that CI_95,k/2 typically dominates
  in records with meaningful age uncertainty.
- **Step 4**: "reference run" corrected to "benchmark run" throughout.
- **Sections 3–4** merged into a unified "Site examples" section (Section 3)
  with a framing paragraph, age-depth comparison figure (Figure 1), and
  summary table. CH10 is described as "unusually well-dated"; SI17 as "less
  precisely dated" (not "low precision"). Subsections 3.1/3.2 follow parallel
  structure: chronology → benchmark run → ensemble behavior and key takeaway.
  Figures renumbered 1–3.

**New: `tests/plot_agedepth_vignette.R`** — age-depth comparison figure:

- Plots median chronology and 95% CI ribbon for CH10 (blue) and SI17 (black)
  on a single panel. Age (cal yr BP) on reversed y-axis (present at top);
  depth on reversed x-axis. Shared axis limits allow direct visual comparison
  of sediment accumulation rates (line slopes) and chronological precision
  (ribbon widths). Legend positioned inside the panel at lower right.
  Saves to `tests/agedepth_vignette.png` (3.5 × 3.5 in, 300 dpi).

AI-assisted development (Claude/Anthropic).

---

### Secondary peak detection and coherence filter (2026-06-22)

**`run_ensemble_analysis.R`** — major update to ensemble-only peak detection:

- **Three peak types** now formally distinguished: *reference* (detected in the
  median-chronology run), *ensemble-only: near-reference* (secondary detections
  within a reference peak's matching window, absent from the reference run), and
  *ensemble-only: independent* (detections outside all reference windows).
- **Section 7 (time-step orphan scan)** now retains independent orphans only;
  near-reference cases are handled exclusively by the new Section 7b.
- **Section 7b (window-level secondary detection)**: for each reference peak *k*,
  scans all iterations for unclaimed detections within ±`match_halfwin_k`. Per
  iteration, the detection furthest from the reference peak age is recorded.
  Candidates present in ≥10% of iterations are flagged.
- **Merging pass**: overlapping secondary candidates (median ages within
  `mFRI_floor` of each other) are collapsed into a single event before summary
  statistics are computed.
- **Coherence filter**: secondary candidates with `ci_frac ≥ 0.5` are discarded
  as scatter. `ci_frac = ci95_width / (2 × match_halfwin)`; values near 1.0
  indicate detections scattered across the full window (noise); values < 0.5
  indicate genuine temporal clustering (real secondary event).
- Section 8 now combines reference peaks, independent orphans, and
  near-reference secondaries into a single `all_summaries` / `out_df`.
- Summary section (e) reports near-reference and independent ensemble-only
  peaks separately with distinct descriptions.

**`plot_ensemble_figure.R`** — updated for three peak types:

- Near-reference ensemble-only peaks now plotted with ±`match_halfwin` as the
  x-axis bar (honest uncertainty bound) rather than the 95% CI (which reflects
  the coherence filter criterion, not chronological uncertainty).
- Legend labels include CI type: "Reference peak (+/- 95% CI)",
  "Ensemble-only: near reference (+/- matching window)",
  "Ensemble-only: independent (+/- 95% CI)".
- Panel title simplified; `secondary_summaries` NULL-checked at load.

**New: `char_extract_bacon_chronologies.R`** — generic function for extracting
*N* chronologies from a completed rbacon MCMC run at user-supplied depths.
Samples one index vector applied consistently across all depths, preserving the
spatial correlation structure of the age-depth model. Returns a data frame
compatible with `char_run_ensemble.R`. Roxygen stubs included; destined for
`R/` with `rbacon` in Suggests.

**New: SI17 (Silver Lake) site files** — `SI17.csv` (Bacon input), 
`SI17_extract_bacon_chronologies.R` (site-specific extraction script),
`SI17_bacon_compare.R` (median age comparison; max diff 25.8 yr, RMSE 6.3 yr),
`SI17_MCAgeDepth_1000_chronologies.csv` (1000-member ensemble),
`SI17_peakAgeUncertainty.csv`, and `SI17_ensemble_figure.pdf/.png`.

**Applied results:**

- *SI17 (Silver Lake)*: 25 reference peaks, all detected in ≥88.2% of
  iterations (median 99.0%); median timing SD 53 yr. Zero ensemble-only peaks —
  wide matching windows (±105–280 yr) and many extra detections per iteration
  (median 32 vs. reference 25) produce only incoherent scatter (all candidates
  ci_frac ≥ 0.76). Correct null result for a high-uncertainty chronology.
- *CH10 (Chickaree Lake)*: 10 ensemble-only peaks total — 6 near-reference
  (ci_frac 0.19–0.29 for the tightest candidates) and 4 independent. High-
  precision chronology (narrow ±52 yr windows) allows genuine secondary
  clustering to be detected.

**Draft vignette** `vignettes/chronological_uncertainty_draft.md` Sections 1–4:
motivation (two-mechanism argument, post-hoc correction critique, individual
peak timing vs. fire-regime characteristics), workflow (4 steps with code
examples), and site examples for CH10 and SI17 with embedded figure references.

**`.gitignore`**: `Bacon_runs/` added (rbacon output directory).

---

## Development scripts (ROADMAP Section 3 — Chronological Uncertainty)

Three exploratory scripts in `tests/` (not yet exported functions), substantially
updated in the 2026-06-19 session:

**`char_run_ensemble.R`** — runs the CharAnalysis pipeline across an N-member
age-depth chronology ensemble and orchestrates the full workflow:

- **Parallel processing** via `parallel::parLapply` (socket clusters; works on
  Windows and Mac). Default: `detectCores() - 1` cores. Benchmark: ~15 min on
  1 core, ~2 min on 19 cores (20-core Windows laptop). Progress is reported
  every 100 iterations in minutes; a startup message describes the speedup and
  how to adjust `n_cores`.
- **`save_results` switch** (0 = interactive prompt, 1 = auto-save) controls
  saving of the ensemble RDS, the peak age uncertainty CSV, and the figure PDF
  together via a single flag — suitable for interactive use or batch/multi-site
  scripts.
- Automatically sources `run_ensemble_analysis.R` and `plot_ensemble_figure.R`
  at the end of each run. All output filenames derived from the site name in the
  params file (generic, not hardcoded to CH10).

**`run_ensemble_analysis.R`** — nearest-neighbor peak matching that computes
detection frequency and timing confidence intervals for each reference peak.
Matching half-window = max(mean Weibull mFRI / 2, `yrInterp`); a
reference-peak-first loop guarantees detection frequency ≤ 100%. Skips the
reference `CharAnalysis()` run and ensemble RDS load if those objects are
already in the workspace (e.g., when sourced from `char_run_ensemble.R`).

**`plot_ensemble_figure.R`** — three-panel ensemble figure (previously four panels):

- Panel (a): CHAR / C_background / peak ID from the reference run.
- Panel (b): SNI time series with the SNI = 3.0 advisory threshold.
- Panel (c, new): **detection frequency (%) on the y-axis** and **95% CI on
  peak timing as horizontal bars** on the x-axis (dot at reference peak age).
  Combines the information previously shown in separate beeswarm (c) and CI-bar
  (d) panels into a single display. Y-axis is data-adaptive: floored to the
  nearest 25% below the minimum observed detection frequency, with a dashed
  reference line at 50%.

Column name normalization: if `peak_summaries` is loaded from the CSV
(charResults-style names), the script renames columns to R-friendly names
automatically, so the fast-path workflow (load CSV, source script) requires
no manual renaming.

**New output file: `{site}_peakAgeUncertainty.csv`** (replaces `CH10_peak_summaries.csv`):
charResults-style column names with units in parentheses; detection frequency
stored as % (1 decimal); mean and median ages rounded to the nearest year;
`n_iter` stored as a column so the file is self-contained (no RDS needed to
regenerate the figure); no row-number column.

Applied to CH10 (Chickaree Lake): 59 reference peaks, detection frequency
58–100% (80% of peaks ≥ 90%), typical 95% CI widths 40–80 yr.

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
