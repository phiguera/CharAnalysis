# *CharAnalysis* Development Roadmap
*Last updated: June 2026*
---
This document describes planned future development of *CharAnalysis*. Items are
listed in approximate priority order. It is not yet determined whether additional
development beyond the R translation will occur in R only or in both R and MATLAB.
---
## 1. R Translation ✓ *Complete — v2.0 (April 2026)*
*CharAnalysis* v2.0 is a direct R translation of the MATLAB v2.0 implementation,
released to CRAN in April 2026. The package is in the
[experimental](https://lifecycle.r-lib.org/articles/stages.html#experimental)
lifecycle stage; the API may evolve as user feedback is incorporated.
The R package (`CharAnalysis_2_0_R/`) reproduces the full five-step analytical
pipeline and nine output figures organized into diagnostic (Figs 1-5) and
analytical (Figs 6-9) categories. Key design choices:
- **Quantitative equivalence**: validated against MATLAB v2.0 on four benchmark
  datasets (Code Lake, Chickaree Lake, Silver Lake, Raven Lake). Known numerical
  differences are documented in `NEWS.md` and `inst/z_Validation_report_R_vs_MATLAB_V_2.0.md`.
- **GMM implementation**: the Gaussian mixture model threshold uses a direct port
  of the MATLAB `GaussianMixture.m` EM algorithm rather than an existing R package,
  preserving numerical comparability with the reference implementation.
- **Figures**: nine publication-quality ggplot2 figures; `char_plot_diagnostic()` (Figs 1-5) and `char_plot_analysis()` (Figs 6-9) separate parameter-evaluation from fire-history interpretation.
- **API**: snake_case function names (`char_plot_peaks()`, `char_write_results()`,
  etc.) following R community conventions; the top-level entry point `CharAnalysis()`
  retains its original name for continuity.
Install from CRAN:
```r
install.packages("CharAnalysis")
```
For the latest in-development version, install from the `dev` branch on
GitHub:
```r
devtools::install_github("phiguera/CharAnalysis",
                         subdir = "CharAnalysis_2_0_R",
                         ref    = "dev")
```
Possible future comparisons with related R packages that were developed based on
*CharAnalysis* Version 1.1:
- `tapas`: https://github.com/wfinsinger/tapas
- `CharcoalFireReconstructionR`: https://github.com/rglueckler/CharcoalFireReconstructionR
---
## 2. Record-quality diagnostics and user guidance ✓ *Complete — v2.0.4 (development)*

Implemented in the `dev` branch (June 2026):

- **SNI advisory**: `CharAnalysis()` now emits a `cli`-styled advisory at the
  end of all console output when SNI falls below 3.0 (Kelly et al. 2011). On
  the local-threshold path it reports the percentage of samples below the
  threshold; on the global path it reports the record-wide SNI value. The
  advisory is rendered with a red bordered block and colour-coded bullets, with
  graceful plain-text fallback. It is placed in `CharAnalysis()` (not inside
  `char_thresh_local/global`) so it does not fire repeatedly during
  `char_bkg_sensitivity()` runs. `cli` added to `Imports` in DESCRIPTION.
- **Vignette guidance**: a "Recommended workflow" section (Step 2) now
  instructs analysts to inspect `char_plot_sni()` (Fig 4) before proceeding to
  any fire-history figures. The rationale — peak detection is only appropriate
  where SNI > 3.0 — is stated explicitly, with guidance on adjusting smoothing
  windows or restricting interpretation when SNI is low.

---
## 3. Chronological Uncertainty — *In progress (development, June 2026)*

Incorporate methods for propagating chronological uncertainty into the
characterization of fire events. Development will take into account existing
approaches formalized in the following R packages, both of which were developed
based on *CharAnalysis* Version 1.1, and will include communication with those
developers:
- `tapas`: https://github.com/wfinsinger/tapas
- `CharcoalFireReconstructionR`: https://github.com/rglueckler/CharcoalFireReconstructionR

### Completed (June 2026, development scripts in `tests/`)

**`char_run_ensemble.R`** — runs the full CharAnalysis pipeline across a
user-supplied N-member age-depth chronology ensemble (default N = 1000).
Each iteration substitutes one ensemble chronology into the pipeline and
records the peak-detection result on a common time grid. Outputs a
`CharEnsemble` list with `peaks_matrix` (time steps × iterations),
`charPeak_matrix`, `charBkg_matrix`, `prob_peak`, and associated metadata.
Results are saved as an RDS for use across sessions. The bottleneck is the
GMM threshold step; 1000 iterations take approximately 16 minutes on a
standard laptop.

**`run_ensemble_analysis.R`** — nearest-neighbor peak matching that assigns
ensemble detections to reference peaks. For each reference peak *k* and each
iteration *i*, the single closest detected age within a half-window of
max(mean Weibull mFRI / 2, `yrInterp`) is assigned to *k*. The
reference-peak-first loop guarantees at most one detection per peak per
iteration, so P(peak) = n_detected / N ≤ 1.0. Outputs `CH10_peak_summaries.csv`
with P(peak), 95% CI on peak timing, n_detected, mean and median detected age,
and ±SD for each reference peak.

**`plot_ensemble_figure.R`** — four-panel ensemble figure. Panels (c) and (d)
are the new analytical panels:

- **Panel (c): P(peak) beeswarm.** Each reference peak is shown as one circle
  at its reference age (x-axis), sized by P(peak) using a universal fixed
  5-bin scale (0–20, 21–40, 41–60, 61–80, 81–100% of iterations). Circles
  are jittered vertically with a greedy beeswarm algorithm to prevent overlap.
  P(peak) quantifies detection robustness: a peak identified in ≥80% of
  ensemble iterations is robust to chronological uncertainty; one identified
  in <60% of iterations may be an artifact of the single best-estimate
  chronology and warrants caution in interpretation.

- **Panel (d): 95% CI on peak timing.** For each reference peak, the 2.5th
  and 97.5th percentiles of detected ages across the ensemble define the
  95% confidence interval on absolute timing. A dot marks the reference
  peak age. Bars are vertically jittered to prevent overlap. Typical CI
  widths for CH10 are 40–80 yr, reflecting a well-constrained age model.

Applied to the CH10 (Chickaree Lake) record: 59 reference peaks identified;
P(peak) range 0.584–1.000; 80% of peaks have P ≥ 0.90. The lowest-probability
cluster (2420–2500 cal yr BP, P = 0.58–0.80) coincides with three closely
spaced events, where chronological uncertainty causes detections to shift
among peaks.

### Also completed (June 2026, Round 8)

- **Per-zone mFRI floor**: `run_ensemble_analysis.R` now uses `mFRI_floor_by_zone[]`
  and `mFRI_floor_for_age()` to apply zone-specific mFRI floors to matching
  windows, replacing the global mean. Single-zone records unaffected; multi-zone
  records use zone-local fire frequency.
- **Zone boundary indicators**: dashed grey vertical lines and centered zone
  labels added to all ensemble figure panels and the methods illustration figure.
- **Methods illustration figure** (`plot_chronUncertainty_methods.R`, renamed
  from `plot_methods_figure.R`): ensemble-only peaks added to bottom panel with
  per-iteration deduplication; CH10-specific y-axis ceiling; multiple bug fixes.
- **Adaptive matching window documentation**: vignette Step 3 now describes the
  window formula, ecological rationale for the mFRI/2 floor, and the observation
  that CI_95,k/2 typically dominates in records with meaningful age uncertainty.
- **Vignette restructured**: CH10 and SI17 examples merged into a unified "Site
  examples" section (Section 3) with age-depth comparison figure, summary table,
  and parallel subsections. CH10 described as "unusually well-dated"; SI17 as
  "less precisely dated." All "reference run" → "benchmark run."
- **Age-depth comparison figure** (`tests/plot_agedepth_vignette.R`): single-panel
  figure with shared axes showing both sites' median chronology and 95% CI ribbon;
  highlights differences in accumulation rate and chronological precision.

### Also completed (June 2026, Round 7)

- **Three peak types**: reference, ensemble-only near-reference, and
  ensemble-only independent are now formally distinguished in analysis,
  output CSV, and figure. Near-reference peaks are detected via a window-level
  secondary scan (Section 7b) rather than the time-step orphan scan, which is
  now reserved for independent peaks only.
- **Window-level secondary detection**: for each reference peak, unclaimed
  detections across iterations are aggregated within ±match_halfwin. A merging
  pass (mFRI_floor criterion) collapses overlapping candidates; a coherence
  filter (ci_frac < 0.5) suppresses scatter and retains only temporally
  coherent secondary events.
- **rbacon interface**: `char_extract_bacon_chronologies()` extracts N
  chronologies from a completed rbacon run, preserving MCMC spatial correlation.
  Roxygen-documented; destined for `R/` with `rbacon` in Suggests.
- **SI17 applied**: null result (0 ensemble-only peaks) confirmed correct for
  a wide-uncertainty chronology; CH10 yields 10 ensemble-only peaks (6
  near-reference, 4 independent).
- **Draft vignette**: `vignettes/chronological_uncertainty_draft.md` Sections
  1–4 written (motivation, workflow, CH10 example, SI17 example).
- **Parallel processing**: implemented in Round 6 (2026-06-22); ✓ complete.

### Pending

- Test ensemble pipeline on additional sites beyond CH10 and SI17; CH10's
  high-precision chronology may make it unusual in yielding coherent
  near-reference secondary peaks.
- Package `char_run_ensemble()` and `char_extract_bacon_chronologies()` as
  exported functions with full Roxygen documentation and input validation.
- Add `rbacon` to DESCRIPTION Suggests; move `char_extract_bacon_chronologies()`
  to `R/`.
- Complete vignette (`vignettes/chronological_uncertainty_draft.md` →
  `chronological_uncertainty.Rmd`): convert to Rmd, add rendered figures,
  finalize for package build.
- Coordinate with `tapas` and `CharcoalFireReconstructionR` developers.
- Validate global-threshold path (`threshType = 1`) against MATLAB.
---
## 4. Regional Synthesis

Add support for synthesizing peak identification across multiple sediment-charcoal
records at regional scales. This will use and generalize methods already developed
and applied in the following publications:
Higuera, P.E., B.N. Shuman, and K.D. Wolf. 2021. Rocky Mountain subalpine forests
now burning more than any time in recent millennia. *Proceedings of the National
Academy of Sciences* 118:e2103135118.
https://www.pnas.org/doi/abs/10.1073/pnas.2103135118
Clark-Wolf, K.D., P.E. Higuera, B.N. Shuman, and K.K. McLauchlan. 2023. Wildfire
activity in northern Rocky Mountain subalpine forests still within millennial-scale
range of variability. *Environmental Research Letters* 18:094029.
https://doi.org/10.1088/1748-9326/acee16

---
*This roadmap reflects current intentions and is subject to change. Feedback and
collaboration are welcome — please use the Issues tab at
https://github.com/phiguera/CharAnalysis/issues to share ideas or express interest
in contributing.*

*Development of Version 2.0 and planning for future updates is being carried out with
the assistance of Claude, an AI assistant by Anthropic. All work has been reviewed
and validated by the author.*
