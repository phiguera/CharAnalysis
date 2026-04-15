# *CharAnalysis* 2.0.0 (beta)

First R implementation of *CharAnalysis*, a direct translation of
*CharAnalysis* v2.0 (MATLAB). Analytical outputs are validated against
four benchmark datasets; user testing is ongoing.

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

**Figures 9 and 10.** The threshold-sensitivity detail plot (Fig. 9) and
multi-site comparison plot (Fig. 10) from the MATLAB interface are not
implemented in this R package.

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
