# CharAnalysis R Translation — Project Brief
*Prepared for Claude Cowork | April 2026 — Updated April 2026 to reflect as-built implementation*

## 1. Project Overview

CharAnalysis is a MATLAB program for reconstructing local fire histories from lake-sediment charcoal records. It implements a peak-detection workflow that decomposes a charcoal accumulation rate (CHAR) time series into low-frequency background and high-frequency peak components, identifies charcoal peaks that likely represent fire events, and summarizes fire-history metrics.

The program has been in use since the mid-2000s and has been used in dozens of published studies analyzing sediment-charcoal records worldwide. The GitHub repository is https://github.com/phiguera/CharAnalysis.

Version 2.0 (March 2026) is the current MATLAB codebase. It is a fully validated modernization of the original code, with no changes to analytical methods. It runs on MATLAB R2019a or higher with no toolbox dependencies. All MATLAB source files are in the `CharAnalysis_2_0_MATLAB/` subfolder of the repository.

The goal of this project is to translate CharAnalysis Version 2.0 from MATLAB into R, producing a well-structured R package that is quantitatively equivalent to the MATLAB implementation on validated reference datasets.

R is the dominant language in the paleoecological community today. An R implementation will significantly broaden access to CharAnalysis for new users and make integration with downstream analyses much easier.

---

## 2. MATLAB Codebase Reference

The MATLAB codebase is organized as a set of `.m` function files called by a single entry point. The analytical call hierarchy is:

```
CharAnalysis.m            (entry point)
  ├── CharParameters.m    Read .csv / .xls input file
  ├── CharPretreatment.m  Interpolate, compute CHAR, transform
  ├── CharSmooth.m        Estimate Cbackground (5 methods)
  ├── CharThreshGlobal.m  Global threshold + noise PDF
  ├── CharThreshLocal.m   Local (sliding-window) threshold
  ├── CharPeakID.m        Peak flagging + minimum-count screen
  ├── CharPostProcess.m   Post-processing: FRIs, Weibull, SNI, GOF
  └── CharPlotResults.m   All output figures (Figures 3–9)

Supporting:
  charLowess.m            Base-MATLAB lowess implementation (no toolbox)
  CharValidateParams.m    Input validation
  GaussianMixture.m       EM algorithm entry point
  EStep.m, MStep.m, EMIterate.m, MDLReduceOrder.m,
  ClusterNormalize.m, initMixture.m, SplitClasses.m,
  GMClassLikelihood.m     EM sub-functions (9 files)
  smoothFRI.m             Smoothed FRI and fire-frequency curves
  bkgCharSensitivity.m    Sensitivity to Cbackground window width
  CharPlotFig3–9.m        Individual figure functions (modular interface)
```

### Key data structures

Data flows through three primary structs and several parameter structs:

| Struct | Contents |
|--------|----------|
| `Charcoal` | Raw and resampled time series: depths, ages, counts, volumes, concentrations, CHAR, smoothed CHAR, peak CHAR, identified peaks, magnitudes, fire frequencies |
| `CharThresh` | Threshold values (pos/neg, 4 columns), noise PDF, SNI, GOF, minCountP |
| `Pretreatment` | zoneDiv, yrInterp, transform |
| `Smoothing` | method (1–5), yr (window width) |
| `PeakAnalysis` | cPeak, threshType, threshMethod, threshValues[4], minCountP, peakFrequ, sensitivity |
| `Results` | save, saveFigures, allFigures |

---

## 3. Reference Datasets

Five dataset/configuration combinations have been validated between MATLAB V1.1 and V2.0 and serve as the quantitative benchmark for the R translation. All datasets are included in the repository under `DataTemplates_and_Examples/`.

| Dataset | Configuration | Notes |
|---------|---------------|-------|
| Code Lake | Local threshold (default params) | Primary test case |
| Code Lake | Global threshold | Tests the global threshold path |
| Chickaree Lake | CH10 (non-default tolerances) | GMM variability near record gaps |
| Thunder Lake | TL06 | — |
| Raven Lake | RA07 | — |

These five combinations are the required validation targets for the R translation. At each phase, R outputs must be compared numerically against the corresponding MATLAB V2.0 CSV output files.

---

## 4. Analytical Workflow Summary

The five steps the program performs are listed here for reference when implementing R equivalents. See the User's Guide (`CharAnalysis_UsersGuide.md`) for the complete description.

1. **Interpolate** — Resample the raw charcoal series (concentration, sediment accumulation rate) to equal time intervals using a proportion-weighted algorithm that preserves sediment structure. Compute CHAR (pieces cm⁻² yr⁻¹). Optionally apply log transformation.

2. **Smooth** — Model low-frequency trends in CHAR (C_background) using one of five methods: lowess, robust lowess, moving average, moving median, or moving mode.

3. **Remove background** — Compute the high-frequency component C_peak as either residuals (C_int − C_back) or ratios (C_int / C_back).

4. **Threshold** — Define a threshold value t and flag samples where C_peak > t as candidate peaks. Threshold can be globally defined (whole record) or locally defined (sliding window). Two noise-distribution approaches: Gaussian assumption or Gaussian mixture model.

5. **Screen peaks** — Apply a minimum-count test (two-sample Poisson test) to remove statistically insignificant peaks. Compute fire-history metrics: peak magnitudes, fire return intervals (FRIs), Weibull fits, bootstrapped CIs, smoothed fire frequency, SNI, and KS goodness-of-fit.

---

## 5. Phased Implementation Plan

The translation is organized into four sequential phases. Phases 1–3 are the core deliverables. Phase 4 is packaging and polish.

### Phase 1 — Parameter Reading and Pretreatment

**Goal:** Read input files and reproduce the interpolated CHAR series exactly.

R functions to implement:

- `char_parameters()` — Read the `.csv` parameter file and charcoal data file. Replicate the parameter-validation logic in `CharValidateParams.m`. Return a named list equivalent to the MATLAB parameter structs.
- `char_pretreatment()` — Reproduce the proportion-weighted resampling algorithm from `CharPretreatment.m`. This is the most arithmetically sensitive function in the program: the nested loop that builds the proportion matrix must be vectorized (as it is in V2.0 MATLAB) and produce identical results. Optionally apply log transformation.

R packages needed: `readxl` or `readr` (input), base R (computation).

**Quantitative validation targets (Phase 1):**

Compare the following columns from the MATLAB V2.0 CSV output against R output on all five reference datasets. Tolerances are absolute differences:

| Column | Variable | Required tolerance |
|--------|----------|--------------------|
| `cmTop_i` | Interpolated depth | Exact match |
| `ageTop_i` | Interpolated age | Exact match |
| `charCount_i` | Interpolated count | < 1e-10 |
| `charVol_i` | Interpolated volume | < 1e-10 |
| `charCon_i` | Interpolated concentration | < 1e-10 |
| `charAcc_i` | Interpolated CHAR | < 1e-10 |

**Validation script (Phase 1):**

```r
# Example validation pattern — apply to each of the five reference datasets
matlab_out <- read.csv("CodeLake_V2_local_MATLAB.csv")
r_out      <- CharAnalysis("CodeLake_charParams.csv")

stopifnot(max(abs(r_out$ageTop_i - matlab_out$ageTop_i)) < 1e-10)
stopifnot(max(abs(r_out$charAcc_i - matlab_out$charAcc_i)) < 1e-10)
```

---

### Phase 2 — Smoothing and Threshold Determination

**Goal:** Reproduce C_background and C_peak, and determine threshold values.

R functions to implement:

- `char_smooth()` — Implement all five smoothing methods from `CharSmooth.m`. Methods 1–2 (lowess, robust lowess) require a custom wrapper that maps MATLAB's span convention (number of points or fraction of data length) to R's `lowess()` convention (fraction only, always in (0,1)). Methods 3–5 (moving average, moving median, moving mode) map to `zoo::rollmean()`, `zoo::rollmedian()`, and a custom `lapply + tabulate` approach.
- `char_thresh_global()` — Implement the global threshold path from `CharThreshGlobal.m`. Two sub-methods: Gaussian assumption (`threshMethod = 2`) and Gaussian mixture model (`threshMethod = 3`). For the GMM, use `mclust::Mclust(G = 2)` as a replacement for the bundled MATLAB EM implementation. See the known differences note in Section 6 below.
- `char_thresh_local()` — Implement the sliding-window local threshold from `CharThreshLocal.m`. This is the most computationally demanding function. Pre-compute window indices before the loop. Apply the same two noise-distribution sub-methods as the global case.

R packages needed: `zoo`, `mclust`, `stats` (base).

**Lowess wrapper note:** MATLAB's `charLowess.m` (V2.0) uses span as the number of points or as a fraction of data length. R's `lowess()` uses `f` as a fraction always. Write a thin wrapper `char_lowess(y, span)` that converts span to `f` before calling `lowess()`. Validate this wrapper carefully against MATLAB output on all five reference datasets before using it in `char_smooth()`.

**Quantitative validation targets (Phase 2):**

| Column | Variable | Required tolerance |
|--------|----------|--------------------|
| `charBkg` | C_background | < 1e-8 |
| `charPeak` | C_peak | < 1e-8 |
| `thresh1`–`thresh3`, `threshFinalPos`, `threshFinalNeg` | Threshold values | < 1e-6 |
| `SNI` | Signal-to-noise index | < 1e-6 |
| `threshGOF` | KS goodness-of-fit p-value | < 1e-4 |

For the Chickaree Lake (CH10) dataset with `threshMethod = 3` (GMM), the allowable tolerance for `threshFinalPos` is relaxed to < 0.015 due to irreducible GMM variability near record gaps (documented in the V2.0 validation report). This relaxed tolerance must be explicitly noted in the Phase 2 validation output.

---

### Phase 3 — Peak Identification and Fire-History Metrics

**Goal:** Reproduce all peak identification and post-processing outputs.

R functions to implement:

- `char_peak_id()` — Implement peak flagging and minimum-count screening from `CharPeakID.m`. Peak flagging is a vectorized comparison. The minimum-count test uses a two-sample Poisson statistic; at 10^10 degrees of freedom, the t-distribution is numerically identical to the standard normal, so use `pnorm()` directly.
- `char_post_process()` — Implement post-processing from `CharPostProcess.m`. Includes: peak magnitudes, smoothed FRI curve, smoothed fire-frequency curve, FRI distributions by zone, Weibull fits (use `MASS::fitdistr()` on raw FRI values rather than histogram-based fitting), bootstrapped 95% CIs (use `boot::boot()`), and two-sample KS tests between zones.

**Weibull parameterization note:** MATLAB uses parameter order (a = scale, b = shape); R's `MASS::fitdistr()` and `pweibull()` / `dweibull()` use (shape, scale). Swap consistently in all downstream calculations and document this explicitly in the code.

R packages needed: `MASS`, `boot`, `stats` (base).

**Quantitative validation targets (Phase 3):**

| Column | Variable | Required tolerance |
|--------|----------|--------------------|
| `peaks1`–`peaks3`, `peaksFinal` | Peak binary flags | Exact match |
| `peaksInsig` | Insignificant peaks | Exact match |
| `peakMag` | Peak magnitude | < 1e-6 |
| `smPeakFrequ` | Smoothed fire frequency | < 1e-6 |
| `smFRIs` | Smoothed FRI curve | < 1e-6 |
| `mFRI` | Mean FRI by zone | < 1e-3 |
| `WBLb`, `WBLc` | Weibull parameters | < 0.5 (see note) |

**Weibull CI note:** Bootstrapped Weibull CIs (`WBLb_uCI`, `WBLb_lCI`, `WBLc_uCI`, `WBLc_lCI`) will differ between MATLAB and R due to random seed differences in bootstrap resampling. Validate point estimates (`WBLb`, `WBLc`) only; treat CI comparisons as qualitative.

---

### Phase 4 — Output, Figures, and Packaging

**Goal:** Produce equivalent output figures and deliver a usable R package.

R functions to implement:

- `char_write_results()` — Write the output data table to CSV (columns matching the MATLAB `CharResults` output exactly, same column names and order).
- **Figure functions** — Build `ggplot2` equivalents for each of the nine standard output figures. Implement as standalone functions analogous to the MATLAB modular figure interface (`CharPlotFig3.m` through `CharPlotFig9.m`). Focus on information content, not pixel-level reproduction of MATLAB figure aesthetics. Use `patchwork` for multi-panel layouts and `ggplot2::sec_axis()` for dual y-axes.
- `CharAnalysis()` — Top-level wrapper that calls the full pipeline and returns a named list of all outputs (equivalent to the MATLAB `results` struct). Call `CharWriteResults()` separately to save the output CSV.

**Package structure (as built):**

```
R/
  CharParameters.R
  CharPretreatment.R
  CharSmooth.R
  CharThreshGlobal.R
  CharThreshLocal.R
  CharPeakID.R
  CharPostProcess.R
  CharWriteResults.R
  CharPlotResults.R       (all five figure functions + char_plot_all wrapper)
  CharAnalysis.R          (top-level pipeline wrapper)
  charLowess.R
  CharValidateParams.R
  GaussianMixture.R       (direct port of MATLAB EM; see Section 6.2)
tests/
  test_phase1.R
  test_phase2.R
  test_phase3.R
  test_phase4.R
vignettes/
  CharAnalysis_intro.Rmd
DESCRIPTION
NAMESPACE
```

R packages needed: `ggplot2`, `patchwork`, `ggtext`, `MASS`, `zoo`, `stats` (base).

**Validation at Phase 4:** Re-run the full Phase 1–3 validation suite through `char_run()` to confirm end-to-end equivalence has not regressed.

---

## 6. Known Translation Issues and Required Workarounds

These are the most important implementation decisions. Address each one explicitly and document the resolution in code comments.

### 6.1 Lowess span convention

MATLAB's `charLowess.m` (V2.0) accepts `span` as either a fraction of data length (if < 1) or a number of points (if ≥ 1). R's `lowess()` accepts `f` as a fraction only. Write a wrapper that converts: `f = span / length(y)` when `span >= 1`, and `f = span` when `span < 1`. Validate the wrapper on all five reference datasets before integrating it into `char_smooth()`.

### 6.2 Gaussian mixture model — direct port of MATLAB EM

The MATLAB codebase bundles a custom EM implementation with MDL-based order selection (`GaussianMixture.m` and eight sub-function files). Rather than replacing this with `mclust`, the R implementation directly ports the MATLAB EM algorithm in `GaussianMixture.R` (and companion files `EStep.R`, `MStep.R`, `EMIterate.R`, `MDLReduceOrder.R`, `ClusterNormalize.R`, `initMixture.R`, `SplitClasses.R`, `GMClassLikelihood.R`). This produces numerically close agreement with MATLAB on the CO reference dataset. Floating-point divergence in low-data windows (e.g., Chickaree Lake CH10) may still occur and is expected; the allowed tolerance for threshold values on CH10 is relaxed to ± 0.015 (see Phase 2 validation targets above). No `mclust` dependency is needed.

### 6.3 Weibull parameterization order

MATLAB's `wblfit()` returns (a = scale, b = shape). R's `MASS::fitdistr()` and distribution functions use (shape, scale). This swap must be applied consistently everywhere Weibull parameters are computed or used downstream.

### 6.4 Minimum-count test: tcdf vs. pnorm

MATLAB's `CharPeakID.m` uses `1 - tcdf(d, 1e10)` for the minimum-count p-value. With 10^10 degrees of freedom, the t-distribution is numerically identical to the standard normal. Use `1 - pnorm(d)` in R. This is mathematically equivalent and requires no tolerance relaxation.

### 6.5 Bootstrap CIs: bootstrp vs. boot

MATLAB's `bootstrp(1000, 'mean', FRI)` maps to `boot::boot(FRI, function(x, i) mean(x[i]), R = 1000)$t`. Bootstrapped CIs will differ between implementations due to random seed differences. Validate point estimates only; treat CI differences as expected.

### 6.6 Anderson-Darling test in smoothFRI

MATLAB calls a custom `AnDarksamtest()` bundled with the distribution. The R implementation uses `ks.test()` (base R) for two-sample comparisons between zones, which avoids an additional package dependency. If Anderson-Darling tests are needed in a future revision, `kSamples::ad.test(list(x1, x2), method = "asymptotic")` is the appropriate replacement.

---

## 7. Quantitative Validation Protocol

### 7.1 MATLAB reference output files

Before beginning R development, generate fresh MATLAB V2.0 CSV output for all five reference datasets and save them as the authoritative reference files. These files contain the ground truth for all numerical comparisons.

Suggested naming convention:

```
CodeLake_local_MATLAB_V2.csv
CodeLake_global_MATLAB_V2.csv
Chickaree_CH10_MATLAB_V2.csv
Thunder_TL06_MATLAB_V2.csv
Raven_RA07_MATLAB_V2.csv
```

### 7.2 Validation script structure

*Note: As of the initial dev-branch release, only the CO (Code Lake, local threshold) dataset has been fully validated through Phase 4. The validation script described below has not yet been written. Creating `tests/validate_all.R` and validating the remaining four reference datasets (CodeLake global, CH10, TL06, RA07) is a planned post-merge task.*

Create a single validation script (`tests/validate_all.R`) that runs all five reference datasets through the R implementation and compares against the MATLAB reference CSVs. The script should:

1. Load the MATLAB reference CSV and the R output for each dataset.
2. For each validated column, compute `max(abs(r_out - matlab_ref))`.
3. Assert the value is within the required tolerance for that column and phase.
4. Print a pass/fail summary table.

Example summary table format:

```
Dataset              Column         Max abs diff   Tolerance   Status
CodeLake (local)     charAcc_i      3.2e-13        1e-10       PASS
CodeLake (local)     charBkg        4.1e-9         1e-8        PASS
CodeLake (local)     threshFinalPos 2.3e-7         1e-6        PASS
CodeLake (local)     peaksFinal     0              exact       PASS
Chickaree (CH10)     threshFinalPos 0.0087         0.015       PASS
...
```

### 7.3 Regression test maintenance

Once Phase 3 validation passes, freeze the reference CSV files and the validation script as permanent regression tests. Any future change to the R code must pass the full validation suite before being merged.

---

## 8. Repository and Documentation

- **GitHub repository:** https://github.com/phiguera/CharAnalysis
- **MATLAB V2.0 source:** `CharAnalysis_2_0_MATLAB/`
- **Reference datasets:** `DataTemplates_and_Examples/`
- **User's Guide:** `CharAnalysis_UsersGuide.md`
- **Modernization and translation planning document:** `CharAnalysis_Modernization_Report.docx`

The R translation should be developed in the same repository, in a new subfolder `CharAnalysis_R/`, and initially on a dedicated branch (e.g. `r-translation`).

### Documentation expectations

- All R functions should have roxygen2 documentation headers.
- A vignette (`vignettes/CharAnalysis_intro.Rmd`) should demonstrate the full workflow on the Code Lake example dataset.
- The validation script and all reference CSV files should be committed to the repository alongside the R source.

---

## 9. R Package Dependencies

| Package | Status | Purpose |
|---------|--------|---------|
| `MASS` | Required | `fitdistr()` for Weibull parameter estimation |
| `zoo` | Required | `rollmean`, `rollmedian`, `rollapply` in `CharSmooth.R` |
| `stats` | Required (base) | `ks.test()`, `lowess()`, `qnorm`, `pnorm`, `dnorm` |
| `graphics` | Required (base) | Legacy plotting fallback |
| `utils` | Required (base) | `read.csv()`, `write.table()` |
| `ggplot2` (≥ 3.4.0) | Suggested | All output figures via `CharPlotResults.R` |
| `patchwork` | Suggested | Multi-panel figure layout |
| `ggtext` | Suggested | HTML/markdown subscripts in figure panel titles |
| `readxl` | Suggested | Read Excel parameter files (if `.xls`/`.xlsx` input) |

*Packages listed in earlier drafts (`mclust`, `boot`, `kSamples`, `openxlsx`) are **not** used in the as-built implementation.*

---

## 10. Summary of Deliverables by Phase

| Phase | Deliverable | Validation status |
|-------|-------------|-------------------|
| 1 | `CharParameters()`, `CharPretreatment()` | CO dataset ✓ — remaining 4 datasets pending |
| 2 | `CharSmooth()`, `CharThreshGlobal()`, `CharThreshLocal()`, `charLowess.R` | CO dataset ✓ — remaining 4 datasets pending |
| 3 | `CharPeakID()`, `CharPostProcess()` | CO dataset ✓ — remaining 4 datasets pending |
| 4 | `CharWriteResults()`, `CharPlotResults.R`, `CharAnalysis()`, vignette, DESCRIPTION, NAMESPACE | CO dataset ✓ — `validate_all.R` and full 5-dataset regression suite pending |

---

*This project was planned with the assistance of Claude, an AI assistant by Anthropic. All code will be reviewed and validated by Philip Higuera (philip.higuera@umontana.edu) against Version 2.0 MATLAB reference outputs before being merged to the main branch.*
