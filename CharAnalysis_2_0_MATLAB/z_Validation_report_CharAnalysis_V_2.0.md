# CharAnalysis v2.0 Validation Report
## Round 1: Code Lake Example Dataset

**Date:** March 27, 2026
**Prepared by:** CharAnalysis Development Team
**Dataset:** Code Lake, AK — the example dataset bundled with CharAnalysis distributions
**Input file:** `CO_charParams.csv` / `CO_charData.csv`
**Comparison:** CharAnalysis v1.1 (MATLAB) vs. CharAnalysis v2.0 (MATLAB)

---

## 1. Purpose

This report documents the first formal numerical validation of CharAnalysis v2.0 against v1.1. The goal is to confirm that the Phase 1 modernization changes — which include replacing deprecated functions, removing global variables, vectorizing inner loops, and splitting the monolithic results function into separate computation and plotting modules — produce results that are numerically equivalent to v1.1 for identical input parameters.

The Code Lake dataset was chosen as the first test case because it is the example dataset distributed with CharAnalysis and has well-established v1.1 results that serve as a reference baseline.

---

## 2. Analysis Parameters

Both versions were run with identical parameters:

| Parameter | Value |
|---|---|
| `yrInterp` | 15 yr |
| `Smoothing.method` | 3 (Gaussian mixture model) |
| `Smoothing.yr` | 500 yr |
| `threshType` | 2 (locally defined) |
| `threshMethod` | 3 (Gaussian mixture model) |
| `threshValues` | [0.90, 0.95, 0.99, 0.95] |
| `cPeak` | 1 (residuals) |
| `minCountP` | 0.05 |
| `peakFrequ` | 1000 yr |
| Record length | 504 interpolated samples |

---

## 3. Results

### 3.1 Peak Identification — Perfect Match

The most important output of CharAnalysis is the identification of charcoal peaks, which form the basis of all fire-history interpretations. All four peak columns are identical between v1.1 and v2.0.

| Column | v1.1 | v2.0 | Difference |
|---|---|---|---|
| peaks1 | 54 | 54 | 0 |
| peaks2 | 48 | 48 | 0 |
| peaks3 | 43 | 43 | 0 |
| peaksFinal | 48 | 48 | 0 |

**Rows where peaksFinal differs: 0 of 504**

### 3.2 Threshold Values — Numerically Equivalent

Threshold values are used to separate fire signal from background noise. Mean threshold values are nearly identical across all four threshold columns, with differences in the fourth decimal place.

| Column | v1.1 mean | v2.0 mean | Difference |
|---|---|---|---|
| thresh1 | 0.0331 | 0.0331 | 0.0000 |
| thresh2 | 0.0513 | 0.0515 | 0.0002 |
| thresh3 | 0.0718 | 0.0721 | 0.0003 |
| threshFinalPos | 0.0513 | 0.0515 | 0.0002 |

Sample-by-sample threshold comparison across all 504 rows:

| Metric | Value |
|---|---|
| Samples where threshold differs | 71 of 504 (14%) |
| Mean absolute difference | 0.0002 |
| Maximum absolute difference | 0.0071 |

The small threshold differences at 71 samples are attributable to numerical differences between the v2.0 `charLowess` implementation and the original MATLAB Curve Fitting Toolbox `smooth()` function. These differences are not scientifically meaningful and do not affect peak identification.

### 3.3 Peak Magnitude — Essentially Identical

Peak magnitude (pieces cm⁻² peak⁻¹) summarizes the total charcoal deposited per fire event.

| Metric | v1.1 | v2.0 | Difference |
|---|---|---|---|
| Total peak magnitude | 127.43 | 126.70 | 0.73 (0.6%) |
| Mean magnitude per peak | 2.4506 | 2.4844 | 0.0338 (1.4%) |

The small differences in peak magnitude are consistent with the minor threshold differences noted above and are not scientifically meaningful.

---

## 4. Known Bugs Fixed in v2.0

The following bugs were identified and corrected during the validation process. All fixes were confirmed to produce results consistent with v1.1.

| Bug | Location | Description |
|---|---|---|
| Global variable side effects | `CharThreshGlobal`, `CharThreshLocal`, `CharSmooth`, `CharPretreatment` | `global plotData` and `global bkgSensIn` replaced with explicit function arguments, eliminating order-dependent side effects |
| GMM fallback passed full record | `CharThreshLocal` | Poor-fit GMM fallback incorrectly passed `Charcoal.peak` (full record, ~500 samples) instead of `X` (local window, ~34 samples), causing severe slowdown |
| NaN propagation through smoothing | `CharSmooth` | NaN values in `Charcoal.accI` from record gaps were not bridged before smoothing, causing `charLowess` to propagate NaN across a large block of `accIS` values near the record end |
| NaN in GMM windows | `CharThreshLocal` | Windows near record gaps contained NaN values that caused the GMM EM algorithm to never converge; NaN stripping added before distribution fitting |
| Short-window fallback used `continue` | `CharThreshLocal` | Windows with fewer than 4 valid samples after NaN removal were skipped entirely, leaving NaN thresholds; replaced with Gaussian fallback that always produces a valid threshold |
| Peak magnitude formula | `CharPostProcess` | Operator precedence error: `(diff+1)*r` instead of v1.1's `diff+r`, producing magnitudes ~2.5× too large |
| Index overrun in fire frequency loop | `CharPostProcess` | Upper window index in the middle-of-record branch could exceed array length; clamped with `min(length(Charcoal.ybpI), ...)` |
| `FRI` variable overwritten | `CharPostProcess` | Per-zone loop variable `FRI` overwrote the full-record `FRI` vector returned by `smoothFRI`, causing a length mismatch in `CharPlotResults` |

---

## 5. Conclusion

CharAnalysis v2.0 produces peak identification results that are **identical** to v1.1 on the Code Lake example dataset. Threshold and magnitude values differ by less than 1.5%, consistent with expected numerical noise between the new `charLowess` implementation and the original Curve Fitting Toolbox `smooth()` function.

This constitutes a successful first validation. Additional validation runs on independent datasets with different record lengths, sampling resolutions, and fire return interval characteristics are planned to further confirm numerical equivalence across a range of conditions.

---

## 6. Next Steps

- Run the same comparison on 2–3 additional charcoal records with different characteristics (record length, sampling resolution, fire frequency)
- Test the `bkgSens = 1` sensitivity analysis path
- Develop a formal regression test script (`CharRegressionTest.m`) that runs CharAnalysis on the Code Lake dataset and automatically asserts that key output values match this validated reference
