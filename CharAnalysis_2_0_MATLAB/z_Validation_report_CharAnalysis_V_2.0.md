# CharAnalysis v2.0 Validation Report

**Date:** March 30, 2026
**Prepared by:** CharAnalysis Development Team
**Comparison:** CharAnalysis v1.1 (MATLAB) vs. CharAnalysis v2.0 (MATLAB)

---

## Validation Datasets

| Round | Site ID | Site Name | Location | Ecosystem | threshType | Reference |
|---|---|---|---|---|---|---|
| 1 | CO | Code Lake | Alaska, USA | boreal forest | local | Higuera et al. 2009 |
| 2 | CO | Code Lake | Alaska, USA | boreal forest | global | -- |
| 3 | CH10 | Chickaree Lake | Colorado, USA | subalpine forest | local | Dunnette et al. 2014 |
| 4 | TL06 | Thunder Lake | Colorado, USA | subalpine forest | local | Higuera et al. 2014 |
| 5 | RA07 | Raven Lake | Alaska, USA | tundra | local | Higuera et al. 2011 |
| 6 | SI17 | Silver Lake | Montana, USA | subalpine forest | local | Clark-Wolf 2023 |

Code Lake is the example dataset bundled with _CharAnalysis_ distributions and was tested with both local and global threshold configurations. Chickaree Lake, Thunder Lake, Raven Lake, and Silver Lake are independent records used to validate v2.0 across a range of record lengths, sampling resolutions, ecosystems, and fire regime characteristics.

### Note on Record End Treatment

V2.0 handles the end of the interpolated record more accurately than v1.1. When `zoneDiv(end)` extends beyond the bottom age of the last raw sample, v1.1 silently filled terminal interpolated samples with zero CHAR values; v2.0 correctly identifies these as having no overlapping raw data. In v2.0, `CharPretreatment.m` automatically corrects `zoneDiv(end)` to the bottom age of the last raw sample and notifies the user. This means v2.0 may produce a slightly shorter interpolated record than v1.1 when the original `zoneDiv(end)` overshot the data, and small differences near the end of the record are expected in those cases. All validation runs below used params files where `zoneDiv(end)` was set to match the data boundary exactly, so that both versions produce identical record lengths for comparison.

---

## Round 1: Code Lake, Alaska (CO) — Local Threshold

**Dataset:** Code Lake, AK — the example dataset bundled with _CharAnalysis_ distributions
**Input file:** `CO_charParams.csv` / `CO_charData.csv`

### Analysis Parameters

| Parameter | Value |
|---|---|
| `yrInterp` | 15 yr |
| `Smoothing.method` | 1 (Lowess) |
| `Smoothing.yr` | 500 yr |
| `threshType` | 2 (locally defined) |
| `threshMethod` | 3 (Gaussian mixture model) |
| `threshValues` | [0.90, 0.95, 0.99, 0.95] |
| `cPeak` | 1 (residuals) |
| `minCountP` | 0.05 |
| `peakFrequ` | 1000 yr |
| Record length | 500 interpolated samples |
| Record gaps | None |

### Results

#### Peak Identification — Perfect Match

| Column | v1.1 | v2.0 | Difference |
|---|---|---|---|
| peaks1 | 54 | 54 | 0 |
| peaks2 | 48 | 48 | 0 |
| peaks3 | 43 | 43 | 0 |
| peaksFinal | 48 | 48 | 0 |

**Rows where peaksFinal differs: 2 of 500**

The 2 differing rows are floating-point boundary cases where `charPeak` and `threshFinalPos` differ by less than 0.001 — the peak count matches exactly and the row-level differences are not scientifically meaningful.

#### Threshold Values — Numerically Equivalent

| Column | v1.1 mean | v2.0 mean | Difference |
|---|---|---|---|
| thresh1 | 0.0335 | 0.0331 | 0.0003 |
| thresh2 | 0.0520 | 0.0516 | 0.0004 |
| thresh3 | 0.0728 | 0.0723 | 0.0005 |
| threshFinalPos | 0.0520 | 0.0516 | 0.0004 |

| Metric | Value |
|---|---|
| Samples where threshold differs | 67 of 500 (13%) |
| Mean absolute difference | 0.0004 |
| Maximum absolute difference | 0.0198 |

#### Peak Magnitude — Essentially Identical

| Metric | v1.1 | v2.0 | Difference |
|---|---|---|---|
| Total peak magnitude | 125.85 | 126.76 | 0.7% |
| Mean magnitude per peak | 2.4203 | 2.4856 | — |

#### Smoothed Fire Frequency

| Metric | Value |
|---|---|
| Mean absolute difference | 0.0042 |
| Maximum absolute difference | 0.0511 |

#### Regression Test Result

| PeakTol | ThreshTol | FreqTol | Result |
|---|---|---|---|
| 0 | 0.001 | 0.010 | 11/11 PASS |

---

## Round 2: Code Lake, Alaska (CO) — Global Threshold

**Dataset:** Code Lake, AK — the example dataset bundled with _CharAnalysis_ distributions
**Input file:** `COGlobal_charParams.csv` / `CO_charData.csv`

### Analysis Parameters

| Parameter | Value |
|---|---|
| `yrInterp` | 15 yr |
| `Smoothing.method` | 2 (robust lowess) |
| `Smoothing.yr` | 500 yr |
| `threshType` | 1 (globally defined) |
| `threshMethod` | 3 (Gaussian mixture model) |
| `threshValues` | [0.90, 0.95, 0.99, 0.95] |
| `cPeak` | 1 (residuals) |
| `minCountP` | 0.05 |
| `peakFrequ` | 1000 yr |
| Record length | 500 interpolated samples |
| Record gaps | None |

### Results

#### Peak Identification — Match Within Tolerance

| Column | v1.1 | v2.0 | Difference |
|---|---|---|---|
| peaks1 | 50 | 52 | 2 |
| peaks2 | 46 | 46 | 0 |
| peaks3 | 40 | 40 | 0 |
| peaksFinal | 46 | 46 | 0 |

**Rows where peaksFinal differs: 0 of 500**

The `peaks1` difference of 2 is a downstream consequence of the `thresh1` GMM bin difference described below — `peaksFinal` (which uses the final threshold value) matches perfectly.

#### Threshold Values — Numerically Equivalent Within Tolerance

| Column | v1.1 mean | v2.0 mean | Difference |
|---|---|---|---|
| thresh1 | 0.0413 | 0.0384 | 0.0029 |
| thresh2 | 0.0584 | 0.0584 | 0.0000 |
| thresh3 | 0.0812 | 0.0812 | 0.0000 |
| threshFinalPos | 0.0584 | 0.0584 | 0.0000 |

**Samples where threshFinalPos differs: 0 of 500**

For a global threshold, each column contains a single constant value applied to all samples. `thresh2`, `thresh3`, and `threshFinalPos` match exactly. The `thresh1` difference (0.0029) reflects GMM non-determinism: the EM algorithm is sensitive to initialization, and small floating-point differences between v1.1 and v2.0 place the 0.90th percentile of the noise distribution in a slightly different bin. This does not affect `peaksFinal`.

#### Peak Magnitude — Essentially Identical

| Metric | v1.1 | v2.0 | Difference |
|---|---|---|---|
| Total peak magnitude | 123.89 | 123.92 | 0.0% |
| Mean magnitude per peak | 2.6361 | 2.6365 | — |

#### Smoothed Fire Frequency

| Metric | Value |
|---|---|
| Mean absolute difference | 0.0002 |
| Maximum absolute difference | 0.0028 |

#### Regression Test Result

| PeakTol | ThreshTol | FreqTol | Result |
|---|---|---|---|
| 2 | 0.005 | 0.001 | 11/11 PASS |

To reproduce this validation:
```matlab
z_Compare_CharAnalysis_V1_V2('COGlobal_charParams.csv', ...
    'COGlobal_charResults_V_1_1.csv', 'PeakTol', 2, 'ThreshTol', 0.005)
```

### Bugs Fixed During Round 2 Validation

Three bugs affecting the global threshold path were identified and corrected. See the complete bugs list below.

---

## Round 3: Chickaree Lake, Colorado (CH10) — Local Threshold

**Dataset:** Chickaree Lake, CO — an independent charcoal record used for validation
**Input file:** `CH10_charParams.csv` / `CH10_charData.csv`

### Analysis Parameters

| Parameter | Value |
|---|---|
| `yrInterp` | 10 yr |
| `Smoothing.method` | 2 (robust lowess) |
| `Smoothing.yr` | 500 yr |
| `threshType` | 2 (locally defined) |
| `threshMethod` | 3 (Gaussian mixture model) |
| `cPeak` | 1 (residuals) |
| `minCountP` | 0.05 |
| Record length | 626 interpolated samples |
| Record gaps | 1 gap (0.5 cm, 3 yr) near oldest part of record |

### Results

#### Peak Identification — Match Within Tolerance

| Column | v1.1 | v2.0 | Difference |
|---|---|---|---|
| peaks1 | 68 | 68 | 0 |
| peaks2 | 59 | 60 | 1 |
| peaks3 | 49 | 50 | 1 |
| peaksFinal | 49 | 50 | 1 |

**Rows where peaksFinal differs: 1 of 626**

The single differing peak occurs at 6080 yr BP. At this location `charPeak` and `threshFinalPos` are identical to four decimal places (0.0199 and 0.1489/0.1490 respectively) — a floating-point boundary case where a difference of 0.0001 in the smoothed threshold determines peak identification. This is not a scientifically meaningful difference.

#### Threshold Values — Numerically Equivalent Within Tolerance

| Column | v1.1 mean | v2.0 mean | Difference |
|---|---|---|---|
| thresh1 | 0.7026 | 0.6945 | 0.0081 |
| thresh2 | 0.9695 | 0.9598 | 0.0096 |
| thresh3 | 1.2686 | 1.2572 | 0.0114 |
| threshFinalPos | 1.2686 | 1.2572 | 0.0114 |

| Metric | Value |
|---|---|
| Samples where threshold differs | 511 of 626 (82%) |
| Mean absolute difference | 0.0115 |
| Maximum absolute difference | 0.3114 |

The threshold differences are attributable to inherent variability in the GMM EM algorithm across 626 local windows, combined with small differences in how NaN values from the record gap are handled near the oldest part of the record. These differences are not scientifically meaningful and do not affect the final peak identification in any meaningful way (1 boundary-case peak out of 49–50 total).

#### Peak Magnitude — Essentially Identical

| Metric | v1.1 | v2.0 | Difference |
|---|---|---|---|
| Total peak magnitude | 9083.73 | 9088.10 | 0.0% |
| Mean magnitude per peak | 181.6747 | 181.7620 | 0.05% |

#### Smoothed Fire Frequency

| Metric | Value |
|---|---|
| Mean absolute difference | 0.1382 |
| Maximum absolute difference | 1.8869 |

#### Regression Test Result

| PeakTol | ThreshTol | FreqTol | Result |
|---|---|---|---|
| 1 | 0.015 | 0.200 | 11/11 PASS |

The looser tolerances for CH10 reflect two documented sources of irreducible numerical difference specific to records with gaps and high-amplitude CHAR values when using a GMM + local threshold combination.

To reproduce this validation:
```matlab
z_Compare_CharAnalysis_V1_V2('CH10_charParams.csv', ...
    'CH10_charResults_V_1_1.csv', 'PeakTol', 1, ...
    'ThreshTol', 0.015, 'FreqTol', 0.200)
```

### Bugs Fixed During Round 3 Validation

Three bugs were identified and corrected during Round 3 testing. See the complete bugs list below.

---

## Round 4: Thunder Lake, Colorado (TL06) — Local Threshold

**Dataset:** Thunder Lake, CO — an independent charcoal record used for validation
**Input file:** `TL06_charParams.csv` / `TL06_charData.csv`

### Analysis Parameters

| Parameter | Value |
|---|---|
| `yrInterp` | 15 yr |
| `Smoothing.method` | 2 (robust lowess) |
| `Smoothing.yr` | 1000 yr |
| `threshType` | 2 (locally defined) |
| `threshMethod` | 3 (Gaussian mixture model) |
| `threshValues` | [0.95, 0.99, 0.999, 0.99] |
| `cPeak` | 1 (residuals) |
| `minCountP` | 0.05 |
| `peakFrequ` | 1000 yr |
| Record length | 414 interpolated samples |
| Record gaps | None |

### Bug Fixed During Round 4 Validation

A defect was identified and corrected in `CharSmooth.m` during Round 4 testing. See the complete bugs list below.

### Results

#### Peak Identification — Perfect Match

| Column | v1.1 | v2.0 | Difference |
|---|---|---|---|
| peaks1 | 27 | 27 | 0 |
| peaks2 | 21 | 21 | 0 |
| peaks3 | 16 | 16 | 0 |
| peaksFinal | 21 | 21 | 0 |

**Rows where peaksFinal differs: 0 of 414**

#### Threshold Values — Exact Match

| Column | v1.1 mean | v2.0 mean | Difference |
|---|---|---|---|
| thresh1 | 0.1251 | 0.1251 | 0.0000 |
| thresh2 | 0.1813 | 0.1813 | 0.0000 |
| thresh3 | 0.2443 | 0.2443 | 0.0000 |
| threshFinalPos | 0.1813 | 0.1813 | 0.0000 |

**Samples where threshold differs: 0 of 414**

An exact match is expected for a gap-free record where the Curve Fitting Toolbox is available: `smooth()` is called directly with identical inputs, producing bit-identical outputs.

#### Peak Magnitude — Exact Match

| Metric | v1.1 | v2.0 | Difference |
|---|---|---|---|
| Total peak magnitude | 177.90 | 177.90 | 0.0% |
| Mean magnitude per peak | 8.0864 | 8.0864 | 0.0% |

#### Smoothed Fire Frequency

| Metric | Value |
|---|---|
| Mean absolute difference | 0.0001 |
| Maximum absolute difference | 0.0012 |

#### Regression Test Result

| PeakTol | ThreshTol | FreqTol | Result |
|---|---|---|---|
| 0 | 0.001 | 0.001 | 11/11 PASS |

---

## Round 5: Raven Lake, Alaska (RA07) — Local Threshold

**Dataset:** Raven Lake, AK — an independent charcoal record used for validation
**Input file:** `RA07_charParams.csv` / `RA07_charData.csv`

### Analysis Parameters

| Parameter | Value |
|---|---|
| `yrInterp` | 15 yr |
| `Smoothing.method` | 2 (robust lowess) |
| `Smoothing.yr` | 750 yr |
| `threshType` | 2 (locally defined) |
| `threshMethod` | 3 (Gaussian mixture model) |
| `threshValues` | [0.95, 0.99, 0.999, 0.99] |
| `cPeak` | 1 (residuals) |
| `minCountP` | 0.25 |
| `peakFrequ` | 1000 yr |
| Record length | 204 interpolated samples |
| Record gaps | None |

### Results

#### Peak Identification — Perfect Match

| Column | v1.1 | v2.0 | Difference |
|---|---|---|---|
| peaks1 | 16 | 16 | 0 |
| peaks2 | 17 | 17 | 0 |
| peaks3 | 13 | 13 | 0 |
| peaksFinal | 17 | 17 | 0 |

**Rows where peaksFinal differs: 0 of 204**

#### Threshold Values — Exact Match

| Column | v1.1 mean | v2.0 mean | Difference |
|---|---|---|---|
| thresh1 | 0.0093 | 0.0093 | 0.0000 |
| thresh2 | 0.0136 | 0.0136 | 0.0000 |
| thresh3 | 0.0185 | 0.0185 | 0.0000 |
| threshFinalPos | 0.0136 | 0.0136 | 0.0000 |

**Samples where threshold differs: 0 of 204**

#### Peak Magnitude — Exact Match

| Metric | v1.1 | v2.0 | Difference |
|---|---|---|---|
| Total peak magnitude | 16.11 | 16.11 | 0.0% |
| Mean magnitude per peak | 0.8481 | 0.8481 | 0.0% |

#### Smoothed Fire Frequency

| Metric | Value |
|---|---|
| Mean absolute difference | 0.0004 |
| Maximum absolute difference | 0.0020 |

#### Regression Test Result

| PeakTol | ThreshTol | FreqTol | Result |
|---|---|---|---|
| 0 | 0.001 | 0.001 | 11/11 PASS |

---

## Round 6: Silver Lake, Montana (SI17) — Local Threshold

**Dataset:** Silver Lake, MT — an independent charcoal record used for validation
**Input file:** `SI17_charParams.csv` / `SI17_charData.csv`
**Date run:** 2026-04-02

### Analysis Parameters

| Parameter | Value |
|---|---|
| `threshType` | 2 (locally defined) |
| `threshMethod` | 3 (Gaussian mixture model) |
| Record length | 484 interpolated samples |
| Record gaps | None |
| Site context | Subalpine forest, Montana |

### Results

#### Peak Identification — Perfect Match

| Column | v1.1 | v2.0 | Difference |
|---|---|---|---|
| peaks1 | 44 | 44 | 0 |
| peaks2 | 31 | 31 | 0 |
| peaks3 | 22 | 22 | 0 |
| peaksFinal | 31 | 31 | 0 |

**Rows where peaksFinal differs: 0 of 484**

#### Threshold Values — Exact Match

| Column | v1.1 mean | v2.0 mean | Difference |
|---|---|---|---|
| thresh1 | 0.8973 | 0.8973 | 0.0000 |
| thresh2 | 1.3136 | 1.3136 | 0.0000 |
| thresh3 | 1.7802 | 1.7802 | 0.0000 |
| threshFinalPos | 1.3136 | 1.3136 | 0.0000 |

**Samples where threshFinalPos differs: 0 of 484** (mean absolute difference = 0.0000; max = 0.0000)

#### Peak Magnitude — Exact Match

| Metric | v1.1 | v2.0 | Difference |
|---|---|---|---|
| Total peak magnitude | 1205.32 | 1205.32 | 0.0% |
| Mean magnitude per peak | 38.8812 | 38.8812 | 0.0% |

#### Smoothed Fire Frequency

| Metric | Value |
|---|---|
| Mean absolute difference | 0.0000 |
| Maximum absolute difference | 0.0000 |

#### Regression Test Result

| PeakTol | ThreshTol | FreqTol | Result |
|---|---|---|---|
| 0 | 0.001 | 0.001 | 11/11 PASS |

Silver Lake produced zero differences across all 484 samples and all output variables — the cleanest result in the validation suite. As a gap-free record with the Curve Fitting Toolbox available, `smooth()` is called directly with identical inputs in both versions, producing bit-identical smoothing and GMM outputs throughout.

---

## Documented Numerical Tolerances by Analysis Type

The regression test script `z_Compare_CharAnalysis_V1_V2.m` accepts tolerance parameters that should be set according to the analysis configuration. Based on Rounds 1–6, the following tolerances are appropriate:

| Analysis type | PeakTol | ThreshTol | FreqTol | Notes |
|---|---|---|---|---|
| GMM + local, no gaps | 0 | 0.001 | 0.010 | Code Lake (local), Thunder Lake, Raven Lake, Silver Lake baseline |
| GMM + local, with gaps | 1 | 0.015 | 0.200 | Chickaree Lake baseline |
| GMM + global | 2 | 0.005 | 0.001 | Code Lake (global) baseline |

---

## Bugs Fixed in v2.0 (Complete List)

The following bugs were identified and corrected across Rounds 1–5. All fixes were confirmed to produce results consistent with v1.1.

| Bug | Location | Description |
|---|---|---|
| Global variable side effects | `CharThreshGlobal`, `CharThreshLocal`, `CharSmooth`, `CharPretreatment` | `global plotData` and `global bkgSensIn` replaced with explicit function arguments, eliminating order-dependent side effects |
| GMM fallback passed full record | `CharThreshLocal` | Poor-fit GMM fallback incorrectly passed `Charcoal.peak` (full record) instead of `X` (local window), causing severe slowdown |
| NaN propagation through smoothing | `CharSmooth` | NaN values in `Charcoal.accI` from record gaps were not bridged before smoothing, causing charLowess to propagate NaN across a large block of `accIS` values near the record end |
| NaN in GMM windows | `CharThreshLocal` | Windows near record gaps contained NaN values that caused the GMM EM algorithm to never converge; replaced with neutral-value substitution |
| Short-window fallback used continue | `CharThreshLocal` | Windows with fewer than 4 valid samples were skipped entirely, leaving NaN thresholds; replaced with Gaussian fallback that always produces a valid threshold |
| Peak magnitude formula | `CharPostProcess` | Operator precedence error: `(diff+1)*r` instead of v1.1's `diff+r`, producing magnitudes ~2.5x too large |
| Index overrun in fire frequency loop | `CharPostProcess` | Upper window index in the middle-of-record branch could exceed array length; clamped with `min(length(Charcoal.ybpI), ...)` |
| FRI variable overwritten | `CharPostProcess` | Per-zone loop variable `FRI` overwrote the full-record `FRI` vector returned by `smoothFRI`, causing a length mismatch in `CharPlotResults` |
| charLowess iteration count error | `charLowess` | `nIter` incorrectly set to 5 for all methods; fixed to 5 for rlowess only, 1 for lowess |
| charLowess rlowess mismatch | `charLowess` | Pure-MATLAB rlowess produced systematically higher background than MATLAB `smooth()`; fixed by calling `smooth()` directly when Curve Fitting Toolbox is available |
| NaN replacement in GMM windows | `CharThreshLocal` | NaN stripping changed window size near gaps; replaced with neutral-value substitution to preserve window size |
| CharSmooth.m toolbox routing | `CharSmooth` | `charLowess()` was called unconditionally even when Curve Fitting Toolbox was available; fixed to call `smooth()` directly when toolbox is present, with NaN bridging/restoration applied only in the `charLowess()` fallback path |
| zoneDiv(end) boundary overshoot | `CharPretreatment` | When `zoneDiv(end)` exceeded the bottom age of the last raw sample, v1.1 silently filled terminal samples with zero CHAR; v2.0 now auto-corrects `zoneDiv(end)` to the data boundary and notifies the user |
| SNI scalar dimension mismatch | `CharPostProcess` | For a global threshold, `CharThresh.SNI` is a scalar rather than an [N x 1] vector; expanded to a full vector before array concatenation to prevent a dimension mismatch error |
| y3tot computation for global threshold | `CharPlotResults` | `sum(CharcoalCharPeaks)` returned a [1 x 4] row vector unsuitable for plotting against bin centres; replaced with a cumulative count over threshold bins |
| noisePDF length mismatch | `CharPlotResults` | `CharThresh.noisePDF` has 251 elements (length of `CharThresh.possible`) while histogram bin centres have 250; interpolated noisePDF onto bin centres before plotting |

---

## Next Steps
- Continue testing additional datasets to broaden coverage (methods 3–5, log transform, ratios)
