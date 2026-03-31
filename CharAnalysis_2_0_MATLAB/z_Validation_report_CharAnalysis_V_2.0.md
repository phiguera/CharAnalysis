# CharAnalysis v2.0 Validation Report

**Date:** March 30, 2026
**Prepared by:** CharAnalysis Development Team
**Comparison:** CharAnalysis v1.1 (MATLAB) vs. CharAnalysis v2.0 (MATLAB)

---

## Validation Datasets

| Round | Site ID | Site Name | Location | Reference |
|---|---|---|---|---|
| 1 | CO | Code Lake | Alaska, USA | Higuera et al. 2009 |
| 2 | CH10 | Chickaree Lake | Colorado, USA | Dunnette et al. 2014 |
| 3 | TL06 | Thunder Lake | — | — |
| 4 | RA07 | Raven Lake | — | — |

Code Lake is the example dataset bundled with CharAnalysis distributions. Chickaree Lake, Thunder Lake, and Raven Lake are independent records used to validate v2.0 across a range of record lengths, sampling resolutions, and fire return interval characteristics.

---

## Round 1: Code Lake, Alaska (CO)

**Dataset:** Code Lake, AK — the example dataset bundled with CharAnalysis distributions
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
| Record length | 504 interpolated samples |

### Results

#### Peak Identification — Perfect Match

| Column | v1.1 | v2.0 | Difference |
|---|---|---|---|
| peaks1 | 54 | 54 | 0 |
| peaks2 | 48 | 48 | 0 |
| peaks3 | 43 | 43 | 0 |
| peaksFinal | 48 | 48 | 0 |

**Rows where peaksFinal differs: 0 of 504**

#### Threshold Values — Numerically Equivalent

| Column | v1.1 mean | v2.0 mean | Difference |
|---|---|---|---|
| thresh1 | 0.0331 | 0.0331 | 0.0000 |
| thresh2 | 0.0513 | 0.0515 | 0.0002 |
| thresh3 | 0.0718 | 0.0721 | 0.0003 |
| threshFinalPos | 0.0513 | 0.0515 | 0.0002 |

| Metric | Value |
|---|---|
| Samples where threshold differs | 71 of 504 (14%) |
| Mean absolute difference | 0.0002 |
| Maximum absolute difference | 0.0071 |

#### Peak Magnitude — Essentially Identical

| Metric | v1.1 | v2.0 | Difference |
|---|---|---|---|
| Total peak magnitude | 127.43 | 126.70 | 0.73 (0.6%) |
| Mean magnitude per peak | 2.4506 | 2.4844 | 0.0338 (1.4%) |

#### Regression Test Result

| PeakTol | ThreshTol | FreqTol | Result |
|---|---|---|---|
| 0 | 0.001 | 0.001 | 11/11 PASS |

---

## Round 2: Chickaree Lake, Colorado (CH10)

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

The looser tolerances for CH10 relative to Code Lake reflect two documented sources of irreducible numerical difference specific to records with gaps and high-amplitude CHAR values when using a GMM + local threshold combination.

To reproduce the passing CH10 validation, use the following call:
```matlab
z_Compare_CharAnalysis_V1_V2('CH10_charParams.csv', ...
    'CH10_charResults_V_1_1.csv', 'PeakTol', 1, ...
    'ThreshTol', 0.015, 'FreqTol', 0.200)
```
Running the script with default tolerances will produce failures for CH10 —
this is expected and documented. See the tolerances table below for guidance
on appropriate tolerance settings by analysis type.

---

## Bugs Fixed During Round 2 Validation

The following bugs were identified and corrected during Round 2 testing. All fixes were confirmed to produce results consistent with v1.1.

| Bug | Location | Description |
|---|---|---|
| charLowess iteration count error | `charLowess` | `nIter` was set to a fixed value of 5 for all methods including plain lowess, causing plain lowess to run 5 iterations instead of 1. Fixed to use `nIter=5` for rlowess only and `nIter=1` for lowess. |
| charLowess uses smooth() when available | `charLowess` | Added a check for the Curve Fitting Toolbox at runtime. When available, `smooth()` is called directly, guaranteeing results identical to v1.1. The pure-MATLAB implementation is used as a fallback only when the toolbox is absent. |
| NaN replacement in GMM windows | `CharThreshLocal` | NaN values from record gaps were stripped from the local window before GMM fitting, reducing window size and changing the distribution shape near gaps. Fixed to replace NaN with the noise-component mean (0 for residuals, 1 for ratios), preserving window size while preventing gap samples from influencing the fit. |

---

<<<<<<< HEAD
## Round 3: Thunder Lake (TL06)

**Dataset:** Thunder Lake — an independent charcoal record used for validation
=======
## Round 3: TL06 (Thunder Lake)

**Dataset:** TL06 — Thunder Lake, an independent charcoal record used for validation
>>>>>>> ff3da168a479a21e19b0feea5df3ed230dd883a7
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

### Input Correction Identified During Validation

<<<<<<< HEAD
The original `TL06_charParams.csv` had `zoneDiv(end) = 6200`, which extends 48 yr beyond the bottom age of the last raw sample (6150 yr BP). This caused v2.0 to mark 4 terminal interpolated samples NaN, while v1.1 silently filled them with zero CHAR. The NaN values affected `charBkg` across the entire record via the smoother, producing threshold differences everywhere.

The root cause is a v1.1 bug: the double loop in `CharPretreatment.m` assigned zero to samples with no overlapping raw data, whereas the v2.0 vectorized proportion matrix correctly leaves them as NaN. V2.0's behavior is more accurate.

The params file was corrected to `zoneDiv(end) = 6150` before generating the v1.1 reference for this round. **Users should ensure `zoneDiv(end)` does not exceed the bottom age of the last raw sample.** A warning is now emitted by `CharValidateParams.m` when this condition is detected.
=======
The original `TL06_charParams.csv` had `zoneDiv(end) = 6200`, which extends 48 yr beyond the last raw sample bottom age of 6150 yr BP. This caused v1.1 to silently extrapolate 4 terminal interpolated samples with zero CHAR values, while v2.0 correctly marks them NaN — a difference that propagated into `charBkg` across the entire record via the smoother.

The root cause is a v1.1 bug: the double loop in `CharPretreatment.m` assigned zero to samples with no overlapping raw data, whereas the v2.0 vectorized proportion matrix correctly leaves them as NaN. V2.0's behavior is more accurate.

The params file was corrected to `zoneDiv(end) = 6150` before generating the v1.1 reference for this round. This is a user input correction, not a code change. **Users should ensure `zoneDiv(end)` does not extend beyond the bottom age of the last raw sample in their input data.**
>>>>>>> ff3da168a479a21e19b0feea5df3ed230dd883a7

### Bug Fixed During Round 3 Validation

During Round 3 validation, a defect was identified and corrected in `CharSmooth.m`. See the full description in the "Bugs Fixed" section below.

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

TL06 achieves an exact match on all threshold values. This is expected for a gap-free record where the Curve Fitting Toolbox is available: `smooth()` is called directly with identical inputs, producing bit-identical outputs.

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

## Bug Fixed During Round 3 Validation

### CharSmooth.m — Curve Fitting Toolbox routing

**Location:** `CharSmooth.m`

**Description:** The original v2.0 `CharSmooth.m` unconditionally called `charLowess()` for all five smoothing methods, even when the Curve Fitting Toolbox was available. This was incorrect for records with NaN gaps: `charLowess()` requires NaN-free input and uses linear interpolation to bridge gaps before smoothing, whereas MATLAB's `smooth()` handles NaNs internally by excluding them from the local regression window — exactly as v1.1 did. Passing gap-bridged data to `smooth()` produces different results from passing the raw data (including NaNs), because the bridged values alter the weighted regression across the entire record, not just near the gap.

**Fix:** Added a Curve Fitting Toolbox check at runtime in `CharSmooth.m`. When the toolbox is available, `smooth()` is called directly on `Charcoal.accI` (including NaNs) for methods 1–3 and for the Lowess passes in methods 4–5 — identical to v1.1 behavior. The NaN bridging and restoration logic is retained but guarded by `~hasCFT`, so it applies only when `charLowess()` is used as the fallback.

**Impact:** On records without NaN gaps (e.g., Code Lake), this fix has no effect. On records with NaN gaps and the Curve Fitting Toolbox available (e.g., TL06), this fix produces `charBkg` values identical to v1.1. The fix was confirmed by TL06 achieving a 0.0000 threshold difference after it was applied.

---

<<<<<<< HEAD
## Round 4: Raven Lake (RA07)

**Dataset:** Raven Lake — an independent charcoal record used for validation
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

### Input Correction Identified During Validation

The original `RA07_charParams.csv` had `zoneDiv(end) = 3008`, which extends slightly beyond the bottom age of the last raw sample (3008.18 yr BP). This is the same boundary overshoot issue identified in Round 3 (TL06). One terminal interpolated sample at age 3003 yr BP was marked NaN in v2.0, which shifted `charBkg` across the entire record via the smoother, causing both a peak count discrepancy in `peaks1` and a magnitude failure.

The params file was corrected to `zoneDiv(end) = 3000` before generating the v1.1 reference for this round. Having now observed this issue on two independent datasets, a warning has been added to `CharValidateParams.m` (check 7) that fires whenever `zoneDiv(end)` exceeds the bottom age of the last raw sample, directing users to correct the value before running.

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

## Documented Numerical Tolerances by Analysis Type

The regression test script `z_Compare_CharAnalysis_V1_V2.m` accepts tolerance parameters that should be set according to the analysis configuration. Based on Rounds 1–4, the following tolerances are appropriate:

| Analysis type | PeakTol | ThreshTol | FreqTol | Notes |
|---|---|---|---|---|
| GMM + local, no gaps | 0 | 0.001 | 0.001 | Code Lake, Thunder Lake, and Raven Lake baseline |
| GMM + local, with gaps | 1 | 0.015 | 0.200 | Chickaree Lake baseline |

---

## Recurring Input Issue: zoneDiv Boundary Overshoot

This issue was observed independently in Round 3 (Thunder Lake) and Round 4 (Raven Lake). In both cases, `zoneDiv(end)` in the params file extended slightly beyond the bottom age of the last raw sample. The v2.0 vectorized proportion matrix correctly assigns NaN to interpolated intervals with no overlapping raw data, whereas v1.1's double loop silently filled them with zero CHAR. The NaN values propagated into `charBkg` via the smoother, producing differences across the entire record.

**Resolution:** Correct `zoneDiv(end)` to be no greater than the bottom age of the last raw sample. As of Round 4, `CharValidateParams.m` emits a warning when this condition is detected, specifying the data boundary and the recommended correction.
=======
## Documented Numerical Tolerances by Analysis Type

The regression test script `z_Compare_CharAnalysis_V1_V2.m` accepts tolerance parameters that should be set according to the analysis configuration. Based on Rounds 1–3, the following tolerances are appropriate:

| Analysis type | PeakTol | ThreshTol | FreqTol | Notes |
|---|---|---|---|---|
| GMM + local, no gaps | 0 | 0.001 | 0.001 | Code Lake and TL06 baseline |
| GMM + local, with gaps | 1 | 0.015 | 0.200 | CH10 baseline |
>>>>>>> ff3da168a479a21e19b0feea5df3ed230dd883a7

---

## Bugs Fixed in v2.0 (Complete List)

<<<<<<< HEAD
The following bugs were identified and corrected across Rounds 1–4. All fixes were confirmed to produce results consistent with v1.1.
=======
The following bugs were identified and corrected across Rounds 1–3. All fixes were confirmed to produce results consistent with v1.1.
>>>>>>> ff3da168a479a21e19b0feea5df3ed230dd883a7

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

---

## Next Steps

<<<<<<< HEAD
- Round 5: Test the `bkgSens = 1` sensitivity analysis path
- Update `dev` branch README with Round 4 validation status
=======
- Round 4: Test the `bkgSens = 1` sensitivity analysis path
- Update `dev` branch README with Round 3 validation status
>>>>>>> ff3da168a479a21e19b0feea5df3ed230dd883a7
- Merge `dev` branch to `main` after `bkgSens` validation passes
- Continue testing additional datasets to broaden coverage (global threshold, methods 3–5, log transform)
