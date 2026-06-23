# Chronological Uncertainty in Fire History Reconstructions

*Draft vignette — CharAnalysis 2.0 | In progress*

---

## 1. Motivation

Lake-sediment charcoal records reconstruct fire history by identifying statistically significant peaks in charcoal accumulation rates (CHAR). Each peak is assigned an age from an age-depth model, and the resulting sequence of peak ages constitutes a fire history: the timing, frequency, and variability of past fire at a site.

Age-depth models are not exact. Bayesian models such as rbacon (Blaauw and Christen 2011) generate a posterior distribution of possible age-depth relationships rather than a single curve. The standard practice of assigning a single "best estimate" age to each sample and treating those ages as fixed discards this distributional information — and the consequences extend beyond peak timing. Chronological uncertainty propagates into the analysis in two coupled ways: changing a sample's age changes the accumulation-rate denominator used to compute CHAR, and it changes the interpolated CHAR series on which background estimation and threshold detection operate. The result is that varying the chronology changes not just *when* peaks are detected, but *whether* they are detected at all. A simple post-hoc correction — shifting detected peak ages by the age-model uncertainty — is therefore insufficient; it treats peak detection as independent of the chronology when the two are coupled throughout the pipeline.

Propagating chronological uncertainty through the full CharAnalysis pipeline addresses this directly. Rather than asking "when did this fire occur?", the ensemble approach asks: "across the plausible range of age-depth histories for this site, how consistently does CharAnalysis detect a fire event at this time, and how precisely can we estimate when it occurred?" The outputs — detection frequency and a 95% confidence interval on peak age — are more honest representations of what the record can support.

This tool is most valuable when the scientific question centers on the timing of individual fire events. Assessments of overall fire-regime characteristics — mean fire-return interval, century-scale fire frequency — are less sensitive to chronological uncertainty in any single peak, because they are ensemble properties derived from many events, and because inherent variability in fire-return intervals is itself a meaningful signal that should not be conflated with chronological uncertainty.

---

## 2. Workflow

The ensemble workflow has four steps. Steps 2–4 are implemented in CharAnalysis; Step 1 requires an external age-depth modeling tool.

### Step 1: Build an age-depth model and extract a chronology ensemble

Run an age-depth model using rbacon (or another tool) and extract *n* = 1000 plausible chronologies from the MCMC posterior. Each chronology is a complete age-depth relationship: a vector of calibrated ages at every sampled depth in the core. Sampling directly from the MCMC posterior preserves the spatial correlation structure of the model — if one depth is assigned an older age in iteration *k*, neighboring depths are also older in that iteration.

The CharAnalysis function `char_extract_bacon_chronologies()` automates this extraction from a completed rbacon run:

```r
library(rbacon)

# Run Bacon in the same R session
Bacon("MySite",
      coredir      = "~/Bacon_runs/",
      acc.mean     = 10,
      acc.shape    = 1,
      mem.mean     = 0.5,
      mem.strength = 10,
      thick        = 5,
      ssize        = 8000,
      ask          = FALSE)

# Extract 1000 chronologies at charData depths
chardata <- read.csv("MySite_charData.csv", check.names = FALSE)

chron <- char_extract_bacon_chronologies(
  depths   = chardata$cmTop,
  n_iter   = 1000,
  out_file = "MySite_MCAgeDepth_1000_chronologies.csv"
)
```

The output CSV has columns `Sample_cm`, `CalAge_1`, ..., `CalAge_1000`. Any age-depth tool can supply this file — rbacon, OxCal, or a Monte Carlo interpolation approach (e.g., `MC_AgeDepth.m` for MATLAB users) — as long as the column structure matches.

### Step 2: Run the CharAnalysis ensemble

`char_run_ensemble()` (script: `char_run_ensemble.R`) reads the chronology ensemble and the site's parameter file, then runs the full CharAnalysis pipeline once per chronology. On each iteration, the charcoal data are re-assigned ages from chronology *k*, interpolated to a common time step, and processed through background estimation, threshold detection, and peak identification. The loop accumulates a peaks matrix (*n* time steps × 1000 iterations), plus matrices for C_peak and C_background.

```r
# Edit params_file and chron_file at the top of char_run_ensemble.R,
# then source:
source("CharAnalysis_2_0_R/tests/char_run_ensemble.R")
```

Parallel processing (via `parallel::parLapply`) reduces run time substantially. On a 20-core laptop, 1000 iterations complete in approximately 2 minutes.

### Step 3: Analyze the ensemble

`run_ensemble_analysis.R` post-processes the peaks matrix to characterize each fire event detected in a reference run (the median chronology). For each reference peak *k*, the script searches each of the 1000 iterations for the nearest detected peak within an adaptive window scaled to the local chronological uncertainty. This nearest-neighbor matching produces, for each reference peak:

- **P(peak)**: the fraction of iterations that detect a peak near the reference peak age
- **95% CI on peak age**: the 2.5th–97.5th percentile of detected ages across iterations
- **n_detected**: the raw count of iterations contributing to P(peak)

The script also identifies "orphan" peaks — detections that cluster across iterations at ages not associated with any reference peak. Orphans are classified as *near-reference* (secondary detection adjacent to a known peak) or *independent* (no reference peak within the search window), providing a conservative estimate of peaks that may have been missed in the single reference run.

### Step 4: Visualize

`plot_ensemble_figure.R` produces a four-panel figure:

- **(a)** CHAR time series with background and peak detections from the reference run
- **(b)** Detection frequency (%) for each peak type (reference, near-reference orphan, independent orphan), with 95% CI bars on timing
- **(c)** SNI (signal-to-noise index) trace
- **(d)** Age-depth model 95% CI ribbon from the chronology ensemble

---

## 3. Site example: CH10 (high-precision chronology)

CH10 is the CharAnalysis validation dataset — a high-resolution lake-sediment record with a well-constrained age-depth model built from closely spaced ²¹⁰Pb and radiocarbon dates (Kelly et al. 2011). The record spans approximately 6,200 cal yr BP and contains 59 peaks in the reference run (median chronology), with a mean fire-return interval of 104 yr.

Because the chronological uncertainty is low, the adaptive matching windows are narrow and uniform across the record: ±52 yr for every reference peak, driven by the mFRI floor rather than the 95% CI width. Peak timing uncertainty is correspondingly small: median 1 SD = 21 yr (range 8–30 yr across peaks). Detection frequencies are high — median 94.5%, with 47 of 59 peaks detected in ≥ 90% of iterations — reflecting that most fire events are robustly identified regardless of which plausible chronology is used. The ensemble mean FRI (median 98 yr, range 85–123 yr) closely tracks the reference value of 104 yr, as expected when peak counts are similar across iterations (ensemble median 62 peaks vs. reference 59).

![CH10 ensemble figure](../tests/CH10_ensemble_figure.png)
*Figure 1. Chronological uncertainty ensemble results for CH10. (a) CHAR time series with background and reference peak detections. (b) Detection frequency and timing uncertainty for each peak type; horizontal bars show ±95% CI (reference and independent peaks) or ±matching window (near-reference peaks). (c) Signal-to-noise index. (d) Chronological uncertainty ribbon (95% CI across the 1000-chronology ensemble).*

The narrow matching windows allow the secondary detection algorithm to distinguish genuine near-reference ensemble-only peaks from scatter. Six near-reference ensemble-only peaks were identified (at 4820, 2210, 1270, 1160, 780, and 480 cal yr BP), each with a 95% CI spanning less than 50% of the matching window (ci_frac < 0.5) — indicating that secondary detections cluster at a consistent age rather than being scattered across the window. Four independent ensemble-only peaks were also detected (5750, 5370, 660, and 460 cal yr BP), each corresponding to a fire event visible in the CHAR record at a lower detection threshold. CH10 illustrates what this tool can identify when chronological precision is high: a well-resolved reference fire history with a small number of additional candidate events surfaced by the ensemble.

---

## 4. Site example: SI17 (lower-precision chronology)

SI17 (Silver Lake) is a lake-sediment record spanning approximately 4,700 cal yr BP, with an age-depth model constructed from ²¹⁰Pb, radiocarbon, and tephra dates using rbacon (v4.1.1). Chronological uncertainty is substantially higher than at CH10, particularly in the middle and lower portions of the record, where the 95% CI on sample ages reaches several hundred years.

The reference run identifies 25 peaks (mean FRI 193 yr). Matching windows are wide and variable: ±105 to ±280 yr (median ±162 yr), scaled to the local 95% CI width at each reference peak age. Despite this uncertainty, detection frequencies remain high — median 99.0%, with all 25 peaks detected in ≥ 75% of iterations — indicating that the reference fire events are robustly identified across chronologies. Peak timing uncertainty is larger than at CH10: median 1 SD = 53 yr (range 18–123 yr).

The most notable difference from CH10 is in ensemble-only peak detection. No near-reference or independent ensemble-only peaks were identified. The wide matching windows absorb potential secondary detections near reference peaks, and any secondary detections that do occur are scattered across the window in different iterations rather than clustering at a consistent age (all candidates had ci_frac ≥ 0.76). This is the correct null result for a record with this level of chronological uncertainty: the ensemble cannot distinguish a real secondary event from the expected scatter of extra detections across a wide window. The ensemble mean FRI (median 149 yr, range 106–223 yr) is notably lower than the reference value (193 yr), reflecting that ensemble iterations detect a median of 32 peaks versus 25 in the reference run — a direct consequence of chronological uncertainty affecting threshold detection, not just peak timing.

SI17 illustrates the limits of secondary event detection under high chronological uncertainty, and demonstrates that the coherence filter (ci_frac < 0.5) correctly suppresses artifactual candidates rather than overclaiming ensemble-only events.

![SI17 ensemble figure](../tests/SI17_ensemble_figure.png)
*Figure 2. Chronological uncertainty ensemble results for SI17 (Silver Lake). See Figure 1 caption.*

---

### References

Blaauw, M., and J. A. Christen. 2011. Flexible paleoclimate age-depth models using an autoregressive gamma process. *Bayesian Analysis* 6:457–474.

Higuera, P. E., L. B. Brubaker, P. M. Anderson, F. S. Hu, and T. A. Brown. 2009. Vegetation mediated the impacts of postglacial climate change on fire regimes in the south-central Brooks Range, Alaska. *Ecological Monographs* 79:201–219.

Kelly, R. F., P. E. Higuera, C. M. Hart, and F. S. Hu. 2011. A signal-to-noise index to quantify the potential for peak detection in sediment-charcoal records. *Quaternary Research* 75:11–17.
