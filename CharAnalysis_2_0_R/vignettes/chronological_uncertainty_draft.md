# Integrating Chronological Uncertainty into *CharAnalysis*

*Draft vignette — CharAnalysis 2.0 \| In progress*

------------------------------------------------------------------------

## 1. Motivation and Rationale

Age uncertainty from age-depth models affects not just sample age, but also charcoal accumulation rate (CHAR) and thus several elements of threshold determination and peak identification executed in *CharAnalysis*.

This tool re-runs *CharAnalysis* across 1,000 plausible chronologies, propagating chronological uncertainty across all steps. The output addresses three questions, each integrating the uncertainty in the underlying age-depth model: (i) how consistently is a peak identified across iterations; (ii) what is the range of ages assigned to an individual peak across iterations; and (iii) are there peaks commonly identified among the iterations that are not identified in the main *CharAnalysis* run?

This tool is most valuable when questions center on the identification and timing of individual peaks. Fire-regime statistics (e.g., mean fire-return interval, century-scale fire frequency) are less sensitive to per-peak chronological uncertainty. For those regime-scale statistics, the main *CharAnalysis* run based on the median age-depth model provides the most plausible outcome.

------------------------------------------------------------------------

## 2. Workflow

The ensemble workflow has four steps. Steps 2–4 are implemented in *CharAnalysis*; Step 1 requires an external age-depth modeling tool.

### Step 1: Build an age-depth model and extract a chronology ensemble

Run an age-depth model (e.g., rbacon, Blaauw & Christen 2011; MCAgeDepth, Higuera et al. 2009) and extract *n* = 1,000 chronologies from the MCMC runs. Each chronology is a complete age-depth model at all sampled depths, preserving the temporal correlation structure of the model and reflecting the uncertainty of the dates (e.g., ¹⁴C, ²¹⁰Pb) informing it.

The *CharAnalysis* function `char_extract_bacon_chronologies()` automates extraction from a completed rbacon run:

``` r
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

The output CSV has columns `Sample_cm`, `CalAge_1`, ..., `CalAge_1000`. Any age-depth tool can supply this file (e.g. rbacon, OxCal, MCAgeDepth) as long as the column structure matches.

### Step 2: Run the *CharAnalysis* ensemble

`char_run_ensemble.R` reads the chronology ensemble and site parameter file, then runs the full *CharAnalysis* pipeline once per chronology. On each iteration, sample ages are reassigned from chronology *k*, charcoal accumulation rates (CHAR) are recalculated using the updated age-depth relationship, and the resulting CHAR series is processed through background estimation, threshold determination, and peak identification. Output is a peaks matrix (*n* time steps × 1,000 iterations) plus matrices for C_peak and C_background. With parallel processing, 1,000 iterations complete in approximately 2 minutes on a modern laptop.

``` r
# Edit params_file and chron_file at the top of char_run_ensemble.R,
# then source:
source("CharAnalysis_2_0_R/tests/char_run_ensemble.R")
```

### Step 3: Analyze the ensemble

`run_ensemble_analysis.R` characterizes each fire event identified in the benchmark *CharAnalysis* run — the run based on the most likely (median) age-depth model. For each benchmark peak, the script searches all 1,000 ensemble iterations for the nearest detected peak within an adaptive window scaled to local chronological uncertainty, producing:

- **Detection frequency**: fraction of iterations detecting a peak near the benchmark peak age
- **95% CI on peak age**: 2.5th–97.5th percentile of detected ages across iterations
- **Ensemble-only peaks**: events detected in ≥10% of iterations but absent from the benchmark run, flagged as candidate fire events that the median chronology may have missed; the 10% threshold is a user-adjustable parameter, and whether flagged peaks warrant interpretation depends on their detection frequency and temporal coherence across the ensemble

**Adaptive matching window.** For each benchmark peak, the search window defines how far — in either direction — an ensemble detection can be from the benchmark peak age and still be attributed to the same peak. Any detection within ±window~k~ of the benchmark peak age is attributed to it. The window is:

> window~k~ = max(mFRI~z~ / 2, CI~95,k~ / 2)

where mFRI~z~ is the mean fire-return interval for the zone containing peak *k* (from the benchmark run) and CI~95,k~ is the 95% width of the age-depth model uncertainty at the age of peak *k*, interpolated from the chronology ensemble. The mFRI~z~/2 floor reflects the expectation that once a peak shifts more than half the mean FRI, it is likely a different event. The CI~95,k~/2 term expands the window where age-depth uncertainty is high, ensuring that ensemble detections of the same peak — shifted in time by a poorly constrained chronology — are still correctly attributed to the benchmark peak. In practice, CI~95,k~/2 is the larger of the two terms in records with meaningful age uncertainty; mFRI~z~/2 serves as a minimum bound in well-constrained records where CI~95,k~ is small. The window for each benchmark peak is reported in the output table.

### Step 4: Visualize

`plot_ensemble_figure.R` produces a four-panel figure:

- **(a)** CHAR time series with background and peak detections from the benchmark run
- **(b)** Detection frequency (%) for each peak type (reference, near-reference orphan, independent orphan), with 95% CI bars on timing
- **(c)** SNI (signal-to-noise index) trace
- **(d)** Age-depth model 95% CI ribbon from the chronology ensemble

`plot_chronUncertainty_methods.R` produces the methods illustration figure described below (Figure 1).

![Methods illustration](../tests/CH10_chronUncertainty_methods.png) **Figure 1.** Methods illustration for the chronological uncertainty ensemble, shown for the most recent 3,000 years of the CH10 record (full record spans ~6,200 cal yr BP). *Row 1 (benchmark run):* CHAR time series (black bars), background estimate (grey line), detection thresholds (red lines), and benchmark peak detections ("+"). *Rows 2–6 (five ensemble iterations):* the same CHAR content recalculated under five randomly drawn chronologies. Detected peaks are shown as filled circles with horizontal bars spanning the adaptive matching window (±window~k~): black circles are matched to a benchmark peak; grey circles are unmatched. Vertical grey lines mark benchmark peak ages for reference. *Bottom panel:* detection frequency (%) for each benchmark peak and ensemble-only peak clusters (grey open circles) identified in ≥40% of the five iterations; horizontal bars show the 95% CI on timing where estimable.

------------------------------------------------------------------------

## 3. Site examples

CH10 and SI17 illustrate how chronological precision shapes ensemble behavior. CH10 is unusually well-dated for a lake-sediment record, with closely spaced dates producing narrow, uniform uncertainty throughout. SI17 is less precisely dated, with wider and more variable uncertainty driven by date spacing and a calibration curve feature near 1,000 cal yr BP — though it remains a well-dated record by typical standards. Together they span a range of precision that captures the experience of most lake-sediment fire history studies, from narrow matching windows and high detection consistency to wide windows and substantial timing uncertainty.

![Age-depth comparison](../tests/agedepth_vignette.png) **Figure 2.** Age-depth models for CH10 (blue) and SI17 (black), each derived from the 1,000-chronology ensemble. Lines show the median chronology; shaded ribbons show the 95% CI. Ribbon width reflects chronological precision: narrow and uniform in CH10, where closely spaced dates constrain the model throughout; wider and more variable in SI17, particularly near 1,000 cal yr BP where calibration curve geometry broadens uncertainty. The shared axes allow direct comparison of sediment accumulation rates (line slopes) and chronological precision (ribbon widths) between sites.

| | CH10 | SI17 |
|:--|--:|--:|
| Record length (cal yr BP) | ~6,200 | ~4,700 |
| Benchmark peaks (*n*) | 59 | 25 |
| Mean FRI (yr) | 104 | 193 |
| Matching window (median ± yr) | ±52 | ±162 |
| Detection frequency (median %) | 95 | 99 |
| Peak timing uncertainty (median 1 SD, yr) | 21 | 53 |
| Ensemble-only peaks detected | Yes | No |

### 3.1 CH10 (unusually well-dated)

CH10 is a lake-sediment record from a Rocky Mountain subalpine watershed spanning ~6,200 cal yr BP, with an age-depth model built from closely spaced ²¹⁰Pb and radiocarbon dates (Dunnette et al. 2014, Leys et al. 2016). Chronological uncertainty is low and nearly uniform across the record.

The benchmark run identifies 59 peaks (mean FRI 104 yr). Matching windows are narrow and uniform (±52 yr), driven by the mFRI floor rather than chronological uncertainty. Detection frequencies are high (median 95%; 47 of 59 peaks detected in ≥90% of iterations) and peak timing uncertainty is small (median 1 SD = 21 yr). The ensemble surfaces a small number of additional candidate fire events absent from the benchmark run. CH10 represents the best case: low uncertainty, high detection consistency, and interpretable ensemble-only candidates.

![CH10 ensemble figure](../tests/CH10_ensemble_figure.png) **Figure 3.** Chronological uncertainty ensemble results for CH10. (a) CHAR time series with background and benchmark peak detections. (b) Detection frequency and timing uncertainty for each peak type; horizontal bars show ±95% CI on peak age. (c) Signal-to-noise index. (d) Chronological uncertainty ribbon (95% CI across the 1,000-chronology ensemble).

### 3.2 SI17 (less precisely dated)

SI17 (Silver Lake, Colorado) spans ~4,700 cal yr BP, with an age-depth model built from ²¹⁰Pb, radiocarbon, and tephra dates (Clark-Wolf et al. 2023). Chronological uncertainty is substantially higher than CH10, particularly near 1,000 cal yr BP where the IntCal calibration curve produces asymmetric age distributions.

The benchmark run identifies 25 peaks (mean FRI 193 yr). Matching windows are wide and variable (±105–280 yr, median ±162 yr), dominated by the CI~95~/2 term. Detection frequencies remain high (median 99%), confirming that fire events are robustly identified even when timing is less precise (median 1 SD = 53 yr). No ensemble-only peaks pass the coherence filter: wide matching windows absorb secondary detections as scatter rather than coherent signal. Ensemble iterations detect more peaks on average (median 32) than the benchmark run (25), a direct consequence of chronological uncertainty affecting CHAR magnitude and threshold determination, not just peak timing.

![SI17 ensemble figure](../tests/SI17_ensemble_figure.png) **Figure 4.** Chronological uncertainty ensemble results for SI17 (Silver Lake). See Figure 3 caption.

------------------------------------------------------------------------

### References

Blaauw, M., and J. A. Christen. 2011. Flexible paleoclimate age-depth models using an autoregressive gamma process. *Bayesian Analysis* 6:457–474.

Higuera, P. E., L. B. Brubaker, P. M. Anderson, F. S. Hu, and T. A. Brown. 2009. Vegetation mediated the impacts of postglacial climate change on fire regimes in the south-central Brooks Range, Alaska. *Ecological Monographs* 79:201–219.

Clark-Wolf, K. D., P. E. Higuera, K. K. McLauchlan, B. N. Shuman, and M. C. Parish. 2023. Fire-regime variability and ecosystem resilience over four millennia in a Rocky Mountain subalpine watershed. *Journal of Ecology* 111:2643–2661.

Dunnette, P. V., P. E. Higuera, K. K. McLauchlan, K. M. Derr, C. E. Briles, and M. H. Keefe. 2014. Biogeochemical impacts of wildfires over four millennia in a Rocky Mountain subalpine watershed. *New Phytologist* 203:900–912.

Leys, B., P. E. Higuera, K. K. McLauchlan, and P. V. Dunnette. 2016. Wildfires and geochemical change in a subalpine forest over the past six millennia. *Environmental Research Letters* 11:125003.
