# *CharAnalysis* Development Roadmap
*Last updated: May 2026*
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
## 2. Record-quality diagnostics and user guidance

Improve guidance on interpreting peak detection results in light of record
quality, focused on the Signal-to-Noise Index (SNI; Kelly et al. 2011):

- Update the R vignette to highlight SNI as a first-step diagnostic to consult
  before moving forward with peak detection and interpretation. The vignette
  should make the rationale explicit: peak detection is only justified where
  the high-frequency component of the record is well separated from background
  noise.
- Have the R implementation return a console warning when portions of the
  record have SNI < 3 (the threshold proposed by Kelly et al. 2011). The
  warning is advisory, not blocking: the analyst retains the final decision.
  Suggested wording: `"x% of the record has SNI < 3; carefully consider
  whether peak analysis is appropriate in these sections."` The warning
  should be issued from within the function that computes SNI (the R
  equivalent of `CharThreshLocal.m`), not from a new diagnostic script.

---
## 3. Chronological Uncertainty
Incorporate methods for propagating chronological uncertainty into the
characterization of fire events. Development will take into account existing
approaches formalized in the following R packages, both of which were developed
based on *CharAnalysis* Version 1.1, and will include communication with those
developers:
- `tapas`: https://github.com/wfinsinger/tapas
- `CharcoalFireReconstructionR`: https://github.com/rglueckler/CharcoalFireReconstructionR
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
