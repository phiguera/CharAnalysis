# *CharAnalysis* Development Roadmap
*Last updated: April 2026*

---

This document describes planned future development of *CharAnalysis*. Items are
listed in approximate priority order. It is not yet determined whether additional
development beyond the R translation will occur in MATLAB, R, or both.

---

## 1. R Translation ✓ *Complete — v2.0 (April 2026)*

*CharAnalysis* v2.0 is a direct R translation of the MATLAB v2.0 implementation,
released to CRAN in April 2026. The package is in the
[experimental](https://lifecycle.r-lib.org/articles/stages.html#experimental)
lifecycle stage; the API may evolve as user feedback is incorporated.

The R package (`CharAnalysis_2_0_R/`) reproduces the full five-step analytical
pipeline and all eight output figures. Key design choices:

- **Quantitative equivalence**: validated against MATLAB v2.0 on four benchmark
  datasets (Code Lake, Chickaree Lake, Silver Lake, Raven Lake). Known numerical
  differences are documented in `NEWS.md` and `inst/z_Validation_report_R_vs_MATLAB_V_2.0.md`.
- **GMM implementation**: the Gaussian mixture model threshold uses a direct port
  of the MATLAB `GaussianMixture.m` EM algorithm rather than an existing R package,
  preserving numerical comparability with the reference implementation.
- **Figures**: eight publication-quality ggplot2 figures mirroring the MATLAB output.
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

## 2. Chronological Uncertainty

Incorporate methods for propagating chronological uncertainty into the
characterization of fire events. Development will take into account existing
approaches formalized in the following R packages, both of which were developed
based on *CharAnalysis* Version 1.1, and will include communication with those
developers:

- `tapas`: https://github.com/wfinsinger/tapas
- `CharcoalFireReconstructionR`: https://github.com/rglueckler/CharcoalFireReconstructionR

---

## 3. Regional Synthesis

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
