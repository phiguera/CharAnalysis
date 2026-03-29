# *CharAnalysis* Development Roadmap
*Last updated: March 2026*

---

This document describes planned future development of *CharAnalysis*. Items are
listed in approximate priority order. It is not yet determined whether additional
development beyond the R translation will occur in MATLAB, R, or both.

---

## 1. R Translation

Translate *CharAnalysis* Version 2.0 into R, the dominant language in the
paleoecological community. The R implementation will use modern packages including
`mclust` for Gaussian mixture modeling and `ggplot2` for figure output. A central
goal of this effort is quantitative equivalence with the MATLAB Version 2.0
results, so that users can move between platforms with confidence that results are
comparable. The R translation will be explicitly benchmarked against MATLAB
Version 2.0 outputs, and will possibly also be compared — qualitatively or
quantitatively — against existing R packages that were developed based on
*CharAnalysis* Version 1.1, including:

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

*Development of Version 2.0 and planning for future updates was carried out with
the assistance of Claude, an AI assistant by Anthropic. All work has been reviewed
and validated by the author.*
