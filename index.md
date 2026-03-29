---
layout: default
title: CharAnalysis
---

# *CharAnalysis*
## Diagnostic and analytical tools for peak detection in sediment-charcoal records

[![CC BY 4.0](https://licensebuttons.net/l/by/4.0/88x31.png)](https://creativecommons.org/licenses/by/4.0/)
© 2004–2026 Philip Higuera  
Professor, Department of Ecosystem and Conservation Sciences  
University of Montana, Missoula, MT, USA  
[philip.higuera@umontana.edu](mailto:philip.higuera@umontana.edu) |
[Faculty page](https://www.umt.edu/people/phiguera) |
[GitHub repository](https://github.com/phiguera/CharAnalysis)

---

## What is *CharAnalysis*?

*CharAnalysis* is a freely available program for reconstructing local fire
histories from lake-sediment charcoal records. It implements a widely applied
approach that decomposes a charcoal record into low- and high-frequency
components, and introduced the use of a locally defined threshold to separate
fire signal from noise (Higuera et al. 2008, 2009). The program is designed to
make the full range of analytical choices explicit, with diagnostic tools to
guide parameter selection and sensitivity analyses to evaluate the robustness of
results.

Since its original development in the mid-2000s, *CharAnalysis* has been used in
dozens of published studies to analyze sediment-charcoal records on six
continents. The entire codebase is distributed and well commented — users are
encouraged to look under the hood, understand what is going on, and modify the
program to suit their needs.

---

## What does it do?

*CharAnalysis* takes a raw sediment-charcoal record and guides the user through
five analytical steps: interpolating the record to equal time intervals,
smoothing to estimate background charcoal accumulation, isolating the
high-frequency peak component, applying a threshold to identify fire events, and
screening peaks using a minimum-count criterion. At each step the analyst makes
explicit parameter choices, informed by diagnostic output from the program.

The figures below illustrate typical program output for the Code Lake record
from the south-central Brooks Range, Alaska (Higuera et al. 2009).

![CharAnalysis output: sensitivity to alternative thresholds](images/fig04_sensitivity_sni.png)
*Figure 1. Sensitivity of peak identification to alternative threshold values
(top), mean fire return intervals by zone for each threshold (second panel),
signal-to-noise index through time (third panel), and boxplot of all SNI values
(bottom). The SNI quantifies the potential for reliable peak detection at each
point in the record.*

![CharAnalysis output: continuous fire history](images/fig07_continuous_fire_history.png)
*Figure 2. Continuous fire history showing peak magnitude (top), fire return
intervals and smoothed FRI curve (middle), and smoothed fire frequency (bottom).*

---

## Getting Started

There are three ways to access *CharAnalysis*, suited to different users and
needs. Full installation and usage instructions are in the
[User's Guide](https://github.com/phiguera/CharAnalysis/blob/master/CharAnalysis_Us
