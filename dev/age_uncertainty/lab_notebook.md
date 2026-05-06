# CharAnalysis age-uncertainty integration : lab notebook

Working document for the rplum integration project. Newest entries at the top.

---

## 2026-05-06 : project kickoff

### Question
How should age-model chronological uncertainty be propagated into CharAnalysis peak detection, given that current peak detection treats sediment ages as point estimates?

### Tool decision
Use `rplum` (Blaauw, Aquino-López, Christen, et al.) rather than `rbacon` alone. rplum is a companion to rbacon, not a successor. It extends the Bacon framework to handle 210Pb chronologies via the Plum method (Aquino-López et al. 2018, *JABES*) and can combine 210Pb, 14C, and other dates in a single Bayesian age-depth model. Most North American lake-sediment charcoal records include 210Pb in the upper sediment, which is exactly the interval where chronological uncertainty is largest and where modern fire-regime comparisons most often live.

### Repo organization
- Branch: `feature/age-uncertainty` off `dev`.
- Working directory: `dev/age_uncertainty/` at repo root, sibling to the package directories. Convention follows tidyverse practice: scratch scripts and notes are version-controlled but isolated from the released package source tree.
- Version-bump decision deferred. If integration is additive (new function or new optional arg), 2.x minor bump. If integration changes existing function signatures or default return structures, 3.0 major bump. Will commit at merge time, not now.

### Open design questions

1. **Posterior ensemble size.** How many age-model posterior realizations to carry forward into peak detection? Common range in published charcoal-uncertainty work is 100 to 1000. Need to evaluate the trade-off between Monte Carlo noise floor and compute cost per record.

2. **Pipeline scope per draw.** Per posterior age-model realization, recompute (a) only age assignments to existing CHAR values, or (b) the full CHAR-to-peaks pipeline including smoothing and local-threshold determination? Option (b) captures feedback between chronology and threshold selection but multiplies cost by the smoothing/thresholding cost.

3. **API design.**
   - *Additive*: new function `char_peak_id_uncertain()` alongside the existing point-estimate function, or new optional `age_posterior` argument to `char_peak_id()`. → 2.x minor bump.
   - *Breaking*: existing peak-detection functions return uncertainty-aware structures by default. → 3.0 major bump.

4. **Output structure.** Per-peak posterior probability of being a peak? Posterior age range per peak? Both? How does this then propagate through to FRI distributions, which CharAnalysis currently summarizes as point intervals?

5. **Validation strategy.** Need at least one test record where age uncertainty is well-characterized for sanity checks. Candidates from the existing validation set: CO, CH10, RA07, SI17. Or simulated data with known peaks and known age structure (cleaner for testing but lower ecological realism).

### Next steps
- Literature review: existing examples of age-uncertainty propagation in charcoal records.
  - Higuera et al. 2009 (*Ecological Monographs*): original CharAnalysis approach, point-estimate ages.
  - Kelly et al. 2013 (*PNAS*): did the Yukon Flats analysis incorporate age uncertainty? Worth checking the methods.
  - Blarquez & Carcaillet 2010: sensitivity-style framing of charcoal records.
  - Dietze et al.: papers using rbacon ensembles to propagate chronology uncertainty into downstream proxies.
- Confirm rplum installs cleanly on the current R setup. Run the package's example dataset to sanity-check the workflow.
- Sketch API options on paper before writing code. Decision criterion: does the additive option leave any awkward gaps where users would have to call two functions to get an uncertainty-aware result?

### References
- Aquino-López, M. A., Blaauw, M., Christen, J. A., Sanderson, N. K. (2018). Bayesian analysis of 210Pb dating. *Journal of Agricultural, Biological, and Environmental Statistics* 23, 317-333.
- rplum on CRAN: https://cran.r-project.org/package=rplum
- rbacon on CRAN: https://cran.r-project.org/package=rbacon
