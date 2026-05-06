# CharAnalysis age-uncertainty integration : lab notebook

Working document for the rplum integration project. Newest entries at the top.

## Project goal

For each sample (depth) in a charcoal record, compute P(peak) = the fraction of N posterior age-model realizations in which CharAnalysis identifies that sample as a peak, where N is on the order of 100 to 1000 draws from the age-depth posterior (rplum or comparable). The output is a sample-level posterior probability of peak rather than a binary classification, allowing users to distinguish robust peaks (high P) from age-marginal ones (low P).

Intended scope of uncertainty propagation: chronological uncertainty as it flows through peak detection, including the threshold-chronology coupling. Different chronologies produce different smoothing window placements, which produce different local thresholds, which produce different peaks. P(peak) therefore reflects both age uncertainty and threshold variability propagated jointly through the algorithm. This is the desired behavior, not a limitation, and methods text describing CharAnalysis-with-uncertainty should be explicit about it.

Out of scope for the initial integration: propagating uncertainty further downstream into fire return interval (FRI) distributions. That cascade (chronology uncertainty -> peak uncertainty -> FRI uncertainty) is a follow-up, with the peak-level posterior matrix providing a clean handoff to whatever the downstream FRI implementation looks like.

(Section added 2026-05-06; revise here as scope evolves.)

---

## 2026-05-06 : architecture decision : source-agnostic core with rplum adapter (Option 3)

### Decision
Adopt a hybrid architecture in which the core peak-identification function takes a posterior age-depth matrix as input and is agnostic to its source.

```r
char_peak_id_uncertain(charcoal_data, age_posterior, params, ...)
# age_posterior: matrix where rows = depths, columns = posterior draws
```

Provide a thin adapter for users who want to drive the workflow from rplum:

```r
age_posterior <- extract_age_posterior_rplum(rplum_object, depths)
```

And optionally a convenience wrapper that chains rplum and CharAnalysis for the common one-shot case:

```r
char_pipeline_with_rplum(rplum_object, charcoal_data, params, n_draws = 1000)
```

### Rationale
Two alternative architectures were considered and rejected:

1. **CharAnalysis wraps rplum end to end.** Adds rplum and rbacon as hard dependencies (both compile C++; rbacon depends on the external Bacon binary), which substantially raises install friction. Couples release cycles. Locks out users who use BChron, Clam, OxCal, or hand-built chronologies. Forces a poor iteration loop because every CharAnalysis parameter tweak triggers a rerun of millions of rplum MCMC iterations. Buries critical rplum priors (flux, supported 210Pb, accumulation rate, memory) behind a CharAnalysis API where they should not live.

2. **Manual chronology export only, no rplum integration.** Reasonable but incomplete. Forces every user into a two-step workflow with no convenience path.

Option 3 dominates: the core function works with any age-model output (including Phil's existing Matlab tool MCAgeDepth.m, which generates Monte Carlo age-depth iterations and predates rplum), the adapter handles the common rplum case, and `rplum` stays a Suggests dependency rather than Imports. This pattern is standard in R; tidymodels and brms use the same shape for source-specific extensions.

### Implications for the kickoff design questions
- Question 3 (additive vs. breaking API) is now resolved via the new function names. Existing `char_peak_id()` is untouched, so the integration is additive. Version bump at merge time will be 2.x minor unless something else in the implementation forces a breaking change.
- Questions 1, 2, 4, and 5 remain open.

### New considerations to track
- `extract_age_posterior_rplum()` needs a clear specification of (a) which subset of MCMC draws to extract (after burn-in and thinning, per rplum's defaults or user override) and (b) how to handle the depth grid. rplum uses 1 cm internal sections by default; charcoal records often have irregular sample depths. Standard answer: interpolate rplum's posterior age-depth surface onto the charcoal sample depths.
- A test using output from Phil's MCAgeDepth.m (or a synthetic posterior) will validate the source-agnostic claim before the rplum adapter is written. This is a good early task because it forces the input contract to be defined first.

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
