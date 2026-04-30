## Resubmission (2.0.3)

Patch release with substantive vignette improvements and several
factual corrections to the worked example. No changes to analytical
behaviour, to any function signature, or to the package's public API.

- The package vignette now live-renders all worked-example figures
  from the bundled Code Lake dataset; previously the figure chunks
  were not evaluated. The chunks are gated on the availability of
  `ggplot2`, `patchwork`, and `ggtext` (all in `Suggests`), so they
  skip cleanly on any test environment lacking those packages.
- Two new bundled files in `inst/validation/`:
  `CO_compensated_charParams.csv` and `CO_compensated_charData.csv`.
  The compensated parameter file is used by the worked example to
  produce results aligned with the published interpretation in
  Higuera et al. (2009). The strict 1-to-1 reference files
  (`CO_charParams.csv`, `CO_charData.csv`) are unchanged.
- Vignette factual corrections: corrected geographic attribution of
  Code Lake (Alaska, not Colorado); corrected the smoothing-method
  label for Code Lake in the comparison table; fixed a stale
  `system.file()` path that referenced `inst/extdata` rather than
  `inst/validation`.
- Version bumped from 2.0.2 to 2.0.3.

## Prior resubmission (2.0.2)

Changes in response to the third round of CRAN reviewer feedback
(Benjamin Altmann, 2026-04-27):

- `Title` shortened from 131 to 62 characters: "Peak Detection and
  Fire History from Sediment-Charcoal Records".
- `char_write_results()` no longer defaults `out_dir` to `"."`. The
  argument is now required; the function `stop()`s with a message
  pointing users to `tempdir()` if it is omitted.
- `char_plot_all()` no longer defaults `out_dir` to `"."`. The argument
  is required when `save = TRUE` and ignored otherwise. The default
  signature is now `out_dir = NULL`.
- The package vignette previously wrote its example output to a
  `Results/` directory in the user's working directory. Both calls
  (`char_plot_all()` and `char_write_results()`) now write to
  `tempdir()` with a note that real users would substitute their own
  path.
- Version bumped from 2.0.1 to 2.0.2.

While preparing this resubmission we also corrected a figure-rendering
bug in three plot functions (`char_plot_peaks()`,
`char_plot_fire_history()`, `char_plot_zones()`) whose axis labels were
not displaying correctly under the ggtext-based theme introduced in
2.0.1. No analytical behaviour changes.

## Prior resubmission (2.0.1)

Changes in response to the second round of CRAN reviewer feedback
(2026-04-22):

- DESCRIPTION: removed the ` .` paragraph separators that rendered as
  "program. . Since its" in CRAN metadata. The Description has been
  condensed to a single paragraph.
- Added `@export` to `char_parameters()` and updated NAMESPACE. The
  function's help page previously documented an example for an
  unexported function.
- Replaced all `\dontrun{}` wrappers with `\donttest{}`, or unwrapped
  them entirely when the example runs in under 5 seconds. The
  `char_parameters()` example now runs on every check. All examples
  reference the bundled validation dataset via `system.file()` and write
  any output to `tempdir()`.
- Version bumped from 2.0.0 to 2.0.1.

## Prior resubmission (2.0.0)

Changes in response to the first round of CRAN reviewer feedback
(Uwe Ligges, 2026-04-16):

- Software name 'CharAnalysis' is now single-quoted in the Description
  field (two occurrences).
- DOIs added to all four citations in Description:
    Higuera et al. (2008) <doi:10.1371/journal.pone.0001744>
    Higuera et al. (2009) <doi:10.1890/07-2019.1>
    Higuera et al. (2010) <doi:10.1071/WF09134>
    Kelly et al. (2011)   <doi:10.1016/j.yqres.2010.07.011>

(Note: for 2.0.1 the Higuera 2008 and Kelly 2011 citations were removed
from DESCRIPTION as part of condensing the paragraph; both remain
referenced in the package vignette and User's Guide.)

## R CMD check results

Win-builder R-devel: 0 errors | 0 warnings | 1 NOTE.

Local R CMD check on Windows: 0 errors | 0 warnings | 0 notes.

The single NOTE on win-builder is "Days since last update: 1", flagged
by the CRAN incoming feasibility check because v2.0.2 was accepted on
2026-04-28 and this v2.0.3 patch was prepared the following day.

This rapid resubmission is intentional. v2.0.3 corrects two
factual errors in the package vignette (the geographic attribution
of the Code Lake worked-example dataset and the smoothing-method
label in the validation table) and substantively improves the
worked example by live-rendering all output figures. None of these
changes alter analytical behavior or any function signature, but
the vignette corrections were significant enough that we preferred
not to leave them in the version users would adopt over the coming
weeks.

## Test environments

- Local: Windows 11 x64 (build 26200), R 4.6.0 (2026-04-24 ucrt)
- win-builder: R Under development (unstable) (2026-04-28 r89972 ucrt),
  Windows Server 2022 x64