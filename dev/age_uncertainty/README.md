# dev/age_uncertainty/

Development workspace for integrating chronological uncertainty into CharAnalysis peak detection using the R package `rplum`.

## Contents

- `lab_notebook.md` : running log of decisions, open questions, and next steps. Newest entries at top.
- (additional scripts, prototypes, and exploratory analyses added as the work progresses)

## Convention

The repo-root `dev/` directory is for exploratory scripts, prototypes, and notes that are version-controlled but not part of any released package or distributable artifact. It sits sibling to `CharAnalysis_2_0_R/` (the released R package) and the MATLAB version directories.

When this work is integrated into a released package version, the production code moves into the package source tree (`CharAnalysis_2_0_R/R/` or its successor), and `dev/age_uncertainty/` retains the historical reasoning trail.

## Branch

Active development branch: `feature/age-uncertainty` (off `dev`).
