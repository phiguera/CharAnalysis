# _CharAnalysis_ 2.0 — MATLAB Implementation

This directory contains the Version 2.0 MATLAB source code for *CharAnalysis*.
For full documentation, installation instructions, and citations, see the
[User's Guide](https://github.com/phiguera/CharAnalysis/blob/master/CharAnalysis_UsersGuide.md)
and the [main README](https://github.com/phiguera/CharAnalysis/blob/master/README.md).

## Getting Started

1. Add this folder to your MATLAB search path:
```matlab
   addpath 'path/to/CharAnalysis_2_0_MATLAB'
   savepath
```
2. Type `CharAnalysis` in the MATLAB Command Window and follow the prompts.
3. Use the template files in
   [`DataTemplates_and_Examples/`](https://github.com/phiguera/CharAnalysis/tree/master/DataTemplates_and_Examples)
   to set up input files for your own data.

## Calling Conventions
```matlab
% Standard run — all figures, existing behavior
CharAnalysis('mysite_charParams.csv')

% Return results to the MATLAB workspace
results = CharAnalysis('mysite_charParams.csv')

% Call any individual figure directly
CharPlotFig7_ContinuousFireHistory(results)

% Interactive figure menu
CharAnalysis('mysite_charParams.csv', 'modular')

% Programmatic figure selection (e.g. Figures 3 and 7 only)
CharAnalysis('mysite_charParams.csv', 'modular', [3 7])

% Programmatic selection with automatic save
CharAnalysis('mysite.csv', 'modular', [3 7], true)

% Run analysis and save data only, no figures generated
CharAnalysis('mysite.csv', 'resultsOnly')
```

See Section 3.1 of the User's Guide for full details on the modular figure interface.

## Directory Contents

**Entry point**
- `CharAnalysis.m` — main entry point; the only `.m` file at this level. Run this to start the program. All supporting functions are in `src/`.

**Source code**
- `src/` — all supporting `.m` files, including analytical functions, figure functions, and helper utilities. Added to the MATLAB path automatically when `CharAnalysis.m` runs; no manual path configuration of this subfolder is needed.

**Version comparison and validation**
- `z_Compare_CharAnalysis_V1_V2.m` — script for comparing V2.0 outputs against
  V1.1 reference results; recommended for users migrating from Version 1.1
- `CO_charResults_V_1_1.csv` and `CH10_charResults_V_1_1.csv` — Code Lake and
  Chickaree Lake results produced by *CharAnalysis* Version 1.1, used as
  reference outputs by the comparison script
- `z_Validation_report_CharAnalysis_V_2.0.md` — documents V2.0 validation
  against V1.1 outputs

**Example datasets** *(for use with the comparison script)*
- `CO_charData.csv`, `CO_charParams.csv`, `CO_charResults.csv` — Code Lake
- `CH10_charData.csv`, `CH10_charParams.csv`, `CH10_charResults.csv` — Chickaree Lake

## Version Comparison

Users migrating from Version 1.1 are encouraged to run
`z_Compare_CharAnalysis_V1_V2.m` on the bundled example datasets to verify that
V2.0 results are consistent with V1.1 outputs. The analytical methods are
unchanged between versions; any differences should be within numerical tolerance.
See `z_Validation_report_CharAnalysis_V_2.0.md` for full documentation of the
validation results.
