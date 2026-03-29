# CharAnalysis 2.0 — MATLAB Implementation

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

## Directory Contents

**Source code**
- `CharAnalysis.m` — main entry point; run this to start the program
- All other `.m` files — supporting functions called by `CharAnalysis.m`

**Version comparison and validation**
- `z_Compare_CharAnalysis_V1_V2.m` — script for comparing V2.0 outputs against
  V1.1 reference results; recommended for users migrating from Version 1.1
- `CO_charResults_V_1_1.csv` and `CH10_charResults_V_1_1.csv` — reference
  outputs for the Code Lake and Chickaree Lake example datasets, when using *CharAnalysis Version 1.1*,  used by the
  comparison script
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
