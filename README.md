## *CharAnalysis*: Diagnostic and analytical tools for peak detection and fire-history interpretations using high-resolution sediment-charcoal records

[![CC BY 4.0](https://licensebuttons.net/l/by/4.0/88x31.png)](https://creativecommons.org/licenses/by/4.0/) © 2004–2026\
Philip Higuera\
Professor, Department of Ecosystem and Conservation Sciences\
University of Montana, Missoula, MT, USA\
https://www.umt.edu/people/phiguera

***CharAnalysis*** is a program for analyzing sediment-charcoal records when the goal is peak detection to reconstruct local fire history. Since its original development in the mid-2000s, the program has been used in dozens of published studies to analyze sediment-charcoal records worldwide. The entire codebase is distributed and well commented — users are encouraged to look under the hood, understand what's going on, and modify the program to suit individual needs.

## Getting Started

**Option 1: Download and run locally in MATLAB** *(recommended)*\
Download the entire *CharAnalysis* program as a `.zip` or `tar.gz` archive from the
project website at https://phiguera.github.io/CharAnalysis/, or clone the repository
directly at https://github.com/phiguera/CharAnalysis. Requires MATLAB R2019a or
higher. No additional toolboxes are required. See the
[User's Guide](CharAnalysis_UsersGuide.md) for installation instructions.

**Option 2: Download and run the standalone Windows application (Version 1.1)**\
A standalone Windows executable (`.exe`) is available for users without a MATLAB
license. Note that this version predates the Version 2.0 update. See the
[standalone application readme](CharAnalysis_1_1_Windows/readme_CharAnalysis_standAlone.md)
for download and installation instructions.

**Option 3: Try it online**\
Not sure if *CharAnalysis* is the right tool for your research? Click the badge
below to run the program instantly in your browser on the bundled Code Lake example
dataset — no installation required. A free MathWorks account is needed; university
users can log in with their institutional email for full access. [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=phiguera/CharAnalysis&branch=main&file=CharAnalysis_2_0_MATLAB/CharAnalysis.m)

## Repository Contents
```
CharAnalysis_2_0_MATLAB/       MATLAB source code (Version 2.0)
CharAnalysis_1_1_MATLAB/       MATLAB source code (Version 1.1)
CharAnalysis_1_1_Windows/      Standalone Windows application (Version 1.1)
DataTemplates_and_Examples/    Template input files and Code Lake example dataset
```

## Documentation

The current User's Guide is available as
[CharAnalysis_UsersGuide.md](CharAnalysis_UsersGuide.md). It covers installation,
data input, parameter selection, and a full description of all analytical choices
and program output. The original guide (v0.9, January 2009) is retained for reference as
[ARCHIVED_UsersGuide_v_0.9_2009.pdf](ARCHIVED_UsersGuide_v_0.9_2009.pdf).

Planned future development is described in [ROADMAP.md](ROADMAP.md).

## Version 2.0 (March 2026)

Version 2.0 is the first major update to *CharAnalysis*, addressing two areas of improvement, with four additional areas planned (see [ROADMAP.md](ROADMAP.md)):

1. **MATLAB Modernization** — eliminated legacy code patterns, vectorized inner loops, and removed deprecated function calls for compatibility with MATLAB R2019a and higher.
2. **Modular Figure Interface** — restructured the codebase so that each output figure (Figures 3–9) is implemented as a standalone function callable independently using the results struct returned by `CharAnalysis`. All `.m` files except `CharAnalysis.m` now reside in a `src/` subfolder, added to the MATLAB path automatically at startup. A new `'modular'` run mode supports interactive and programmatic figure selection:
```matlab
% Standard run — all figures, existing behavior unchanged
CharAnalysis('mysite.csv')

% Return results to workspace for further use
results = CharAnalysis('mysite.csv')

% Call any individual figure directly
CharPlotFig7_ContinuousFireHistory(results)

% Interactive figure menu — select figures from a command-window prompt
CharAnalysis('mysite.csv', 'modular')

% Programmatic figure selection — generate figures 3 and 7 only
CharAnalysis('mysite.csv', 'modular', [3 7])

% Programmatic selection with automatic save
CharAnalysis('mysite.csv', 'modular', [3 7], true)
```

The modular architecture also lays the groundwork for the planned R translation by cleanly separating computation from visualization.

*Version 2.0 was developed with the assistance of Claude, an AI assistant by Anthropic. Claude assisted with code modernization, bug fixes, architecture redesign, and documentation. All code was reviewed and validated by the author against Version 1.1 reference outputs.*

## Citation

If you use *CharAnalysis* in a publication, please cite Higuera et al. (2009),
the first study to apply the core analytical tools implemented in *CharAnalysis*.
If you used Version 2.0 specifically, please also cite the software:

[Higuera, P.E., L.B. Brubaker, P.M. Anderson, F.S. Hu, and T.A. Brown. 2009.
Vegetation mediated the impacts of postglacial climate change on fire regimes
in the south-central Brooks Range, Alaska. *Ecological Monographs*
79:201–219.](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/07-2019.1)

[Higuera, P.E. 2026. *CharAnalysis*: Diagnostic and analytical tools for peak detection and fire-history interpretations using high-resolution sediment-charcoal records (Version 2.0).
Zenodo.](https://doi.org/10.5281/zenodo.19304064)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19304064.svg)](https://doi.org/10.5281/zenodo.19304064)

## Issues and Troubleshooting

If you encounter problems running *CharAnalysis*, please use the Issues tab at
https://github.com/phiguera/CharAnalysis/issues to search for existing reports or
register a new issue. A GitHub account is required. Resolved issues from before
April 2014, when the project was hosted on Google Code, are archived at
http://code.google.com/p/charanalysis/.

## Disclaimer

THIS SOFTWARE PROGRAM AND DOCUMENTATION ARE PROVIDED "AS IS" AND WITHOUT
WARRANTIES AS TO PERFORMANCE. THE PROGRAM *CharAnalysis* IS PROVIDED WITHOUT ANY
EXPRESSED OR IMPLIED WARRANTIES WHATSOEVER. BECAUSE OF THE DIVERSITY OF CONDITIONS
AND HARDWARE UNDER WHICH THE PROGRAM MAY BE USED, NO WARRANTY OF FITNESS FOR A
PARTICULAR PURPOSE IS OFFERED. THE USER IS ADVISED TO TEST THE PROGRAM THOROUGHLY
BEFORE RELYING ON IT. THE USER MUST ASSUME THE ENTIRE RISK AND RESPONSIBILITY OF
USING THIS PROGRAM.
