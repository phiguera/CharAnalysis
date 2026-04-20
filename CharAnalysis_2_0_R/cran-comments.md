## Resubmission

Changes in response to CRAN reviewer feedback (Uwe Ligges, 2026-04-16):

- Software name 'CharAnalysis' is now single-quoted in the Description
  field (two occurrences).
- DOIs added to all four citations in Description:
    Higuera et al. (2008) <doi:10.1371/journal.pone.0001744>
    Higuera et al. (2009) <doi:10.1890/07-2019.1>
    Higuera et al. (2010) <doi:10.1071/WF09134>
    Kelly et al. (2011)   <doi:10.1016/j.yqres.2010.07.011>

## R CMD check results

0 errors | 0 warnings | 1 NOTE

The single NOTE is:

  checking for future file timestamps ... NOTE
  unable to verify current time

This is a transient network connectivity issue during the check (R is
unable to reach a time-verification server). It is not related to the
package and is expected on restricted or firewalled machines.

## Spell-check false positives in DESCRIPTION

The CRAN incoming feasibility check on win-builder flagged the following
as possibly misspelled:

  CharAnalysis — the name of this package (now single-quoted)
  Higuera      — author surname in citations
  et, al       — parts of "et al." in citations

## Test environments

- Local: Windows 11 x64 (build 26200), R 4.5.3 (2026-03-11 ucrt)
- win-builder: R-devel (4.6.0 beta, 2026-04-15 r89885 ucrt),
  Windows Server 2022 x64
