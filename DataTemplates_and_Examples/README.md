# DataTemplates and Examples

This directory contains everything needed to get started with *CharAnalysis*:
two worked example datasets you can run immediately, and blank template files
for setting up your own analysis. For full documentation see the
[User's Guide](https://github.com/phiguera/CharAnalysis/blob/master/CharAnalysis_UsersGuide.md).

## Example Datasets

Two example datasets are provided. Both can be run immediately using the
instructions in the User's Guide.

**Code Lake, Alaska** (`CO_` files)
- `CO_charData.csv` — charcoal data input
- `CO_charParams.csv` — analysis parameter selections
- `CO_charResults.csv` — expected output for reference

Data source: Higuera, P.E., L.B. Brubaker, P.M. Anderson, F.S. Hu, and T.A.
Brown. 2009. Vegetation mediated the impacts of postglacial climate change on
fire regimes in the south-central Brooks Range, Alaska. *Ecological Monographs*
79:201–219. https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/07-2019.1

---

**Chickaree Lake, Colorado** (`CH10_` files)
- `CH10_charData.csv` — charcoal data input
- `CH10_charParams.csv` — analysis parameter selections
- `CH10_charResults.csv` — expected output for reference

Data source: Dunnette, P.V., P.E. Higuera, K.K. McLauchlan, K.M. Derr, C.E.
Briles, and M.H. Keefe. 2014. Biogeochemical impacts of wildfires over four
millennia in a Rocky Mountain subalpine watershed. *New Phytologist* 203:900–912.
https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.12828

---

These same datasets are used in `CharAnalysis_2_0_MATLAB/` alongside Version 1.1
reference outputs to support version comparison. See the
[CharAnalysis_2_0_MATLAB README](https://github.com/phiguera/CharAnalysis/blob/master/CharAnalysis_2_0_MATLAB/README.md)
for details.

## Template Files

To analyze your own data, start with these blank template files:

- `template_charData.csv` — enter your charcoal data here (depths, ages, volumes,
  counts)
- `template_charParams.csv` — enter your parameter selections here

Copy and rename each file to identify your site (e.g. `MyLake_charData.csv`),
then follow the data input instructions in Part II of the
[User's Guide](https://github.com/phiguera/CharAnalysis/blob/master/CharAnalysis_UsersGuide.md).
