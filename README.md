## *CharAnalysis*: Diagnostic and analytical tools for peak detection in sediment-charcoal records

**_CharAnalysis_** is a program for analyzing sediment-charcoal records when the goal is peak detection to reconstruct local fire history. Since its original development in the mid-2000s, the program has been used in dozens of published studies across six continents. The entire codebase is distributed and well commented — users are encouraged to look under the hood, understand what's going on, and modify the program to suit individual needs.

###### [![CC BY 4.0](https://licensebuttons.net/l/by/4.0/88x31.png)](https://creativecommons.org/licenses/by/4.0/) © 2004–2026
Philip Higuera\
Professor, Department of Ecosystem and Conservation Sciences\
University of Montana, Missoula, MT, USA\
https://www.umt.edu/people/phiguera

## Getting Started

**Option 1: Download and run locally in MATLAB**\
Download the entire *CharAnalysis* program as a .zip or tar.gz archive from the 
project website at https://phiguera.github.io/CharAnalysis/. Users familiar with 
GitHub can also download individual files or clone the repository directly at 
https://github.com/phiguera/CharAnalysis. Requires MATLAB R2019a or higher.

**Option 2: Download and run the standalone Windows application (Version 1.1)**\
A standalone Windows executable (.exe) is available for users without MATLAB. 
Note that this version (1.1) predates the Version 2.0 update. See the 
[standalone application readme](https://github.com/phiguera/CharAnalysis/blob/master/CharAnalysis_1_1_Windows/readme_CharAnalysis_standAlone.md) 
for download and installation instructions.

**Option 3: Try it online**\
Not sure if *CharAnalysis* is the right tool for your research? Click the badge 
below to run the program instantly in your browser on the bundled Code Lake 
example dataset — no installation required. A free MathWorks account is needed; 
university users can log in with their institutional email for full access. [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=phiguera/CharAnalysis&branch=main&file=CharAnalysis_2_0_MATLAB/CharAnalysis.m)


## Issues and troubleshooting

If you are having problems running *CharAnalysis*, please use the Issues tab on GitHub to 
look for other users who may have had the same problem, or to register a new 
issue. To register an issue you must have a GitHub account. Resolved issues from 
before April 2014, when the project was hosted on Google Code, are archived at 
http://code.google.com/p/charanalysis/.

## Version 2.0 (March 2026)

Version 2.0 is the first major update to *CharAnalysis*, addressing five areas of improvement:

1. **MATLAB Modernization** — eliminated legacy code patterns, vectorized inner loops, and removed deprecated function calls for compatibility with MATLAB R2019a and higher.
2. **Architecture Improvements** — separated computational logic from visualization, introduced formal parameter objects, and added robust input validation.
3. **Chronological Uncertainty** — integrated methods for incorporating chronological uncertainty into the characterization of fire events.
4. **Regional Synthesis** — added support for synthesizing peak identification across multiple records at regional scales.
5. **R Translation** — established a roadmap for a native R implementation using modern packages including `mclust` and `ggplot2`, to be quantitatively compared with the `tapas` R package (https://github.com/wfinsinger/tapas).

*Version 2.0 was developed with the assistance of Claude, an AI assistant by Anthropic. Claude assisted with code modernization, bug fixes, architecture redesign, and documentation. All code was reviewed and validated by the author against Version 1.1 reference outputs.*

## Understanding and citing the program

**The following citations provide the most thorough background for the analyses employed in _CharAnalysis_, including the assumptions inherent in charcoal peak identification using a local threshold. These may be downloaded from the author's web site, linked to above.**

> Kelly, R.F., P.E. Higuera, C.M. Barrett, and F.S. Hu. 2011. A signal-to-noise index to quantify the potential for peak detection in sediment-charcoal records. _Quaternary Research_ 75: 11-17.

> Higuera, P.E., D.G. Gavin, P.J. Bartlein and D.J. Hallett. 2010. Peak detection in sediment-charcoal records: impacts of alternative data analysis methods on fire-history interpretations. _International Journal of Wildland Fire_ 19: 996-1014. 

**Please use the following citation (also available via the author's web site) if you use the _CharAnalysis_ program in a publication, and please note that it is freely available at https://phiguera.github.io/CharAnalysis/**

> Higuera, P.E., L.B. Brubaker, P.M. Anderson, F.S. Hu, T.A. Brown. 2009. Vegetation mediated the impacts of postglacial climatic change on fire regimes in the south-central Brooks  Range, Alaska. _Ecological Monographs_ 79: 201-219. 

**Other studies using _CharAnalysis_ (incomplete list, through 2014)**

> Dunnette, P.V., P.E. Higuera, K.K. McLauchlan, K.M. Derr, C.E. Briles, and M.H. Keefe. 2014. Biogeochemical impacts of wildfires over four millennia in a Rocky Mountain subalpine watershed. _New Phytologist_ 203: 900-912.

> Brossier, B., F. Oris, W. Finsinger, H. Asselin, Y. Bergeron, and A. A. Ali. 2014. Using tree-ring records to calibrate peak detection in fire reconstructions based on sedimentary charcoal records. _The Holocene_ 24: 635–645.

> Kelly, R., M.L. Chipman, P.E. Higuera, I. Stefanova, L.B. Brubaker, and F.S. Hu. 2013. Recent burning of boreal forests exceeds fire regime limits of the past 10,000 years. _Proceedings of the National Academy of Sciences_ 110:13055-13060.

> Mustaphi, C.J.C. and M.F.J. Pisaric. 2013. Varying influence of climate and aspect as controls of montane forest fire regimes during the late Holocene, south-eastern British Columbia, Canada. _Journal of Biogeography_ 40:1983-1996.

> Blarquez, O., M. P. Girardin, B. Leys, A. A. Ali, J. C. Aleman, Y. Bergeron, and C. Carcaillet. 2013. Paleofire reconstruction based on an ensemble-member strategy applied to sedimentary charcoal. _Geophysical Research Letters_ 40:2667-2672.

> Barrett, C.M., R.F. Kelly, P.E. Higuera, and F.S. Hu. 2013. Climatic and land cover influences on the spatiotemporal dynamics of Holocene boreal fire regimes. _Ecology_ 92:389-402.

> Ali, A. A., O. Blarquez, M. P. Girardin, C. Hely, F. Tinquaut, A. El Guellab, V. Valsecchi, A. Terrier, L. Bremond, A. Genries, S. Gauthier, and Y. Bergeron. 2012. Control of the multimillennial wildfire size in boreal North America by spring climatic conditions. _Proceedings of the National Academy of Sciences_ 109:20966-20970.

> Marlon, J.R., P.J. Bartlein, D.G. Gavin, C.J. Long, R.S. Anderson, C.E. Briles, K.J. Brown, D. Colombaroli, D.J. Hallett, M.J. Power, E.A. Scharf, and M.K. Walsh. 2012. Long-term perspective on wildfires in the western USA. _Proceedings of the National Academy of Sciences_ doi:10.1073/pnas.1112839109 

> van Bellen, S., M. Garneau, A.A. Ali, and Y. Bergeron. 2012. Did fires drive Holocene carbon sequestration in boreal ombrotrophic peatlands of eastern Canada? _Quaternary Research_ 78:50-59.

> Higuera, P.E., M.L. Chipman, J.L. Barnes, J.L., M.A. Urban, and F.S. Hu. 2011. Variability of tundra fire regimes in Arctic Alaska: millennial scale patterns and ecological implications. _Ecological Applications_ 21: 3211-3226.
     
> Briles, C.E., C. Whitlock, C.N. Skinner, and J. Mohr. 2011. Postglacial forest development on different substrates in the Klamath mountains, northern California, USA. _Ecology_ 92: 590-601. 
     
> Whitlock, C., C.E. Briles, M.C. Fernandez, J. Gage. 2011. Holocene vegetation, fire, and climate history of the Sawtooth Range, central Idaho, U.S.A. _Quaternary Research_ 75: 114-124.
     
> Higuera, P.E., C. Whitlock, and J. Gage. 2011. Linking tree-ring and sediment-charcoal records to reconstruct fire occurrence and area burned in subalpine forests of Yellowstone National Park, U.S.A. _The Holocene_ 21: 327-341. 

> Hu, F.S., P.E. Higuera, J.E. Walsh, W.L. Chapman, P.A. Duffy, L.B. Brubaker and M.L. Chipman. 2010. Tundra burning in Alaska: linkages to climatic change and sea-ice retreat. _Journal of Geophysical Research - Biogeosciences_ 115, G04002, doi:10.1028/2009JG001270

> Walsh, M. K., C. Whitlock, and P. J. Bartlein. 2010. 1200 years of fire and vegetation history in the Willamette Valley, Oregon and Washington, reconstructed using high-resolution macroscopic charcoal and pollen analysis. Palaeogeography Palaeoclimatology Palaeoecology 297:273-289.

> Higuera, P.E., L.B. Brubaker, P.M. Anderson, F.S. Hu, T.A. Brown. 2009. Vegetation mediated the impacts of postglacial climatic change on fire regimes in the south-central Brooks  Range, Alaska. _Ecological Monographs_ 79: 201-219. 

> Chepstow-Lusty, A. J., M. R. Frogley, B. S. Bauer, M. J. Leng, K. P. Boessenkool, C. Carcaillet, A. A. Ali, and A. Gioda. 2009. Putting the rise of the Inca Empire within a climatic and land management context. _Climate of the Past_ 5:375-388.

> Marlon, J.R., P.J. Bartlein, M.K. Walsh, S.P. Harrison, K.J. Brown, M.E. Edwards, P.E. Higuera, M.J. Power, R.S. Anderson, C. Briles, A. Brunelle, C. Carcaillet, M. Daniels, F.S. Hu, M. Lavoie, C. Long, T. Minckley, P.J.H. Richard, S.L. Shafer, W. Tinner, C. Umbanhowar, and C. Whitlock. 2009. Wildfire responses to abrupt climate change in North America. _PNAS_ 106: 2519-2524. 

> Walsh, M. K., C. Whitlock, and P. J. Bartlein. 2008. A 14,300-year-long record of fire-vegetation-climate linkages at Battle Ground Lake, southwestern Washington. _Quaternary Research_ 70:251-264.

> Huerta, M. A., C. Whitlock, and J. Yale. 2009. Holocene vegetation-fire-climate linkages in northern Yellowstone National Park, USA. _Palaeogeography, Palaeoclimatology, Palaeoecology_ 271:170-181.

> Briles, C.E., C. Whitlock, P.J. Bartlein, and P.E. Higuera. 2008. Regional and local controls on postglacial vegetation and fire in the Siskiyou Mountains, northern California, USA. _Palaeogeography Palaeoclimatology Palaeoecology_. 265: 159-169. 

> Higuera, P. E., L. B. Brubaker, P. M. Anderson, T. A. Brown, A. T. Kennedy, and F. S. Hu. 2008. Frequent Fires in Ancient Shrub Tundra: Implications of Paleorecords for Arctic Environmental Change. _PLoS ONE_ 3:e0001744.

## Acknowledgments

Many features in _CharAnalysis_ are based on the programs CHAPS, by Patrick Bartlein (U of OR), and Charster, by Daniel Gavin (U of OR). The Gaussian mixture model used in _CharAnalysis_ was created by Charles Bouman (Purdue). _CharAnalysis_ was written in and compiled by Matlab 7.0 with resources from the University of Washington, Montana State University, the University of Illinois, the University of Idaho, and the University of Montana. Development of the program has benefited greatly from discussions with and testing by members of the Whitlock Paleoecology Lab at Montana State University, Dan Gavin, and Ryan Kelly. Please see the _CharAnalysis_ User's Guide for more details.

## Disclaimer

THIS SOFTWARE PROGRAM AND DOCUMENTATION ARE PROVIDED “AS IS” AND WITHOUT WARRANTIES AS TO
PERFORMANCE. THE PROGRAM _CharAnalysis_ IS PROVIDED WITHOUT ANY EXPRESSED OR IMPLIED WARRANTIES WHATSOEVER. BECAUSE OF THE DIVERSITY OF CONDITIONS AND HARDWARE UNDER WHICH THE PROGRAM MAY BE USED, NO WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE IS OFFERED. THE USER IS ADVISED TO TEST THE PROGRAM THOROUGHLY BEFORE RELYING ON IT. THE USER MUST ASSUME THE ENTIRE RISK AND RESPONSIBILITY OF USING THIS PROGRAM.
