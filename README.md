## _CharAnalysis_: Diagnostic and analytical tools for sediment-charcoal analysis

**Note: Current software updates and the development of Version 2.0 are being conducted with the assistance of AI tools, including Gemini and Claude.**

_CharAnalysis_ is a program for analyzing sediment-charcoal records, when the goal is peak detection to reconstruct "local" fire history. Originally developed in the mid-2000s, the program is currently undergoing a significant update to **Version 2.0**. This evolution aims to modernize the existing MATLAB codebase—eliminating legacy patterns and improving computational efficiency—while establishing a formal bridge to the R programming environment.

###### (c) 2004-2026

Philip Higuera, Professor, Department of Ecosystem and Conservation Sciences, University of Montana, Missoula, MT, USA

https://www.umt.edu/people/phiguera


## The Road to Version 2.0
The development of _CharAnalysis_ 2.0 is focused on several core pillars of improvement to better serve the paleoecological community:

1.  **MATLAB Modernization:** Eliminating legacy code patterns, vectorizing inner loops for faster processing, and removing deprecated function calls for compatibility with modern MATLAB versions.
2.  **Architecture Improvements:** Separating computational logic from visualization, introducing formal parameter objects for better record-keeping, and adding robust input validation.
3.  **Chronological Uncertainty:** Integrating methods for incorporating chronological uncertainty directly into the characterization of fire events.
4.  **Regional Synthesis:** Supporting new methods for synthesizing peak identification from multiple sediment-charcoal records at regional spatial scales.
5.  **R Translation Strategy:** Establishing a roadmap for a native R implementation, utilizing modern packages like `mclust` and `ggplot2`.


## Using this site 

### Downloads
Download the current stable _CharAnalysis_ program (V 1.1) as a .zip or tar.gz archive by clicking on the appropriate icons on this page. Alternatively, download individual files by visiting the GitHub page. Development branches for Version 2.0 and the R-translation modules will be posted as they reach beta status.

### Updates and Comments
Updates regarding the V 2.0 transition and legacy maintenance are described in the **Wiki** tab. We encourage users to "look under the hood" of the new modules to provide feedback on the modernization of specific diagnostic and analytical tools.

### Issues
If you are having problems running _CharAnalysis_ or have suggestions for the Version 2.0 architecture, please use the **Issues** tab to (1) search for similar problems/solutions, or (2) post a new issue.


## Getting Started with _CharAnalysis_ (V 1.x)
To use the current version of the program, follow these general steps. Detailed instructions are provided in the _CharAnalysis_ User's Guide (included in the download).

1.  **Prepare Input Data:** Create a `.xls`, `.xlsx`, or `.csv` file following the format specified in the User's Guide. This includes columns for depth, thickness, volume, and charcoal counts.
2.  **Set Parameters:** Modify the `CharParameters.m` file (or the input spreadsheet) to define the site name, interpolation resolution, and smoothing windows for background CHAR.
3.  **Run Analysis:** Execute `CharAnalysis.m` in MATLAB. The program will perform data pretreatment, background charcoal estimation, and peak identification.
4.  **Review Diagnostics:** Examine the generated figures (Figures 1-9) to evaluate the reliability of peak detection, including the Signal-to-Noise Index (SNI).


## Acknowledgments
Many features in _CharAnalysis_ are based on the programs CHAPS, by Patrick Bartlein (U of OR), and Charster, by Daniel Gavin (U of OR). The Gaussian mixture model used in _CharAnalysis_ was created by Charles Bouman (Purdue). 

_CharAnalysis_ was developed with resources from the University of Washington, Montana State University, the University of Illinois, the University of Idaho, and the University of Montana. Development has benefited from testing by members of the Whitlock Paleoecology Lab, Dan Gavin, and Ryan Kelly. Please see the _CharAnalysis_ Modernization Report for the full technical roadmap of the V 2.0 update.


## Disclaimer
THIS SOFTWARE PROGRAM AND DOCUMENTATION ARE PROVIDED “AS IS” AND WITHOUT WARRANTIES AS TO PERFORMANCE. THE DEVELOPMENT OF VERSION 2.0 IS ONGOING; USERS ARE ADVISED TO TEST NEW MODULES THOROUGHLY AGAINST THE LEGACY VERSION BEFORE RELYING ON THEM FOR PUBLICATION. 

THE USE OF THE SOFTWARE DOWNLOADED FROM THIS WEBSITE IS DONE AT YOUR OWN RISK AND WITH AGREEMENT THAT YOU ARE SOLELY RESPONSIBLE FOR ANY DAMAGE TO YOUR COMPUTER SYSTEM OR LOSS OF DATA THAT RESULTS FROM SUCH ACTIVITIES.
