## _CharAnalysis_ 2.0: Modernization and R-Translation

**Note: This software update and the development of Version 2.0 are being conducted with the assistance of AI tools, including Gemini and Claude.**

_CharAnalysis_ is a program for analyzing sediment-charcoal records, when the goal is peak detection to reconstruct "local" fire history. Originally developed in the mid-2000s, the program is currently undergoing a significant update to Version 2.0. This evolution aims to modernize the existing MATLAB codebase—eliminating legacy patterns and improving computational efficiency—while establishing a formal bridge to the R programming environment.

###### (c) 2004-2026

Philip Higuera, Professor
Department of Ecosystem and Conservation Sciences 
University of Montana
Missoula, MT, USA
https://www.umt.edu/people/phiguera

## The Road to Version 2.0
The development of _CharAnalysis_ 2.0 is focused on several core pillars of improvement to better serve the paleoecological community:

1. **MATLAB Modernization:** We are eliminating legacy code patterns, vectorizing inner loops for faster processing, and removing deprecated function calls to ensure compatibility with modern MATLAB versions.
2. **Architecture Improvements:** To increase stability, we are separating computational logic from visualization, introducing formal parameter objects for better record-keeping, and adding robust input validation.
3. **Chronological Uncertainty:** A primary goal of this update is to integrate methods for incorporating chronological uncertainty directly into the characterization of fire events.
4. **Regional Synthesis:** Version 2.0 will support new methods for synthesizing peak identification from multiple sediment-charcoal records at regional spatial scales.
5. **R Translation Strategy:** Recognizing the dominance of R in paleoecology, Version 2.0 serves as a module-by-module roadmap for a native R implementation, utilizing modern packages like `mclust` and `ggplot2`.

## Using this site 
### Downloads
Download the current stable _CharAnalysis_ program as a .zip or tar.gz archive by clicking on the appropriate icons on this page. Development branches for Version 2.0 and the R-translation modules will be posted to the GitHub page as they reach beta status.

### Updates and Comments
Progress on the Version 2.0 action plan is described in the Wiki tab. We encourage users to "look under the hood" of the new modules to provide feedback on the modernization of specific diagnostic and analytical tools.

### Issues
If you encounter bugs while running the legacy code or have suggestions for the Version 2.0 architecture, please use the "Issues" tab. This collaborative feedback is vital as we transition to a more efficient, multi-language platform.

## Acknowledgments
Many features in _CharAnalysis_ are based on the programs CHAPS, by Patrick Bartlein (U of OR), and Charster, by Daniel Gavin (U of OR). The Gaussian mixture model used in _CharAnalysis_ was created by Charles Bouman (Purdue). 

_CharAnalysis_ was developed with resources from the University of Washington, Montana State University, the University of Illinois, the University of Idaho, and the University of Montana. The Version 2.0 modernization efforts are driven by the need for increased portability and performance in modern research workflows. Please see the _CharAnalysis_ Modernization Report for the full technical roadmap.

## Disclaimer
THIS SOFTWARE PROGRAM AND DOCUMENTATION ARE PROVIDED “AS IS” AND WITHOUT WARRANTIES AS TO PERFORMANCE. THE DEVELOPMENT OF VERSION 2.0 IS ONGOING; USERS ARE ADVISED TO TEST NEW MODULES THOROUGHLY AGAINST THE LEGACY VERSION BEFORE RELYING ON THEM FOR PUBLICATION. THE USER MUST ASSUME THE ENTIRE RISK AND RESPONSIBILITY OF USING THIS PROGRAM.

THE USE OF THE SOFTWARE DOWNLOADED FROM THIS WEBSITE IS DONE AT YOUR OWN RISK AND WITH AGREEMENT THAT YOU ARE SOLELY RESPONSIBLE FOR ANY DAMAGE TO YOUR COMPUTER SYSTEM OR LOSS OF DATA THAT RESULTS FROM SUCH ACTIVITIES.
