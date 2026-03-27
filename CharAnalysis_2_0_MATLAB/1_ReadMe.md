# CharAnalysis 2.0: Modernization & Functional Expansion

**Note: This software update and the development of Version 2.0 are being conducted with the assistance of AI tools, including Gemini and Claude.**

This directory contains the development version of _CharAnalysis_ 2.0. This update represents a major overhaul of the mid-2000s MATLAB codebase, focused on increasing computational efficiency, improving software architecture, and expanding the analytical scope to meet the needs of modern paleoecological research.

## Core Development Goals

### 1. MATLAB Modernization & Efficiency
* **Legacy Cleanup:** Removing deprecated function calls and modernizing syntax to ensure long-term compatibility with current and future MATLAB releases.
* **Vectorization:** Re-engineering inner loops into vectorized operations to significantly reduce processing time for high-resolution and long-duration records.
* **Architecture:** Separating computational logic from visualization routines to facilitate "headless" batch processing and more robust data handling.

### 2. New Analytical Capabilities
* **Chronological Uncertainty:** Integration of methods to formally incorporate age-model uncertainty into the characterization of fire events.
* **Regional Synthesis:** Development of standardized modules for synthesizing peak identification results across multiple records at regional spatial scales.

### 3. R-Translation Roadmap
* A primary objective of Version 2.0 is to serve as a module-by-module blueprint for a native R implementation. 
* **Redundancy Reduction:** This work is being integrated with—or quantitatively compared to—the charcoal peak analysis tools in the `tapas` R package ([https://github.com/wfinsinger/tapas](https://github.com/wfinsinger/tapas)).

## Directory Structure (Ongoing)
* `/Functions`: Modernized MATLAB `.m` files with improved input validation.
* `/Documentation`: Updated User's Guide and the *CharAnalysis Modernization Report*.
* `/Testing`: Regression test scripts to ensure V 2.0 outputs remain consistent with legacy V 1.x results.

## Contributing & Feedback
As this is an active development branch, users are encouraged to:
1. Compare V 2.0 outputs against legacy results.
2. Report any dimension mismatches or syntax errors via the **Issues** tab on GitHub.
3. Review the *Modernization Report* for detailed technical specifications.

## License & Citation
Usage of this development version should cite the original methodology (Higuera et al., 2009) and acknowledge the use of the Version 2.0 development branch.

**Disclaimer:** This version is under active development. Users should assume the risk of using beta-status code and are advised to verify all critical results against the stable Version 1.1.

---
**Philip Higuera** 

Department of Ecosystem and Conservation Sciences  

University of Montana  

[https://www.umt.edu/people/phiguera](https://www.umt.edu/people/phiguera)