## Flux Clamacitor 2.0.0

This release accompanies the Addition and Correction to Hartland et al. (2026), *Environ. Sci. Technol.* 60 (38), 27204−27214, doi:10.1021/acs.est.6c03421.

The calculations now sit in a small Python package, `clamacitor`, with two settings. `legacy` reproduces the published values from the published workbook, and `correction` gives the corrected values, in which Ca at each site and date is the median of the filtered and unfiltered ICP-MS replicates and the shell bulk density is 0.632 g cm⁻³ from the corrected specimen weighings. For the study period (17 October 2024 to 23 July 2025), Ca-based removal is 9.1 ± 7.8 t CaCO₃ d⁻¹ and the combined density 550 ± 496 individuals m⁻², against a survey mean of 557 ± 435 individuals m⁻².

The replicate ICP-MS data and the corrected specimen measurements are in `data/correction/`, and `scripts/30_correction_figures.py` redraws Figures 3 and 4 and Table 1. Input templates and a lake configuration file allow the method to be applied elsewhere. The browser tool is unchanged.
