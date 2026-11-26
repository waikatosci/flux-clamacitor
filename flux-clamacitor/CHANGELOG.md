# Changelog

## 2.0.0

Accompanies the Addition and Correction to Hartland et al. (2026), *Environ. Sci. Technol.* 60 (38), 27204−27214.

- New `clamacitor` Python package, with a `legacy` setting that reproduces the published values and a `correction` setting that gives the corrected values. Both are checked by the test suite.
- Ca at each site and date is taken as the median of the filtered and unfiltered ICP-MS replicates (four measurements), which removes an instrumental artefact affecting about one sample in ten. The replicate data are in `data/correction/`.
- Shell bulk density is total dry mass divided by total ellipsoid volume for the 29 specimens, using weighings corrected from the laboratory notes: 0.632 g cm⁻³ (bootstrap 95% interval 0.57–0.70), in place of 1.15.
- `scripts/30_correction_figures.py` redraws Figures 3 and 4 and Table 1 with the corrected values, using the published figure scripts.
- Lake settings in YAML, input templates and a data dictionary for applying the method to other lakes, and a command line (`python -m clamacitor`).

## 1.0.0

Code and data accompanying the article (doi:10.5281/zenodo.18854151).
