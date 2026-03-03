[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18854151.svg)](https://doi.org/10.5281/zenodo.18854151)

# Flux Clamacitor — *Corbicula fluminea* CaCO₃ Mass Balance Model

**Companion code and interactive tool for:**

> Hartland, A., Melchior, M., Hamilton, D., Lehto, N.J., Mullarney, J., Sandwell, D., Robb, L., Jabbari, A., Lang, J., Clague, J., Özkundakci, D., and Hofstra, D. (*in review*). **Massive calcium drawdown by invasive bivalves undermines arsenic removal in water supply from a major river system.** *Environmental Science & Technology.*

---

## Overview

This repository contains the Python analysis notebooks, input datasets, and an interactive HTML tool used to derive CaCO₃ biomineralization fluxes and population-scale density estimates for *Corbicula fluminea* (Asian clam) in Lake Karāpiro, Waikato River, New Zealand. The study links a novel *C. fluminea* invasion to large-scale calcium depletion and downstream arsenic treatment failures at a major drinking water supply.

The **Flux Clamacitor** (`FluxClamacitor.html`) is a self-contained browser-based GUI allowing users to run the removal flux model interactively, adjust parameters, and compare model outputs against field survey densities — with no Python installation required.

---

## Repository Structure

```
├── FluxClamacitor.html                               # Interactive GUI (open in any browser)
├── LaunchFluxClamacitor.bat                          # Windows launcher for the GUI
├── README.md                                         # This file
│
├── notebooks/
│   └── C_fluminea_mass_balance_study.ipynb           # Full analysis notebook (5 cells)
│
└── data/                                             # Input data (see Data section below)
    ├── Density_and_shell_mass_for_AH.xlsx
    ├── Full_data_for_Corbicula_Calcs.xlsx
    ├── Full_data_for_Corbicula_Calcs_Landscape_Updated_16_1_26.xlsx
    ├── processed_corbicula_biovolumes.xlsx
    ├── Monthly_Water_Quality_Narrows_Summary.csv
    └── Revised_Corbicula_Mass_Balance.csv
```

> **Note:** All data files are provided as supplementary material accompanying the published paper. Download them from the journal supplement and place them in the `data/` directory before running the notebooks.

---

## The Model

The mass balance framework attributes net dissolved Ca and alkalinity depletions to biogenic CaCO₃ precipitation during *C. fluminea* shell growth, scaled through four sequential equations:

| Step | Equation | Output |
|------|----------|--------|
| **Eq. 1** | Longitudinal Ca flux: *J*ₗₒₙ𝓰 = ([Ca]ᵤₚ·Qᵤₚ − [Ca]ᵈₒ𝓌ₙ·Qᵈₒ𝓌ₙ) × 86400 / M_Ca | mol CaCO₃ d⁻¹ |
| **Eq. 2** | Longitudinal Alk flux (×½ stoichiometric factor) | mol CaCO₃ d⁻¹ |
| **Eq. 3** | Vertical depletion flux (CSTR): *J*ᵥₑᵣₜ = ΔC × (V/τ) / M | mol CaCO₃ d⁻¹ |
| **Eq. 4** | Mass fixed: *m*_CaCO₃ = *J* × M_CaCO₃ × 10⁻⁶ | t CaCO₃ d⁻¹ |
| **Eq. 5** | Biovolume rate: *V*_bio = *m*_CaCO₃ × shell ratio / bulk density | m³ d⁻¹ |
| **Eq. 6** | Population: *N* = *V*_bio / (d*V*/d*t*) | individuals |
| **Eq. 7** | Density: *D* = *N* / colonizable area | ind m⁻² |

Negative fluxes (dissolution, die-off events, or measurement artefacts) are excluded from population scaling. Abiotic calcite precipitation is ruled out by systematic calcite undersaturation (SI_cc < 0) throughout the study period.

---

## Notebook Contents

The analysis notebook (`C_fluminea_mass_balance_study.ipynb`) is structured in five cells:

1. **Cell 1 — Biovolume factor derivation:** Calculates per-individual biovolume (ellipsoid formula) and dry density from shell morphometrics by size class (small, medium, large). Exports `processed_corbicula_biovolumes.xlsx`.

2. **Cell 2 — Biomineralization and population metrics:** Implements Eqs. 1–7 across all four flux methods (Ca spatial, Alk spatial, Ca residence, Alk residence). Exports `Revised_Corbicula_Mass_Balance.csv`.

3. **Cell 3 — Residence time correlations:** Tests Pearson correlations between hydraulic residence time (τ) and removal fluxes / concentration gradients for growth periods (positive flux only). Produces Fig. S3.

4. **Cell 4 — Formation rate visualisation:** Six-panel publication figure (Fig. 4) showing Ca and alkalinity time series against historical baselines, and CaCO₃ removal fluxes by method with residence time overlay.

5. **Cell 5 — Survey density comparison:** Box plots comparing modelled density distributions against 2025 field survey data from five sites (Waipuke, Horahora, Moana Roa, Bobs Landing, Little Waipa). Produces Fig. 5.

---

## Using the Flux Clamacitor GUI

### Quick Start (Windows)
1. Download `FluxClamacitor.html` and `LaunchFluxClamacitor.bat` to the same folder.
2. Double-click `LaunchFluxClamacitor.bat`.
3. The tool opens in your default browser — no installation required.

### Manual Launch
Open `FluxClamacitor.html` directly in any modern browser (Chrome, Firefox, Edge).

### Editable Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| Lake volume | Dynamic lake volume (*V*) | 85 × 10⁶ m³ |
| Colonizable benthic area | Area < 1 m depth | 7.41 × 10⁶ m² |
| Individual growth rate (d*V*/d*t*) | Steady-state median from biovolume data | 0.0055 cm³ d⁻¹ ind⁻¹ |
| Total shell ratio | g total shell / g CaCO₃ shell | 1.13 |
| Shell bulk density | Packing density | 1.15 t m⁻³ |
| Flux concentrations | Ca (mg/L), Alk (mg/L as CaCO₃), Q (m³/s) | From ECP monitoring |
| Survey site densities | Field-observed densities per site | May 2025 survey |

### Fixed Constants (non-editable)

| Constant | Value | Reason fixed |
|----------|-------|--------------|
| Molar mass Ca | 40.08 g mol⁻¹ | Periodic table |
| Molar mass CaCO₃ | 100.09 g mol⁻¹ | Stoichiometric |
| Alk equivalent weight | 50.043 g eq⁻¹ | CaCO₃ equivalence |

---

## Requirements

### Python Notebooks
```
python >= 3.9
pandas
numpy
matplotlib
scipy
openpyxl
```

Install dependencies:
```bash
pip install pandas numpy matplotlib scipy openpyxl
```

### Flux Clamacitor GUI
No installation required. Runs entirely in-browser using vanilla JavaScript and [Chart.js 4.4.1](https://www.chartjs.org/) (loaded from CDN).

---

## Data

All input data files are provided in the supplementary material of the companion paper. The datasets include:

- **Density and shell mass** — *C. fluminea* density surveys and individual shell morphometrics from Lake Karāpiro sites (2024–2025).
- **Full water quality data** — Fortnightly Ca, alkalinity, and flow records from ECP monitoring stations ECP-S1 (Aniwaniwa Reserve), ECP-S2 (Lake Karāpiro profiler), and ECP-S3 (Cambridge Golf Course), October 2024–May 2025.
- **Monthly Narrows summary** — Historical water quality statistics (Ca, alkalinity) from Waikato Regional Council monitoring (1995–2021) used as pre-invasion baseline.
- **Depth–residence time data** — Processed depth profile data with computed hydraulic residence times.
- **Revised mass balance CSV** — Full model output from Cell 2, used as input for Cells 3–5.

---

### Citation

If you use this code or the Flux Clamacitor tool, please cite:

Hartland, A., Melchior, M., Hamilton, D., Lehto, N.J., Mullarney, J., Sandwell, D.,
Robb, L., Jabbari, A., Lang, J., Clague, J., Özkundakci, D., and Hofstra, D. (in review).
Massive calcium drawdown by invasive bivalves undermines arsenic removal in water supply
from a major river system. *Environmental Science & Technology*.

**Code:**

Hartland, A.; Melchior, M.; Hamilton, D.; Lehto, N.J.; Mullarney, J.; Sandwell, D.;
Robb, L.; Jabbari, A.; Lang, J.; Clague, J.; Özkundakci, D.; Hofstra, D. (2025).
*Flux Clamacitor v1.0.0 — C. fluminea CaCO₃ mass balance and population scaling
model (Waikato River, NZ)*. Zenodo. https://doi.org/10.5281/zenodo.18854151

---

## Funding and Permissions

This work was funded by the New Zealand Ministry of Business, Innovation and Employment
(MBIE) Endeavour programme grants to Lincoln Agritech Ltd (LVLX2302) and Earth Sciences
New Zealand (NIW2475), with in-kind and cash contributions from Waikato Regional Council.
All sampling and analytical procedures complied with MPI biosecurity protocols for
containing the *C. fluminea* invasion (MPI permit P0134).

We acknowledge Ngāti Kōroki Kahukura, kaitiaki of Maunga Maungatautari and Lake Karāpiro,
and Waikato-Tainui iwi.

---

## AI Assistance

The `FluxClamacitor.html` interactive GUI was developed with the assistance of Claude Sonnet 4.6 (Anthropic), using the underlying Python codebase and mass balance framework as the source logic. All scientific content, equations, parameter values, and analytical methodology were authored by the research team and verified against the companion paper. AI assistance was used solely for GUI development and did not contribute to data analysis, interpretation, or manuscript preparation.

This disclosure is provided in accordance with the AI use policy of the *American Chemical Society* and *Environmental Science & Technology*.

---

## Contact

**Adam Hartland** — adam.hartland@lincolnagritech.co.nz  
Lincoln Agritech Ltd, Ruakura, Hamilton, New Zealand
