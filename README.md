# Purple sea urchin energetics under ocean warming, acidification, and thermal variability

**R analysis code and supporting data for the 2021 Hakai–Marna laboratory mesocosm experiment**

[GitHub repository](https://github.com/nspindel/2021_Hakai_Marna_Urchin_Energetics) · [Citation metadata](CITATION.cff) · [Software license](LICENSE)

## Overview

This repository contains the analytical R code, experimental datasets, environmental context files, and figure-generation scripts associated with the manuscript:

> **Contemporary ocean acidification and marine heatwaves shape individual energetics and rates of herbivory in a dominant ecosystem engineer.**  
> Nathan B. Spindel and coauthors. Manuscript prepared for submission to *Limnology and Oceanography*.

The study quantifies how seawater temperature, elevated partial pressure of carbon dioxide (pCO₂), and temporally variable thermal exposure affect the physiology and energetics of the purple sea urchin (*Strongylocentrotus purpuratus*). The central responses include resting metabolic rate, per-capita food consumption, assimilation efficiency, gonadal biomass production, food conversion efficiency, and skeletal (demipyramid) growth. Experiments were conducted in autumn 2021 at the Hakai Institute's Quadra Island Ecological Observatory, British Columbia, Canada.

**Repository scope:** This is a research analysis repository, not a general-purpose R package. The supplied source snapshot contains source datasets and scripts but is **not yet a fully self-contained, one-command reproduction archive** (see [Reproducibility status](#reproducibility-status-and-known-limitations)). The accompanying curated biological dataset is intended for Dryad; a versioned snapshot of this code will be archived by Zenodo on GitHub release. Permanent DOI links will be added when assigned.

## Study design and observational units

- **Collection:** Adult urchins were collected by SCUBA near Ucluelet, British Columbia (approximately 48.94° N, 125.56° W), in September 2021.
- **Mesocosms:** 32 independently controlled 340-L flow-through seawater mesocosms, stocked initially with approximately 15 PIT-tagged urchins each.
- **Constant-temperature treatments:** 10, 13, 16, 17, 18, and 20 °C crossed with nominal pCO₂ levels of 600 and 1,200 µatm.
- **Variable-temperature treatments:** La Niña-like (18–14 °C; mean approximately 16 °C) and El Niño-like (21–18 °C; mean approximately 20 °C) regimes, both at nominal 600 µatm pCO₂.
- **Food supply:** Formulated macroalgal pellets supplied *ad libitum*.
- **Experimental unit:** The **mesocosm** is the independently assigned treatment replicate. Tagged individuals and repeated sampling events are nested within mesocosms; they must not be analyzed as independent treatment replicates.

Nominal treatment values are not interchangeable with measured daily seawater conditions. See the supplementary manuscript methods for details of thermoregulation and carbonate chemistry monitoring.

## Repository organization

The following paths were verified in the source ZIP provided for preparation of this documentation. The table describes intended use, not a claim that every script can currently execute from beginning to end without intervention.

| File or directory | Contents / purpose |
| --- | --- |
| `2021_Hakai_Marna_Urchin_Energetics.Rproj` | RStudio project; use the project root as the working directory. |
| `Utility_Functions.R` | Shared plotting, model-prior, and data-loading helpers. |
| `Biometrics_Merge.R` | Initial/final individual matching and derived body-size change metrics. |
| `Metabolic_Rate_Estimation.R` | Respirometry processing, Bayesian metabolic-rate models, contrasts, and figures. |
| `Consumption_Estimation.R` | Group food-consumption estimates, thermal/pCO₂ models, contrasts, and figures. |
| `Assimilation_Efficiency_Estimation.R` | Individual feeding/egestion calculations and assimilation-efficiency models. |
| `Gonad_Estimation.R` | Wild-urchin initial-gonad allometry, derived gonad production, and models. |
| `Food_Conversion_Efficiency.R` | Estimated intake, gonad gain, food-conversion efficiency, and models. |
| `Skeletal_Growth_Estimation.R` | Tetracycline-marked demipyramid growth and related models. |
| `Sensor_Plots.R` | Mesocosm and regional sensor time-series graphics. |
| `Mapping.R` | Environmental context, mapping, and additional figure production. |
| `Data/` | Biological observations, derived analytical tables, field-lab workbooks, and supporting environmental data (inventory below). |
| `Figures/` | Existing figure exports and supplementary graphics (PDF/PNG/JPG). |
| `LICENSE` | Original code's GNU General Public License, version 3. |
| `CITATION.cff` | Machine-readable citation metadata for GitHub and Zenodo. |

### Principal biological files in `Data/`

| Input file | Unit / content |
| --- | --- |
| `Biometrics_Initial.xlsx` | Initial per-individual PIT-tagged measurements: date, tank, treatment, test dimensions, and wet mass. |
| `Biometrics_Final.xlsx` | Final individual measurements and dissection records, including gonad and other tissue masses. |
| `Respirometry.csv` | 843 **processed** oxygen-consumption observations; contains treatment, tank, PIT tag, body mass/ash-free dry mass, and oxygen-consumption rate, **not** full optical-sensor time series. |
| `Aggregated_Feeding.csv` | 256 observations from repeated group feeding assays; offered/remaining food dry masses, group sizes, dates, and elapsed time. |
| `Assimilation_Efficiency.xlsx` | Individual ingestion/egestion assay measurements (80 records); primary workbook read by the assimilation script. |
| `Assimilation_Efficiency.csv` | Alternative CSV export of the preceding workbook (80 records); two columns have empty headings. |
| `Initial_Wild_Morphometrics.csv` | Allometry calibration sample from 35 wild urchins: test volume and dry gonad mass. |
| `Derived_Gonad_Production.csv` | 470 derived individual-level records combining dissection values and estimated initial gonad masses. **Validation caveat below.** |
| `FCE_Data.csv` | 466 individual-level records combining wet gonad masses, modeled initial gonad mass, and estimated food intake. |
| `Demipyramid_Growth.csv` | 464 individual records of tetracycline mark position and final skeletal extension. |
| `Skeletal_Growth.csv` | 464 records joining external biometrics and demipyramid growth; substantially overlaps the preceding table. |

**Supporting environmental inputs.** `Data/Sensor_Data/Marna_Lab_Mesocosm_Temps.xlsx` contains experimental temperature measurements. In `Data/Sensor_Data/Urchins Fall 2021 - Adult & Larvae PCO2 - Data.xlsx`, the sheet **`Adult Urchin 2021 PCO2 Data Fin`** is associated with this adult-urchin experiment; other sheets relate to a separate larval experiment and should not be pooled. `Data/Bamfield/`, `Data/Newport_OR/`, and the `Data/Sensor_Data/QuadraBoL/` and `0208638/` directories contain contextual records from external monitoring sources. Multiple QuadraBoL files are duplicated across directories. They are **not** experimental replicates, and original provider terms must be checked before redistribution.

## Key variables and derived responses

| Name or common fields | Meaning / units |
| --- | --- |
| `tank_id`, `tank` | Mesocosm/treatment replicate identifier. |
| `pit_tag_id`, `pit_tag`, `pit_id` | Individual passive integrated transponder identifier; names vary among source files. |
| `temp_deg_c`, `temp` | Assigned thermal treatment. `18-14` and `21-18` are **character labels for trajectories**, not numeric values or arithmetic expressions. |
| `pco2_uatm`, `pco2_h_l`, `pco2` | Assigned nominal pCO₂ treatment (µatm); in some files `l`/`h` denotes lower/higher pCO₂. Recode explicitly before merging. |
| `trial_id`, `trial_day` | Feeding trial occasion or experimental day (not necessarily calendar date). |
| `diam_mm`, `ht_mm`, `initial_*`, `final_*` | Test diameter/height (mm) and time-specific measurements. |
| `wet_mass_g`, `afdm_g`, `*_wet_*_g`, `*_dry_*_g` | Body/tissue mass (g); AFDM is ash-free dry mass. |
| `absolute_rate_vo2_mg_o2_per_hour` | Oxygen-consumption rate in mg O₂ h⁻¹ before model-based body-size standardization. |
| `initial_food_dry_mass_g`, `leftover_food_dry_mass_g` | Group food offered and uneaten (g dry mass). |
| `single_gonad_dry_mass_g` | Measured dry mass of one gonad lobe (g); scripts multiply by five to estimate total dry gonad mass. |
| `init_nad_est`, `estimated_initial_dry_gonad` | **Estimated** baseline gonad quantities; these are not direct measurements and differ in units/processing across analysis tables. |
| `daily_consumption_g_dry_est` | Model-estimated daily dry food consumption (g individual⁻¹ day⁻¹). |
| `demi_tag_um`, `demi_final_um`, `total_days` | Tetracycline reference/final positions (µm) and elapsed duration (days). |

**Core calculations implemented in code:**

- **Test volume** (mL), approximating an oblate spheroid: `(4/3) × π × (diameter_mm/2)^2 × (height_mm/2) / 1000`.
- **Per-capita consumption** (g dry food individual⁻¹ day⁻¹): `(initial dry food − remaining dry food) / number of urchins / (feeding hours / 24)`.
- **Assimilation efficiency** (unitless): `(ingested dry food − egested dry feces) / ingested dry food`.
- **Estimated gonad production:** final gonad mass minus the modeled initial gonad mass. Experimental initial gonad biomass was **not** measured non-destructively; the baseline model was calibrated against 35 separate wild urchins.
- **Food conversion efficiency (FCE):** estimated gonad gain divided by estimated cumulative dry food consumption. The `Food_Conversion_Efficiency.R` implementation uses **wet gonad-mass gain** in the numerator and **dry food consumption** in the denominator; do not interpret this as a dry-to-dry tissue conversion ratio.
- **Demipyramid extension rate** (µm day⁻¹): `(demi_final_um − demi_tag_um) / total_days`. This measures jaw extension, not calcium carbonate accretion or dissolution directly.

The analysis uses Bayesian regression in `brms`/Stan, including body-size covariates, tank-level random effects, and Gaussian-process temperature terms for selected responses. Candidate models and posterior contrasts are described in the manuscript and Supporting Information.

## Software requirements and getting started

- **R** with **RStudio** recommended; use the included `.Rproj` to establish the project working directory.
- **Stan** and the `cmdstanr` backend for the Bayesian fits. Install and configure CmdStan separately according to the current `cmdstanr` instructions.
- R dependencies used across the scripts include `tidyverse`, `janitor`, `readxl`, `brms`, `rstan`, `cmdstanr`, `posterior`, `tidybayes`, `emmeans`, `glmmTMB`, `DHARMa`, `ggrepel`, `scales`, `Hmisc`, `boot`, `marmap`, `terra`, `raster`, `sf`, `rnaturalearth`, `rnaturalearthdata`, `lubridate`, `here`, `forcats`, `reshape2`, `viridis`, `zoo`, and `plyr`. Base/recommended R packages such as `parallel` are also used.

To inspect or resume the analyses:

1. Clone or download the repository and open `2021_Hakai_Marna_Urchin_Energetics.Rproj` in RStudio.
2. Install the dependencies relevant to the response variable you are analyzing. Install CmdStan if fitting Bayesian models using `cmdstanr`.
3. Verify that the required `Data/` file for the selected analysis exists and that its column names and treatment encodings match the script.
4. Read the script's data-preparation and model-fitting sections **before** running it. The scripts are intended as research workflows, not parameterized functions; some sections assume previously saved model fits.
5. Rerun or supply the missing model fits and revise output paths as necessary before attempting to recreate downstream figures.

There is no `renv.lock`, package version lockfile, automated test suite, or single entry-point runner in the supplied snapshot. Dependencies and execution order are consequently not pinned for exact environment recreation.

## Reproducibility status and known limitations

The following items were identified directly in the supplied repository snapshot and should be resolved, documented, or handled explicitly before describing the Zenodo release as fully reproducible:

1. **Missing model objects:** Multiple scripts call `readRDS('./Model_Output/*.rds')`, but the supplied ZIP contains **no** `Model_Output/` files. Fitting code is present in the scripts, yet some plotting and contrast sections require fitted objects to be regenerated or supplied.
2. **Missing derived biometrics output:** `Biometrics_Merge.R` writes/reads `Data/Data_Final/Biometrics_Growth_PercentChange.csv`, but the supplied archive contains neither that directory nor the derived CSV. Create the output directory and inspect unmatched and duplicate PIT-tag records before running this step.
3. **QC of gonad derivatives:** `Derived_Gonad_Production.csv` contains historical model-output columns, including anomalously large values in `init_nad_est`. They should **not** be treated as verified gonad masses without reconciling the derived table to the latest corrected initial/final biometrics and allometric model.
4. **Analysis-specific filtering:** `Consumption_Estimation.R` replaces negative calculated consumption rates with `0.0001`; `Food_Conversion_Efficiency.R` uses `fce > 0` for model fitting; `Skeletal_Growth_Estimation.R` restricts demipyramid rates to 0.1–5 µm day⁻¹; `Assimilation_Efficiency_Estimation.R` adjusts values at the boundaries of (0, 1) for beta regression. These are modeling decisions and should **not** silently replace original observations in a curated data deposit.
5. **Input provenance and format differences:** Some source Excel-to-CSV exports retain unnamed fields or numeric Excel dates. External environmental archives are substantially larger than the core experiment files and may have independent redistribution and acknowledgment requirements.
6. **Interpretive scope:** Feeding was measured with formulated food under food-replete conditions. The measurements do not directly quantify removal of attached kelp under field conditions, nor do resting oxygen-consumption measurements establish the complete energy budget.

A curated Dryad dataset should retain observed records and auditable transformations, with a data dictionary and documented exclusions. Uncertainty in the wild-urchin allometric model used to estimate initial gonad mass is not fully propagated through all subsequent response models, as noted in the manuscript.

## Study team and contributions

The **scientific author list below follows the current manuscript title page** (including affiliation numbering). All research authors contributed to manuscript refinement; **the primary R analysis scripts identify Nathan B. Spindel as their author**. Including the scientific team in citation metadata acknowledges the broader empirical contribution and does not imply that each person wrote analytical source code.

| Author | Affiliation(s) from title page |
| --- | --- |
| Nathan B. Spindel (corresponding) | Florida State University; Moss Landing Marine Laboratories, San Jose State University; California State University Monterey Bay |
| Maya Munstermann | Florida State University; University of California, Berkeley |
| Sam Karelitz | Florida State University; University of California, Berkeley |
| Rachele Ferraro | Florida State University |
| Ondine Pontier | Hakai Institute |
| Brenna Collicutt | Hakai Institute |
| Kate Rolheiser | Hakai Institute |
| Iria Gimenez | Hakai Institute |
| Evan Cronmiller | Hakai Institute |
| Megan Foss | Hakai Institute |
| Natalie Mahara | Hakai Institute |
| Daniel S. Swezey | Bodega Marine Laboratory, University of California Davis |
| Laura Rogers-Bennett | Bodega Marine Laboratory, University of California Davis |
| Daniel K. Okamoto | Florida State University; University of California, Berkeley |

**Roles as documented on the title page:** N. B. Spindel conceived the study, designed methods, led the mesocosm work, curated and analyzed data, produced visualizations, and drafted the manuscript. M. Munstermann, S. Karelitz, R. Ferraro, O. Pontier, B. Collicutt, K. Rolheiser, I. Gimenez, E. Cronmiller, M. Foss, and N. Mahara carried out field collections and experimental work. D. S. Swezey, L. Rogers-Bennett, and D. K. Okamoto contributed to study conceptualization, experimental design, and interpretation. L. Rogers-Bennett and D. K. Okamoto also secured funding, supplied resources, and oversaw project administration and supervision.

**Correspondence:** Nathan B. Spindel, [nspindel@csumb.edu](mailto:nspindel@csumb.edu).

**Support:** NSF OCE 2023649 (D. K. Okamoto), NSF OCE 2023664 (L. Rogers-Bennett), Tula Foundation, and a PADI Foundation Research Grant (N. B. Spindel). The Hakai Institute provided research facilities and technical support.

## How to cite and archive

- **Code:** Cite the archived **version-specific Zenodo DOI** for the GitHub release used. `CITATION.cff` contains manuscript-aligned citation metadata but intentionally has **no DOI or software version** until the first release is assigned and finalized.
- **Data:** Cite the curated Dryad dataset and its **Dryad DOI** once deposited. The GitHub repository's historical data copies are not a substitute for the curated dataset of record.
- **Article:** Cite the associated *Limnology and Oceanography* article after publication; do not cite a nonexistent journal DOI while the paper remains under submission.

To archive a release, connect the public GitHub repository to Zenodo, enable archiving, then create a tagged GitHub release (e.g., `v1.0.0`). Once Zenodo has minted the DOI, update the GitHub README with links to the archived release and Dryad dataset, and add the stable Zenodo concept DOI to future `CITATION.cff` versions as desired. Do not put invented or provisional DOI strings into metadata.

## Licensing and data provenance

**Software:** The repository already includes the **GNU General Public License, version 3**. This documentation preserves that existing license for original analysis source code; see [`LICENSE`](LICENSE). License identifier in citation metadata: `GPL-3.0-only`.

**Data:** Dryad-published primary research data will be subject to Dryad's CC0 public-domain dedication. That designation **does not automatically cover third-party source files stored in this GitHub repository**. Hakai, OOI, NCEI, and other external sensor datasets, base maps, photographs, and externally sourced materials retain their own provenance and any applicable terms. Verify redistribution rights and attribution before packaging them into a Zenodo release or public repository. Research involving nonhuman marine invertebrates has no human-subject records in the datasets described here.

The source code is distributed without warranty under the terms of the GPLv3.
