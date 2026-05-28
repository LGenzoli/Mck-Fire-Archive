# Script and Data Overview

This README outlines the scripts used to complete the analysis for:

> **"Debris flows suppressed riverine productivity and respiration following high-severity wildfire on the Klamath River, California"**

## Notes

- With the exception of the 6 days immediately before, during, and after the debris flow, high-frequency data is excluded from public data repositories by request of the Karuk Tribe. This high-frequency data (used to model metabolism and calculate daily water quality stats including turbidity and dissolved oxygen) can be requested directly from the Karuk Tribe.
- Raw data used in the `Depth_model` is from Bradley (2021) and can be requested from the author.
- Scripts were checked and last run on **24 Feb 2024**, with the exception of the metabolism models, which take too long to run to justify a casual re-check.
- For questions about the data or scripts, contact **Laurel Genzoli**: laurel.genzoli@gmail.com

---

## Scripts

### `SV-MetabCalcs.R`
Runs streamMetabolizer at Seiad Valley from 2018–2023. The first run uses the Bayesian framework where K600 is estimated as a free parameter; the second run uses fixed K600 values estimated from discharge. This script works iteratively with `Final_Metab_ModelOut.qmd`.

**Input data:**
- `SV_clean-23072025.csv` *(raw data available from Karuk Tribe on request)*
- `SV_clean-MLE23072025.csv` *(raw data available from Karuk Tribe on request)*

**Output:**
- `Metab_Output/Binned/SV_GPP_ER_BIN.csv`
- `Metab_Output/Binned/SV_RDA_BIN.rda`
- `Metab_Output/Binned/SV_fits_BIN.csv`
- `Metab_Output/MLE/SV_GPP_ER_BIN.csv`
- `Metab_Output/MLE/SV_RDA_BIN.rda`
- `Metab_Output/MLE/SV_fits_BIN.csv`

---

### `Depth_model.qmd`
Calculates the slope and intercept needed to estimate reach depth based on daily discharge at Seiad, using depth estimates from the Bradley Model.

**Input data:**
- `LongProfiles/DenseMesh` *(data available from Bradley 2021)*
- `LongProfiles/LongProfiles` *(data available from Bradley 2021)*

**Output:** Model coefficients to estimate reach depth based on daily flow.

---

### `Final_Metab_ModelOut.qmd`
Prepares the final metabolism data frame. Examines output from the binned model, builds the K600–discharge relationships for the MLE run, then examines the MLE output for model fit and realistic values. Builds the final metabolism data frame and the high-frequency data frame. Includes reach length calculations.

**Input data:**
- `SV_clean-23072025.csv` *(raw data available from Karuk Tribe on request)*
- `SV_turbidity_update_2025121700173883.csv` *(raw data available from Karuk Tribe on request)*
- `SV_2018-2024_AllParams.csv` *(raw data available from Karuk Tribe on request)*
- `Metab_Output/Binned/SV_GPP_ER_BIN.csv`
- `Metab_Output/Binned/SV_fits_BIN.csv`
- `Metab_Output/Binned/SV_RDA_BIN.rda`
- `Metab_Output/MLE/SV_fits_BIN.csv`
- `Metab_Output/MLE/SV_GPP_ER_BIN.csv`
- `Metab_Output/MLE/SV_RDA_BIN.rda`
- `SV_light_2018-2023-20250723.csv`

**Output:**
- `met_mle_final.csv`
- `sub-daily-FIRE.csv`

---

### `Final_models.qmd`
Runs the time series model and the breakpoint model. **This script must be run before `Final_Figures.qmd`**, as Figure 2 depends on model output generated in the same session.

**Input data:**
- `met_mle_final.csv`

**Output:** Model coefficients for final tables and figures.

---

### `Final_Figures.qmd`
Produces all final figures, including most supplementary figures. Figure 1 and Figure S6 require post-processing layout.

**Input data:**
- `met_mle_final.csv`
- `sub-daily-FIRE.csv`

**Output:** Figures
