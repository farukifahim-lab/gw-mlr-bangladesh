# Interpretable MLR equations for groundwater level in Bangladesh

Code and processed data for: *Groundwater Level Forecasting in Bangladesh Using an Interpretable Multiple Linear Regression Equation* (WRR submission 2026WR046512).

## Contents
- `GW_MLR_reproducible.R` – full workflow (cleaning, 70/30 split, full MLR, VIF screening, standardized-beta ranking, Top-N models, equations, extreme-value metrics).
- `data/Final-ALL-MSD.csv` – 15,000 station-month records, GWL plus 57 candidate predictors (add after BWDB permission is confirmed; see below).

## Run
`Rscript GW_MLR_reproducible.R` from this folder. Required packages: caTools, hydroGOF, car, lmtest.

## Data notes
- GWL: Bangladesh Water Development Board (BWDB); quality-controlled and gap-filled by the authors.
- Predictors: ERA5, MERRA-2 (NASA POWER), GLDAS, GRACE and static geographic layers (see manuscript Table 1).
- `GWETPROF` and `GWETROOT` contain the missing-value code -999 in 20 records of one well; the script removes them.
- Station statistics (`Mean-GWL`, `Var-GWL`, `Range-GWL`, `MX_GWL`, `Mean-2SD`) are computed from each well's record.

## License
Code: MIT. Data: CC BY 4.0, subject to BWDB redistribution terms.
