# IvermectinCattle

This repository contains the average and population physiologically based pharmacokinetic (PBPK) model code files used in our ivermectin PBPK study in cattle. The model was developed to estimate withdrawal intervals (WDIs) following labeled (0.2 mg/kg) and extra-label (e.g., 0.3 mg/kg) drug use.

The manuscript title of this study is “Development of a physiologically based pharmacokinetic model for ivermectin in cattle.” This manuscript is under review.

## 📁 Repository Contents

- 📂 **Model Calibration and Evaluation** — R code for the general PBPK model used for model calibration and evaluation.
- 📂 **Model Calibration and Evaluation_ModFit** — R code for sensitivity analysis and parameter optimization.
- 📂 **Population PBPK** — R code for the population PBPK model used to estimate withdrawal intervals. Please refer to Supplementary Table S1 for the parameter distributions used in the population simulations. WDI predictions are based on three regulatory criteria: historical FDA tolerance (1995) for liver, current FDA tolerances (2020) for liver and muscle, and Codex maximum residue limits (MRLs) for liver, kidney, muscle, and fat.
- 📂 **Pharmacokinetic dataset** — Pharmacokinetic datasets used for model calibration and evaluation.
