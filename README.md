# Multimodal Independent Vector Analysis

This repository contains code for the paper: [A Method for Multimodal IVA Fusion Within a MISA Unified Model Reveals Markers of Age, Sex, Cognition, and Schizophrenia in Large Neuroimaging Studies](https://onlinelibrary.wiley.com/doi/10.1002/hbm.70037)

Note that you would need to clone [MISA](https://github.com/rsilva8/MISA.git) repository under MMIVA folder.

```
git clone https://github.com/trendscenter/MMIVA.git
cd MMIVA
git clone https://github.com/rsilva8/MISA.git
```

## UK Biobank Data
### MMIVA
  [runMMIVA_C030_regSite_N2907.m](UKB/MMIVA/runMMIVA_C030_regSite_N2907.m)

### MANCOVA
  [runStats_UKB_MIVA2step_wICApre_preregSite_C030_N2907.m](/UKB/MANCOVA/runStats_UKB_MIVA2step_wICApre_preregSite_C030_N2907.m)

## Patient Data
### MMIVA
GICA initialization: 
  [runMMIVA_allHCSZ.m](/HCSZ/MMIVA/GICA/runMMIVA_allHCSZ.m),
  [runMMIVA_allHCSZ_preregSite.m](/HCSZ/MMIVA/GICA/runMMIVA_allHCSZ_preregSite.m)

UKB initialization:
  [runMMIVA_allHCSZ_w_init.m](/HCSZ/MMIVA/UKB/runMMIVA_allHCSZ_w_init.m),
  [runMMIVA_allHCSZ_w_initUKB2907_preregSite.m](/HCSZ/MMIVA/UKB/runMMIVA_allHCSZ_w_initUKB2907_preregSite.m)

### MANCOVA
GICA initialization: 
  [runStats_HCSZ_preregSite_C030_gicainit.m](/HCSZ/MANCOVA/runStats_HCSZ_preregSite_C030_gicainit.m)

UKB initialization: 
  [runStats_HCSZ_preregSite_C030_ukbinit.m](/HCSZ/MANCOVA/runStats_HCSZ_preregSite_C030_ukbinit.m)

## Reference
If you find this repository useful, please cite the following paper:
```
@article{https://doi.org/10.1002/hbm.70037,
author = {Silva, Rogers F. and Damaraju, Eswar and Li, Xinhui and Kochunov, Peter and Ford, Judith M. and Mathalon, Daniel H. and Turner, Jessica A. and van Erp, Theo G. M. and Adali, Tulay and Calhoun, Vince D.},
title = {A Method for Multimodal IVA Fusion Within a MISA Unified Model Reveals Markers of Age, Sex, Cognition, and Schizophrenia in Large Neuroimaging Studies},
journal = {Human Brain Mapping},
volume = {45},
number = {17},
pages = {e70037},
doi = {https://doi.org/10.1002/hbm.70037},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/hbm.70037},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/hbm.70037},
year = {2024}
}
```
