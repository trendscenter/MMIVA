# Multimodal Independent Vector Analysis

This repository contains code for the paper: [A method for multimodal IVA fusion within a MISA unified model reveals markers of age, sex, cognition, and schizophrenia in large neuroimaging studies](https://www.biorxiv.org/content/10.1101/2021.12.13.472507)

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
  [runMMIVA_allHCSZ.m](/HCSZ/MMIVA/runMMIVA_allHCSZ.m),
  [runMMIVA_allHCSZ_preregSite.m](/HCSZ/MMIVA/runMMIVA_allHCSZ_preregSite.m)

UKB initialization:
  [runMMIVA_allHCSZ_w_init.m](/HCSZ/MMIVA/runMMIVA_allHCSZ_w_init.m),
  [runMMIVA_allHCSZ_w_initUKB2907_preregSite.m](/HCSZ/MMIVA/runMMIVA_allHCSZ_w_initUKB2907_preregSite.m)

### MANCOVA
GICA initialization: 
  [runStats_HCSZ_preregSite_C030_gicainit.m](/HCSZ/MANCOVA/runStats_HCSZ_preregSite_C030_gicainit.m)

UKB initialization: 
  [runStats_HCSZ_preregSite_C030_ukbinit.m](/HCSZ/MANCOVA/runStats_HCSZ_preregSite_C030_ukbinit.m)

## Reference
If you find this repository useful, please cite the following paper:
```
@article{silva2021direct,
  title={Direct linkage detection with multimodal IVA fusion reveals markers of age, sex, cognition, and schizophrenia in large neuroimaging studies},
  author={Silva, Rogers F and Damaraju, Eswar and Li, Xinhui and Kochunov, Peter and Belger, Aysenil and Ford, Judith M and Mathalon, Daniel H and Mueller, Bryon A and Potkin, Steven G and Preda, Adrian and others},
  journal={bioRxiv},
  pages={2021--12},
  year={2021},
  publisher={Cold Spring Harbor Laboratory}
}
```