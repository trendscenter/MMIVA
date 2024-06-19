import os
import h5py
import numpy as np
import scipy.io as sio

ukb_smri_data_path = '/data/users1/xinhui/MISA/MMIVA/UKB/postMMIVAresults/mancovanouts/UKB_MMIVA_C30_preregSite_SMRI_MancovanOuts_wX_FINAL.mat'
ukb_fmri_data_path = '/data/users1/xinhui/MISA/MMIVA/UKB/postMMIVAresults/mancovanouts/UKB_MMIVA_C30_preregSite_ALFF_MancovanOuts_wX_FINAL.mat'
ukb_dmri_data_path = '/data/users1/xinhui/MISA/MMIVA/UKB/postMMIVAresults/mancovanouts/UKB_MMIVA_C30_preregSite_DMRI_MancovanOuts_wX_FINAL.mat'

# load actual age
ukb_smri_data = sio.loadmat(ukb_smri_data_path)['MODELUKB0s_ful']
ukb_smri_data_array = ukb_smri_data[0][0][0]
ukb_smri_data_key = np.squeeze(ukb_smri_data[0][0][3])
age_idx = np.where(ukb_smri_data_key==['age_when_attended_assessment_centre_f21003_2_0'])[0][0]
sex_idx = np.where(ukb_smri_data_key==['sex_f31_0_0'])[0][0]
Y = ukb_smri_data_array[:, age_idx]

ukb_fmri_data = sio.loadmat(ukb_fmri_data_path)['MODELUKB1s_ful']
ukb_fmri_data_array = ukb_fmri_data[0][0][0]
ukb_fmri_data_key = np.squeeze(ukb_fmri_data[0][0][3])
fd_idx = np.where(ukb_fmri_data_key==['meanFD'])[0][0]

# load expression level
miva_data_path = '/data/users2/eswar/temp/MISA_analysis/UKB/postMMIVAresults/C030_UKB_regSite/UKBMRI_MMIVAfull_w0GICA_C030_results_goodSub2907_new_P000.mat'
arrays = {}
f = h5py.File(miva_data_path)
for k, v in f.items():
    arrays[k] = np.array(v)

smri_ref = arrays['icasig'][0][0]
fmri_ref = arrays['icasig'][1][0]
dmri_ref = arrays['icasig'][2][0]

smri_data = np.array(f[smri_ref])
fmri_data = np.array(f[fmri_ref])
dmri_data = np.array(f[dmri_ref])

significant_component=np.array([3, 5, 8, 10, 12, 16, 17])
significant_component=significant_component-1 # change MATLAB indexing to Python indexing

X = np.concatenate( (smri_data[:, significant_component], fmri_data[:, significant_component], dmri_data[:, significant_component]), axis=1 ) 
Y_mean_rm = Y-np.mean(Y)
beta1 = np.linalg.inv( X.T @ X ) @ X.T @ Y

a = 1 / X.shape[1]
delta1 = X @ beta1 - a * Y

Y2 = np.concatenate((np.expand_dims(Y, axis=1), np.expand_dims(Y**2, axis=1), np.expand_dims(Y**3, axis=1), np.expand_dims(ukb_smri_data_array[:, sex_idx], axis=1), np.expand_dims(ukb_fmri_data_array[:, fd_idx], axis=1)), axis=1)
beta2 = np.linalg.inv( Y2.T @ Y2 ) @ Y2.T @ delta1
delta2 = Y2 @ beta2 - delta1
error = np.mean(np.abs(delta2))
print("Brain age error: " + str(round(error, 2)))

np.save('./X.npy', X)
np.save('./Y.npy', Y)
np.save('./Y2.npy', Y2)