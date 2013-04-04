# Look at specific ROIs within the intersubject cross-correlation matrix
#
#   M Chow, Nov 2012
#   machow@princeton.edu
#
#NOTE: THIS HAS BEEN PULLED INTO THE pieman_intersubj.py file

import nibabel as nib
import numpy as np
from scipy import io as sio

base = 'pieNDiv/data/'
mat = sio.loadmat('pieNDiv/intersubj/ccr_intact01.mat')
CCR, sub_indx = mat['CCR'], mat['sub_indx']
N = CCR.shape[-1]

ROI_files = dict(precun = base + 'precun_3mm_thr50.nii',
                 bhipp = base + 'BHipp_3mm_thr30.nii',
                 mpfc = base + 'MPFC_150.nii',
                 rsc = base + 'RSC.nii',
                 aud = base + 'audenv_3mm_gdboth_wordscram_mean_thr0.2.nii')

ROI_maps = {roiname : nib.load(fname).get_data() for roiname, fname in ROI_files.items()}

ROI_corr = {}
for key, roi in ROI_maps.items():
    roi_indx = np.nonzero(roi)
    ccr_roi = (CCR[roi_indx].sum(axis=-1) - 1) / (N-1)           #subtract 1 for diagonal, take 1 from N for corr w/self
    ROI_corr[key] = ccr_roi

