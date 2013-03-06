#   Create a matrix which holds the correlation matrices for each voxel
#
#   M Chow, Oct 2012
#   machow@princeton.edu
#

import numpy as np
from scipy import io as sio
from scipy.stats.stats import nanmean
import nibabel as nib
from funcs_correlate import crosscor, intersubcorr, load_nii_or_npy, standardize
from pietools import loadwith
import yaml

#LOAD CONFIG
condname = ['intact01', 'scram01', 'intact02', 'scram02'][0]
cnfg = yaml.load(open('config.yaml'))
cond = cnfg['cond'][condname]

#INPUTS
subdir = cond['subdir']
nii_synced = cond['filedir'] + cond['nii_synced']       #synced nifti for each sub
roi_files = cnfg['roi']                                 #rois for xcorrs
#OUTPUTS
roi_out = 'intersubj/isc_roi_%s.csv'%condname
ccr_out = 'intersubj/ccr_'+condname
mean_c_out = 'intersubj/ccr_mean_'+condname+'.nii'

#LOAD SUB DATA
ddict = loadwith(subdir, nii_synced, load_nii_or_npy)
sub_indx = sorted(ddict.keys())
dlist = [ddict[sub] for sub in sub_indx]            #list ordered by sub name

#FIND SUBS WITH MEAN TCs < 6000, STANDARDIZE DATA
below_thresh = map(lambda M: M.mean(axis=-1) < 6000, dlist)
for data in dlist: standardize(data)

#ROIs
roi_maps = {roiname : nib.load(fname).get_data() for roiname, fname in roi_files.items()}

ROI_corr = {}
ROI_corrmat = {}
for key, roi in roi_maps.items():
    roi_indx = np.nonzero(roi)                                                      #indices where roi is nonzero
    print 'TC length:', '\t\t', dlist[0][roi_indx].mean(axis=0).shape
    ccr_roi = crosscor([subentry[roi_indx].mean(axis=0) for subentry in dlist], standardized=True)     #crosscors from single ROI timecourse
    ROI_corr[key] = intersubcorr(ccr_roi)
    ROI_corrmat[key] = ccr_roi

#SAVE ROIs
ROI_corr['Folder'] = sub_indx
keys = sorted(ROI_corr.keys())
hdr = np.vstack([key +'.'+condname for key in keys]).T
body = np.vstack([ROI_corr[key] for key in keys]).T
with open(roi_out, 'w') as outfile:
    np.savetxt(outfile, hdr, delimiter=',', fmt='%s')
    np.savetxt(outfile, body, delimiter=',', fmt='%s')

#WHOLE BRAIN SUBxSUB CORRELATION MATRIX
C = crosscor(dlist, standardized = True)
for ii in range(len(dlist)):
    C[below_thresh[ii], ii, :] = np.NaN
    C[below_thresh[ii], :, ii] = np.NaN

#Intersubject Correlations
mean_C = nanmean(intersubcorr(C), axis=-1)
mean_C_nii = nib.Nifti1Image(mean_C, affine=None)
print mean_C.shape
sio.savemat(ccr_out, dict(CCR = C, sub_indx=sub_indx))
nib.save(mean_C_nii, mean_c_out)
