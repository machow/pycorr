#   Create a matrix which holds the correlation matrices for each voxel
#
#   M Chow, Oct 2012
#   machow@princeton.edu
#

import numpy as np
from scipy import io as sio
import nibabel as nib
from funcs_correlate import loadData, crosscor, intersubcorr

#load data
cond = 'intact01'
fname = 'trans_filtered_func_data_sync.npy'
pipeargs = dict(subdir = 'pieNDiv/subjects',
                mrpaths = ['analysis/preproc/preproc01.feat/'+fname],
                condnames = [cond],
                behpaths = [])

olgaargs = dict(subdir = 'pieNDiv/subsnotpipe',
                mrpaths = ['intact1.feat/'+fname],
                condnames = [cond],
                behpaths = [])

subdata = loadData(**pipeargs) + loadData(**olgaargs)       #list of subject dicts

#Get a dictionary with sub : data dict
ddict = {sub['Sub'] : sub['Ni'].values()[0] for sub in subdata}
sub_indx = sorted(ddict.keys())
dlist = [ddict[sub] for sub in sub_indx]

#ROIs
base = 'pieNDiv/data/'
ROI_files = dict(precun = base + 'precun_3mm_thr50.nii',
                 bhipp = base + 'BHipp_3mm_thr30.nii',
                 mpfc = base + 'MPFC_150.nii',
                 rsc = base + 'RSC.nii',
                 aud = base + 'audenv_3mm_gdboth_wordscram_mean_thr0.2.nii')
ROI_maps = {roiname : nib.load(fname).get_data() for roiname, fname in ROI_files.items()}

ROI_corr = {}
ROI_corrmat = {}
for key, roi in ROI_maps.items():
    roi_indx = np.nonzero(roi)
    print dlist[0][roi_indx].mean(axis=0).shape
    ccr_roi = crosscor([subentry[roi_indx].mean(axis=0) for subentry in dlist])     #crosscors from single ROI timecourse
    N = ccr_roi.shape[-1]
    ROI_corr[key] = intersubcorr(ccr_roi)
    ROI_corrmat[key] = ccr_roi


ROI_corr['Folder'] = sub_indx
#Save ROIs
keys = sorted(ROI_corr.keys())
hdr = np.vstack([key +'.'+cond for key in keys]).T
body = np.vstack([ROI_corr[key] for key in keys]).T
with open('pieman_intersubj_%s_out.csv'%cond, 'w') as outfile:
    np.savetxt(outfile, hdr, delimiter=',', fmt='%s')
    np.savetxt(outfile, body, delimiter=',', fmt='%s')



#Whole Brain Correlation matrix
#This will likely run on the cluster without taking up too much memory, but can be segmented spatially if needed
C = crosscor(dlist)

#Intersubject Correlations (mean off-diagonal crosscorrelation)
mean_C = intersubcorr(C).mean(axis=-1)
mean_C_nii = nib.Nifti1Image(mean_C, affine=None)
print mean_C.shape
sio.savemat('pieNDiv/intersubj/ccr_'+cond, dict(CCR = C, sub_indx=sub_indx))
nib.save(mean_C_nii, 'pieNDiv/intersubj/ccr_mean_'+cond+'.nii')
