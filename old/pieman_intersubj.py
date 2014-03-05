#   Create a matrix which holds the correlation matrices for each voxel
#
#   M Chow, Oct 2012
#   machow@princeton.edu
#

import numpy as np
from scipy import io as sio
from scipy.stats.stats import nanmean
import nibabel as nib
from funcs_correlate import crosscor, intersubcorr, load_nii_or_npy, corsubs, sum_tc, standardize
from pietools import loadwith
import yaml

def rois(dlist, sub_indx, roi_maps, roi_out, roi_mat_out, ind_tc, save_rois):
    """Calculate ISC and xcorr matrices for all ROIs, across last dimension

    Parameters:
    dlist - list of numpy arrays of equal dimensions
    sub_indx - names corresponding to arrays in dlist (for output)
    roi_maps - dict containing regions of interest
    roi_out - name for ISC output (csv)
    roi_mat_out - TODO name for xcorr output (csv)
    ind_tc - optional, time course with same dim as in dlist for ISC and xcorr
    """
    ROI_corr = {}
    ROI_corrmat = {}
    for key, roi in roi_maps.items():
        roi_indx = np.nonzero(roi)                                                      #indices where roi is nonzero
        roi_tcs = [nanmean(subentry[roi_indx], axis=0) for subentry in dlist]
        if np.any(ind_tc):
            roi_ind_tc = nanmean(ind_tc[roi_indx], axis=0)
            ROI_corr[key] = np.array([corsubs(sub_tc, roi_ind_tc, standardized=False) for sub_tc in roi_tcs])
        else:
            ccr_roi = crosscor(roi_tcs, standardized=False)     #crosscors from single ROI timecourse
            ROI_corr[key] = intersubcorr(ccr_roi)
            ROI_corrmat[key] = ccr_roi
        if save_rois:
            sio.savemat(save_rois+key, {sub : sub_roitc for sub, sub_roitc in zip(sub_indx, roi_tcs)})

    #SAVE ROIs
    ROI_corr['Folder'] = sub_indx
    keys = sorted(ROI_corr.keys())
    hdr = np.vstack([key +'.'+condname for key in keys]).T
    body = np.vstack([ROI_corr[key] for key in keys]).T
    with open(roi_out, 'w') as outfile:
        np.savetxt(outfile, hdr, delimiter=',', fmt='%s')
        np.savetxt(outfile, body, delimiter=',', fmt='%s')

    #SAVE ROI CORR MATS
    for key, roi in ROI_corrmat.items():
        with open('%s_%s.csv'%(roi_mat_out, key), 'w') as outfile:
            np.savetxt(outfile, np.vstack(sub_indx).T, delimiter=',', fmt='%s')
            np.savetxt(outfile, roi, delimiter=',', fmt='%s')

def isc(dlist, sub_indx, thresh=6000, mustpassprop=.7, ccr_out="", mean_c_out="", mean_out="", supersubj=None, super_out=""):
    """Calculate the [subject x subject] xcorr matrix for each voxel, as well as mean ISC per voxel.

    Parameters:
    dlist - list of numpy arrays of equal dimension
    sub_indx - names corresponding to arrays in dlist (for output)
    ccr_out - name for xcorr matrix output (.mat)
    mean_c_out - name for mean ISC output (.nii)
    mean_out - name for mean timecourse at each voxel
    supersubj - nifti for supersubject to do between group ISC
    """
    #FIND SUBS WITH MEAN TCs < 6000, STANDARDIZE DATA
    n_max = (mustpassprop) * len(dlist)
    above_thresh = map(lambda M: M.mean(axis=-1) >= thresh, dlist)
    thresh_fail = np.sum(above_thresh, axis=0) < n_max 			#sum num of failed subjects per voxel
    print 'n_max: ', n_max
    for entry in above_thresh:
        print entry.shape
    thresh_fail.shape

    if mean_c_out or ccr_out:
        #WHOLE BRAIN SUBxSUB CORRELATION MATRIX
        C = crosscor(dlist, standardized = False)
        for ii in range(len(dlist)):
            C[~above_thresh[ii], ii, :] = np.NaN                 #make entries that didn't meet threshold NA
            C[~above_thresh[ii], :, ii] = np.NaN

        #Intersubject Correlations
        mean_C = nanmean(intersubcorr(C), axis=-1)
        mean_C[thresh_fail] = np.NaN
        mean_C_nii = nib.Nifti1Image(mean_C, affine=None)
        print mean_C.shape
        if ccr_out: sio.savemat(ccr_out, dict(CCR = C, sub_indx=sub_indx))
        if mean_c_out: nib.save(mean_C_nii, mean_c_out)

    if mean_out:
        standardize(dlist, inplace=True)
        mean = sum_tc(dlist, standardize_out=True)
        mean[thresh_fail] = np.NaN
        nib.save(nib.Nifti1Image(mean, affine=None), mean_out)


    #CORRELATE EACH SUB WITH SUPERSUB TC
    if np.any(supersubj):
        ISC_list = map(lambda X: corsubs(X, supersubj, standardized=False), dlist)
        ISC_out = nanmean(ISC_list, axis=0)             #Mean ISC per voxel
        ISC_out[thresh_fail] = np.NAN
        nib.save(nib.Nifti1Image(ISC_out, affine=None), super_out)

if __name__ == '__main__':
    import sys
    #LOAD CONFIG
    condname = sys.argv[1]
    cnfg = yaml.load(open('config.yaml'))
    cond = cnfg['cond'][condname]

    #INPUTS
    subdir = cond['subdir']
    nii_synced = cond['filedir'] + cond['nii_synced']           #synced nifti for each sub

    #LOAD SUB DATA
    ddict = loadwith(subdir, nii_synced, load_nii_or_npy)
    sub_indx = sorted(ddict.keys())
    dlist = [ddict[sub] for sub in sub_indx]            #list ordered by sub name

    #LOAD SUPERSUBJECT (if using)
    if cond['supersubj']:
        ind_tc = load_nii_or_npy(cond['supersubj'])         	    #independent timecourse if specified
        is_ind = '_ind'
    else:
        ind_tc = False
        is_ind = ''

    #ROI Analysis
#    roi_map = {roiname : nib.load(fname).get_data() for roiname, fname in cnfg['roi'].items()}
#    roi_out = 'intersubj/isc_roi_%s%s.csv'%(condname, is_ind)
#    roi_mat_out = 'intersubj/isc_roi_mat_%s%s'%(condname, is_ind)
#    rois(dlist, sub_indx, roi_map, roi_out, roi_mat_out, ind_tc, False)#, 'analysis/alph_rois/')

    #Whole-Brain Analysis
    out = dict(
                ccr_out = 'intersubj/ccr_'+condname,
                mean_c_out = 'intersubj/ccr_mean_'+condname+'.nii',
                mean_out = 'data/supersubj/mean_'+condname+'.nii',
                super_out = 'intersubj/ccr_mean_supersubj_'+condname+'.nii')
    isc(dlist, sub_indx, supersubj=ind_tc, **out)

    #save mean time course for each voxel
#    mean_nii = nib.Nifti1Image(sum_tc(dlist), affine=None)              #modifies dlist in place
#    nib.save(mean_nii, 'data/whole_brain_tc_%s.nii.gz'%condname)
