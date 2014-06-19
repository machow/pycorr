#   Functions for loading, standardizing, and calculating correlations (including lagged)
#
#   M Chow, Oct 2012
#   machow@princeton.edu
#

import numpy as np
import nibabel as nib
from itertools import combinations, product
from scipy.stats.stats import nanmean

def standardize(A, axis = -1, demean = True, devar = True, inplace=False):
    """Subtract mean, divide standard deviation (z-scoring).  Operate in-place.
    """
    m  = np.expand_dims(np.mean(A, axis),         axis)
    sd = np.expand_dims( np.std(A, axis, ddof=1), axis)
    if inplace and 'f' in A.dtype.str:                   #make sure A is float
        #inplace operations
        if demean: A -= m
        if devar:  A /= sd
    else:
        #copy operations
        if demean: A = A - m
        if devar:  A = A / sd
    return A



def corsubs(A, B, axis = -1, standardized = False):
    """Return correlation mat of time series along specified axis.

    Parameters:
    axis - along which to correlate
    standardized - specify if sequences are already demeaned, with variance of one

    TODO: generalize axis argument

    """
    n_df = A.shape[axis] - 1
    if standardized: return np.sum(A * B, axis) / n_df
    else:
        demean = lambda M: M - M.mean(axis = axis, keepdims=True)
        std    = lambda M: M.std(axis, ddof=1)
        return np.sum(demean(A) * demean(B), axis) / (std(A)*std(B) * n_df)

def corcomposite(sub, dsummed, **kwargs):
    """Convenience function for correlating subject with a composite it is part of"""
    return corsubs(sub, dsummed-sub, **kwargs)

def sum_tc(dlist, nans = True, standardize_subs=False, standardize_out=False, shape=None):
    """Returns new timecourse from sum of all timecourses.""" 
    newA = np.zeros(shape or dlist[0].shape)
    if standardize_subs: 
        dlist = (standardize(sub, inplace=False) for sub in dlist)
    if nans:
        for entry in dlist: 
            tmp_entry = entry.copy()
            tmp_entry[np.isnan(tmp_entry)] = 0
            newA += tmp_entry
    else:
        for entry in dlist: newA += entry
    if standardize_out: standardize(newA, inplace=True)
    return newA

def intersubcorr(C_all, excludeself = True):
    """
    Returns an array of intersubj correlations.  Last two dims must be the cov matrix.

    This function uses that the var( sum of R.V.'s ) is the sum of their cov matrix.  If entire row / col
    is NaN, then returns ISC as if that subject was removed (to allow dropping voxels with mean activity < 6000.

    Parameters:
    excludeself -- is diagonal corr(x_i, x_i)?  remove from calculation?

    """
    covttl = np.apply_over_axes(np.nansum, C_all, [-1,-2]).reshape(C_all.shape[:-2])
    covcolsums = np.nansum(C_all, axis=-1)
    N = C_all.shape[-1]
    if excludeself:
        return (covcolsums - 1) / np.sqrt(covttl[...,np.newaxis] - 2*covcolsums + 1)        #TODO replace final +1 with +diagonal of last 2 dims (np.diagonal(C_
    else:
        return covcolsums / np.sqrt(covttl)


def crosscor(dlist, B=None, standardized = True):
    """Takes list of subject data, returns matrix of correlation matrices at each voxel."""
    if not np.any(B):
        N = len(dlist)
        dims = list(dlist[0].shape[:-1]) + [N,N]   #collapsed across time, expanded #subs x #subs
        C_all = np.zeros(dims)
        for ii, jj in combinations(range(N), 2):   #just need to loop through upper triangle, then use transpose
            C_all[..., ii, jj] = corsubs(dlist[ii], dlist[jj], standardized = standardized)    #diagonal not filled

        C_all += C_all.transpose(range(len(dims)-2) + [-1,-2]) + np.eye(C_all.shape[-1])       #Fill in symmetry and diag
        return C_all

    else:
        #Second list is provided, calculate all cross-correlations
        #can't use symmetry, so must loop through the cartesian product
        N_a, N_b = len(dlist), len(B)
        dims = list(dlist[0].shape[:-1]) + [N_a, N_b]
        C_all = np.zeros(dims)
        for ii, jj in product(range(N_a), range(N_b)):
            C_all[..., ii, jj] = corsubs(dlist[ii], B[jj], standardized=standardized)
        return C_all




def lagcor(A, B, h, axis=-1, standardized=False, offset=0):
    """Return the corr of A_t+h with B_t.  If dim length differs, cut from end.

    Parameters:
    A -- n-dim array with time as final dim
    B -- same shape array
    axis -- time axis (passed to corsubs)
    standardized -- already demained, and var == one? (passed to corsubs)
    offset -- A is offset number of timepoints passed B

    """

    h -= offset               #time-lock series
    if h < 0:                 # corr(A_x-h, B) == corr(B_x+h, A)
        A,B = B,A
        h = abs(h)

    b_max = A.shape[-1]-h
    a_max = min(B.shape[-1]+h, A.shape[-1])
    if h >= A.shape[-1] or b_max <= 0: raise BaseException('Lag causes zero-length or negative length arrays!')

    return corsubs(A[..., h:a_max], B[...,:b_max], axis, standardized)


def trim(A, ends=(0, 0), h=None):
    """Returns Array with last dim trimmed

    Parameters:
    A -- n-dim array
    ends -- how much trim from each end, e.g. (10, -10)
    h -- specify how much to trim from one side. Overrides that argument in ends.

    """
    ends = list(ends)
    ends[h < 0] = h
    return A[..., ends[0] : ends[1] or None]

def shift(A, h, outlen=None, offset=0):             #TODO remove offset arg?
    """Shifts entire time series by h"""
    h -= offset
    if h < 0: 
        #raise BaseException('after adjusting for offset, h is negative.  Cannot keep outlen.')
        #pad = np.vstack([np.nan]*(-h))
        pad = [np.nan] * (-h)
        A = np.insert(A, 0, pad, axis=-1)
        h = 0
    return A[..., h : outlen and outlen+h]


def roimask(data, roi, filter_func = None, proc_func = None, mean_ts = False):
    """Mask data using values in roi > 0

    Parameters:
    data -- numpy array to be masked
    roi -- array to be checked for nonzero values
    filter_func -- function to further subset roi (e.g. remove time courses with mean < 6000)
    proc_func -- TODO

    """
    if type(roi) == str: roi = nib.load(roi).get_data()
    if len(roi.shape) > 3: roi = roi[:,:,:,0]
    shapes = roi.shape, data.shape
    if roi.shape[:3] != data.shape[:3]: raise BaseException('roi shape: %s \ndata shape: %s'%shapes)
    roi_indx = np.nonzero(roi)
    roi = data[roi_indx]
    if filter_func: roi[filter_func(roi)] = np.nan
    if mean_ts: roi = nanmean(roi, axis=0)
    return roi


def load_roi(nii, thresh=6000, standardize=True, outlen=None, padding=np.mean, meantc=False):
    """Load a nifti file and apply several functions to it.

    nii -- name of nifti to load or ndarray
    standardize -- return all tcs with mean = 0 and stdev = 1
    outlen -- specify how long resultant tcs should be, pad if necessary
    padding -- function to apply to ndarray, must take axis=-1 argument
    meantc -- return the meantc over all voxels (before evaluating standardize)

    TODO: fix default outlen arg

    """
    if not hasattr(nii, 'shape'): nii = nib.load(nii).get_data()
    nii[nii.mean(axis=-1) < thresh] = np.nan

    if meantc: out = nanmean(nii.reshape([-1, nii.shape[-1]]), axis=0)  #USE apply_over_axes instead?
    else: out = nii

    if standardize: out = (out - out.mean(axis=-1)[..., np.newaxis]) / out.std(axis=-1, ddof=1)[..., np.newaxis] #this is unruly, use the current standardize func
    outdiff = outlen - out.shape[-1]
    if outdiff >= 0:
        if hasattr(padding, 'func_name'): padding = padding(out, axis=-1)[..., np.newaxis]    #call padding if it is a function
        out = np.concatenate([out] + [padding]*outdiff, axis=-1)
    else: 
        out = out[..., :outlen] 
    return out
