#   Functions for loading, standardizing, and calculating correlations (including lagged)
#
#   M Chow, Oct 2012
#   machow@princeton.edu
#

import numpy as np
import nibabel as nib
import os
from itertools import combinations

def standardize(A, axis = -1, demean = True, devar = True):
    """Subtract mean, divide standard deviation (z-scoring).  Operate in-place.

    TODO: fix so that axis doesn't have to be last
          return SD, so that diagonals can be replaced with meaningful information

    """
    if demean:
        A -= np.mean(A, axis = axis)[..., np.newaxis]
    if devar:
        A /= np.std(A, axis = axis)[..., np.newaxis]



def corsubs(A, B, axis = -1, standardized = False):
    """Return correlation mat of time series along specified axis.

    Parameters:
    axis - along which to correlate
    standardized - specify if sequences are already demeaned, with variance of one

    TODO: generalize axis argument

    """
    n = A.shape[axis]
    if standardized: return np.sum(A * B, axis) / n
    else:
        demean = lambda M: M - np.mean(M, axis = axis)[..., np.newaxis]         #TODO here
        return np.sum(demean(A) * demean(B), axis) / ( np.std(A, axis)*np.std(B,axis) * n )


def intersubcorr(C_all, excludeself = True):
    """
    Returns an array of intersubj correlations.  Last two dims must be the cov matrix.

    This function uses that the var( sum of R.V.'s ) is the sum of their cov matrix.

    Parameters:
    excludeself -- is diagonal corr(x_i, x_i)?  remove from calculation?

    """
    covttl = np.apply_over_axes(np.sum, C_all, [-1,-2]).reshape(C_all.shape[:-2])
    covcolsums = C_all.sum(axis=-1)
    N = C_all.shape[-1]
    if excludeself:
        return (covcolsums - 1) / np.sqrt(covttl[...,np.newaxis] - 2*covcolsums)        #TODO replace final +1 with +diagonal of last 2 dims
    else:
        return covcolsums / np.sqrt(covttl)


def crosscor(dlist):
    """Takes list of subject data, returns matrix of correlation matrices at each voxel."""
    map(standardize, dlist)
    N = len(dlist)
    dims = list(dlist[0].shape[:-1]) + [N,N]   #collapsed across time, expanded #subs x #subs
    C_all = np.zeros(dims)
    for ii, jj in combinations(range(N), 2):
        C_all[..., ii, jj] = corsubs(dlist[ii], dlist[jj], standardized = True)

    C_all += C_all.transpose(range(len(dims)-2) + [-1,-2]) + np.eye(C_all.shape[-1])       #Fill in symmetry
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


def trim(A, ends=[0, 0], h=None):
    """Returns Array with last dim trimmed

    Parameters:
    A -- n-dim array
    ends -- how much trim from each end, e.g. (10, -10)
    h -- specify how much to trim from one side.

    """

    if h: ends[h < 0] = h
    return A[..., ends[0] : ends[1] or None]

def shift(A, h, outlen, offset=0):
    """Shifts entire time series by h"""
    h -= offset
    if h < 0: raise BaseException('after adjusting for offset, h is negative.  Cannot keep outlen.')
    return A[..., h : outlen+h]


def load(fname):
    """convenience function to load nib or np filetypes"""
    if os.path.splitext(fname)[1] == '.npy': return np.load(fname)
    else: return nib.load(fname)


def loadData(subdir, mrpaths = [], behpaths = [], condnames = None, subs = None):
    """Returns list with dictionary entry for each subject in subdir

    Parameters:
    subdir -- directory containing subject folders
    mrpaths -- list with paths to fmri nifti files.  Subs in subdir load each path.
    behpath -- list with paths to behavioral data.  Subs in subdir load each path.
    condnames -- (optional) list of names, corresponding to mrpath conditions.
    subs -- (optional) list which subfolders to use instead of taking all in subdir (default)

    Doesn't get data from nifti files initially, in case many subjects are collected (FACTCHECK?)

    """
    if subs is None: subs = os.listdir(subdir)
    if not condnames: condnames = mrpaths
    mrdict = dict(zip(condnames, mrpaths))
    dlist = [dict(Ni = {cond : load(os.path.join(subdir, sub, mrdict[cond])) for cond in mrdict},
                  Sub = sub,
                  Dir = os.path.join(subdir, sub),
                  beh = {behfile : open(os.path.join(subdir, sub, behpaths), 'r').readlines() for behfile in behpaths},
                  mri = None,
                  mrpaths = mrdict,
                  behpaths = behpaths,
                  basedir = os.path.join(subdir, sub))
                for sub in subs]
    return dlist

def loadROI(roifile, binary=True):
    pass
