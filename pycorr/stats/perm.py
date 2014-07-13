import numpy as np
from pycorr.funcs_correlate import corcomposite, sum_tc, intersubcorr
from scipy.stats import nanmean

def isc_within_diff(A, B, standardized=False):
    """Contrast within-group subject-total correlation for A and B.

    This function operates on the timecourse data, so is slower
    than isc_corrmat_within_diff. Inputs may be multi-dimensional.
    The last dimension is used for correlations (e.g. time should be last).

    Arguments:
        A (list): List of timecourse data for each member of group A.
        B (list): Timecourses of same length as A.

    Returns:
        ndarray with isc for A minus isc for B.

    """

    isc = lambda L, ttl: nanmean([corcomposite(dat, ttl, standardized=standardized) for dat in L],
                                 axis=0)
    A_composite = sum_tc(A)
    B_composite = sum_tc(B)
    A_mean_isc = isc(A, A_composite)
    B_mean_isc = isc(B, B_composite)
    return A_mean_isc - B_mean_isc

def isc_corrmat_within_diff(indxA, indxB, C):
    """Faster within-group subject-total correlation contrast using correlation matrix

    Arguments:
        indxA (list): list of indices corresponding to group A members.
        indxB (list): likewise for group B (should be no overlap)

    Returns:
        ndarray with isc for A minus isc for B.

    """
    C_A = C[..., np.vstack(indxA), np.hstack(indxA)]   # last rows and columns using indxA
    C_B = C[..., np.vstack(indxB), np.hstack(indxB)]
    return nanmean(intersubcorr(C_A), axis=-1) - nanmean(intersubcorr(C_B), axis=-1)

def perm(A, B, fun, nreps = 1, out = None,  **kwargs):
    """Permutation test. Randomly shuffles group labels, then runs fun. Group
    sizes are preserved.

    Parameters:
        A: lists with elements (or indices) to permute across groups
        B: similar to A, but other group
        fun: function of form fun(new_A, new_B, [opt1, ... ,])
        nreps: number of repetitions
        out: optional container for results (e.g. numpy array with dtype)
        kwargs: optional parameters passed to fun
    
    """
    if out is None: out = [None] * nreps
    AB = A + B
    for ii in range(nreps):
        np.random.shuffle(AB)
        out[ii] = fun(AB[:len(A)], AB[len(A):], **kwargs)

    return out

def run_perm(indx_A, indx_B, C, n_reps):
    out = {}

    out_shape = (n_reps, ) + C.shape[:-2]      #n_reps x spatial_dims
    out_arr = np.zeros(out_shape, dtype=float)
    swap_dims = range(1,len(out_shape)) + [0]                        #list with first and last dims swapped

    out['null'] = perm(indx_A, indx_B, isc_corrmat_within_diff, C = C,
                    nreps=n_reps, out=out_arr)
    out['null'] = out['null'].transpose(swap_dims)                            #put corrs on last dim
    out['r'] = isc_corrmat_within_diff(indx_A, indx_B, C)[..., np.newaxis] #since 1 corr, add axis for broadcasting
    out['p'] = np.mean(np.abs(out['r']) <= np.abs(out['null']), axis=-1)
    return out
