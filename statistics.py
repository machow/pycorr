import numpy as np
from funcs_correlate import corcomposite, sum_tc, intersubcorr
from scipy.stats import nanmean

def isc_within_diff(A, B, standardized=False):
    isc = lambda L, ttl: nanmean([corcomposite(dat, ttl, standardized=standardized) for dat in L],
                                 axis=0)
    A_composite = sum_tc(A)
    B_composite = sum_tc(B)
    A_mean_isc = isc(A, A_composite)
    B_mean_isc = isc(B, B_composite)
    return A_mean_isc - B_mean_isc

def isc_corrmat_within_diff(indxA, indxB, C):
    C_A = C[indxA][:,indxA]
    C_B = C[indxB][:,indxB]
    return nanmean(intersubcorr(C_A), axis=-1) - nanmean(intersubcorr(C_B), axis=-1)

def perm_test(A, B, fun, nreps = 1, out = None,  **kwargs):
    """Permutation test. Randomly shuffles group labels, then runs fun. Group
    sizes are preserved.

    Parameters:
        A, B - lists with elements to permute across groups
        fun  - function of form fun(new_A, new_B, [opt, ... ,])
        nreps- number of repetitions
        out  - optional container for results (e.g. numpy array with dtype)
        kwargs - passed to fun

    
    """
    if out is None: out = [None] * nreps
    AB = A + B
    for ii in range(nreps):
        np.random.shuffle(AB)
        out[ii] = fun(AB[:len(A)], AB[len(A):], **kwargs)

    return out


