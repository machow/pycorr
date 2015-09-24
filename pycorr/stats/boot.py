from numpy.random import randint
from scipy.stats import norm
import numpy as np

def circ(start, l, maxlen):
    """Change start index in-place, so end indexes can grab next l elements in circular manner"""
    indx_edge = maxlen - start < l
    start[indx_edge] = start[indx_edge] - maxlen

def bootstrap_ts_indexes(data, l, n_samples=10000, method='circular'):
    """
    Given data points data, where axis 0 is considered to delineate points, return
    an array where each row is a set of bootstrap indexes. This can be used as a list
    of bootstrap indexes as well.

    Parameters:
        data (ndarray): data to be resampled on axis 0
        l:              time window
        n_samples:      number of bootstrap samples to take
        method:         timeseries method (fixed, moving, circular). Only circular is implemented.

    Returns: index array with dim (n_samples X data.shape[0])
    """
    npoints = data.shape[0]
    if method == 'circular': 
        N = int(np.ceil(npoints / float(l)))
        blocks = randint(npoints, size=(n_samples, N))
        circ(blocks, l, maxlen=npoints)
        expand_block = lambda block: np.array([block + ii for ii in xrange(l)]).flatten(order='F')
        indexes = np.array([expand_block(block) for block in blocks])
    return indexes[..., :npoints]


from statsmodels.tsa.stattools import acf, ccovf

def lam(s):
  return (np.abs(s)>=0)*(np.abs(s)<0.5) + 2*(1-np.abs(s))*(np.abs(s)>=0.5)*(np.abs(s)<=1)

# PORT of b.star function from R
def b_star(data, Kn = None, mmax = None, Bmax = None, c = None,
           round = True):
    if len(data.shape) == 1: data = data.reshape([-1,1])
    elif len(data.shape) > 2: raise Exception('data may not be greater than 2d')

    n, k = data.shape
    if Kn   is None: Kn   = max(5, np.ceil(np.log(n)))
    if mmax is None: mmax = np.ceil(np.sqrt(n)) + Kn
    if Bmax is None: Bmax = np.ceil(min(3 * np.sqrt(n), n/3.))
    if c    is None: c    = norm.ppf(0.975)
    BstarCB = np.empty(k)         # TODO: Risky to initialize as empty?
    for ii in range(k):
        rho_k = acf(data[:,ii])[1:]
        rho_k_crit = c * np.sqrt(np.log10(n)/n)
        insig = np.abs(rho_k) < rho_k_crit
        num_insig = np.array([np.sum(insig[jj:jj + Kn]) for jj in range(int(mmax-Kn+1))])
        
        # Determine mhat
        if np.any(num_insig == Kn): 
            # mhat is indx of first with Kn insig
            mhat = np.where(num_insig == Kn)[0][0] + 1   #[indx_tuple][arr_pos]
        elif (np.abs(rho_k) > rho_k_crit).sum() == 1:
            # mhat is indx of max rho_k greater than rho_k_crit
            mhat = np.where(np.abs(rho_k) > rho_k_crit)[0][0] + 1
        else: mhat = 1
        
        M = min(2*mhat, mmax)
        kk = np.arange(-M, M+1.)
        
        ccov_right = ccovf(data[:,ii], data[:,ii], unbiased=False)[:M+1]       # M+1 to include lag 0
        R_k = np.concatenate([ccov_right[-1:0:-1], ccov_right])
        Ghat = np.sum(lam(kk / M) * np.abs(kk) * R_k)
        DCBhat = 4./3 * np.sum(lam(kk/M) * R_k)**2

        BstarCB[ii] = (2 * Ghat**2 / DCBhat)**(1/3.) * n**(1/3.)

    if round: 
        BstarCB[BstarCB > Bmax] = Bmax
        BstarCB[BstarCB < 1] = 1
        np.round(BstarCB, out=BstarCB)
    
    return BstarCB

def ts_boot(dlist, func, l, n_samples=10000, method='circular', out=None, indx_file='', **kwargs):
    """Return func results from bootstrapped samples of data
    
    Parameters:
        data (ndarray): data to draw bootstrap samples from last dim
        func:           func to run over new samples
        n_samples:      number of bootstrapped samples to draw
        kwargs:         additional arguments for func

    Returns: 
        array with results of func call over each bootstrapped sample,
        has dimension [n_samples x dim(func(dlist))].
    """
    assert len(np.unique([data.shape[-1] for data in dlist])) == 1    #All have same number timepoints
    
    example_tc = dlist[0].T
    if out is None: out = [None] * n_samples
    
    if not indx_file: indexes = bootstrap_ts_indexes(example_tc, l, n_samples)
    else: indexes = np.load(indx_file)[:n_samples]

    for ii, boot_indx in enumerate(indexes):
        sample = [data[..., boot_indx] for data in dlist]
        out[ii] = func(sample, **kwargs)
    return np.array(out, copy=False)

from pycorr.funcs_correlate import intersubcorr, crosscor, corcomposite, sum_tc, standardize
from scipy.stats import nanmean
def calc_isc_cormat(dlist, mean=True): 
    """calculate within group ISC by first deriving cross correlation matrix"""
    dlist = [standardize(sub, inplace=False) for sub in dlist]
    isc = intersubcorr(crosscor(dlist, standardized=True))
    if mean: return nanmean(isc, axis=-1)
    else:    return isc

def calc_isc_subttl(dlist, mean=True): 
    """standardize subjects, then correlate against sum of others."""
    dlist = [standardize(sub, inplace=False) for sub in dlist]
    dsummed = sum_tc(dlist, nans = True, standardize_subs = False)
    isc = [corcomposite(sub, dsummed) for sub in dlist]
    if mean: return nanmean(isc, axis=0)
    else:
        # Make a copy and transpose array. This is not very efficient,
        # and should not be used for extremely time-sensitive computation.
        tmp = np.array(isc)
        return tmp.transpose(range(1,tmp.ndim) + [0])  # put first axis last

def run_boot_within_isc_diff(A, B, l, n_samples, out_arr=None, indx_file=''):
    out = {}

    out_shape = (n_samples, ) + A[0].shape[:-1]      #n_reps x spatial_dims
    out_arr = np.zeros(out_shape, dtype=float)
    swap_dims = range(1,len(out_shape)) + [0]                        #list with first and last dims swapped

    out['distA'] = ts_boot(A, calc_isc_subttl, l, n_samples=n_samples, out = out_arr.copy(), indx_file=indx_file)
    out['distB'] = ts_boot(B, calc_isc_subttl, l, n_samples=n_samples, out = out_arr.copy(), indx_file=indx_file)
    # swap axis with correlations to be last dim
    for k in ['distA', 'distB']: out[k] = out[k].transpose(swap_dims)
    # since 1 corr, add axis for broadcasting
    out['r'] = (calc_isc_subttl(A) - calc_isc_subttl(B))[..., np.newaxis]
    # calc one-sided p-values against null that diff is 0 
    out['p_ltq'] = (out['distA'] - out['distB'] >= 0).mean(axis=-1)
    out['p_gt'] = (out['distA'] - out['distB'] < 0).mean(axis=-1)
    return out
