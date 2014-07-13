from pycorr.stats.boot import circ, bootstrap_ts_indexes, b_star, ts_boot
import numpy as np

def test_circ_same_val():
    a = np.arange(10)
    circ_indx = a.copy()
    circ(circ_indx, 4, len(circ_indx))
    assert np.all(a[a] == a[circ_indx])

def test_circ_window_1():
    a = np.arange(10)
    circ(a, 1, len(a))
    assert np.all(np.unique(a) == a)

def test_b_star_nan_entry():
    assert np.isnan(b_star(np.array([1,2,3,4,np.nan])))

def test_b_star_no_var():
    assert np.isnan(b_star(np.array([1,1,1,1,1])))

def bootstrap_ts_indexes_uniform(npoints, l):
    """Ensure that items sampled are uniform for each position in bootstrapped samples"""
    np.random.seed(1)
    a = np.arange(npoints)
    uniform = 1. / a.size
    err = .01
    indx = bootstrap_ts_indexes(a, l, n_samples=20000)
    S = a[indx]
    for pos in range(S.shape[-1]):
        proportions = np.bincount(S[:,pos]) / float(S[:,pos].size)
        assert np.all(np.abs(proportions - uniform) < err)

def test_bootstrap_ts_indexes_uniform():
    npoints = 4
    for window_length in range(1, npoints + 1):
        yield bootstrap_ts_indexes_uniform, npoints, window_length

def test_ts_boot_mean():
    npoints = 200
    dvs = 2
    data = np.random.normal(loc=1, size=[2, npoints])
    out = np.zeros([10000, dvs])
    ts_boot(data, np.mean, dvs, out=out, axis=-1)
    assert np.all(np.abs(out.mean(axis=0) - 1) < .1)

from pycorr.stats.boot import run_boot_within_isc_diff, calc_mean_isc, calc_mean_isc2
from pycorr.gen_corrmat import corr_eig
def test_ts_boot_isc_within():
    """Test A_isc > B_isc w/data that yields exact correlation matrix."""

    _, A = corr_eig(.64, 4, 100)
    _, B = corr_eig(.49, 4, 100)
    A, B = A.T, B.T
    res = run_boot_within_isc_diff(A, B, l=2, n_samples=100, 
                             out_arr=None)
    assert np.allclose(calc_mean_isc(A), calc_mean_isc2(A))
    assert res['p_gt'] < .05
    

#rpy2 tests for b_star
    # 1 dimensional
    # 2 dimensional
