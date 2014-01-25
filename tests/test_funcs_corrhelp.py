import numpy as np
from funcs_correlate import standardize
from numpy.testing import assert_almost_equal

"""
Test funcs_correlate.standardize
"""

S1 = np.array([1.,2,3])                         #test float array
S1_inplace = S1.copy()
standardize(S1_inplace, inplace=True)

def test_returns_copy_by_default():              #returns copy by default
    assert S1 is not standardize(S1)

def test_returns_inplace_arg():                  #can work inplace
    tmp_S1  = S1.copy()
    assert tmp_S1 is standardize(tmp_S1, inplace=True)

def test_inplace_equals_copy():                  #equal results for copy and inplace
    print S1
    assert_almost_equal(standardize(S1), S1_inplace)

def test_inplace_int():                          #int arrays never operate inplace
    tmp_S1 = S1.astype(int)
    assert tmp_S1 is not standardize(tmp_S1, inplace=True)
    assert_almost_equal(standardize(tmp_S1), standardize(tmp_S1, inplace=True))

def test_axis_arg():                             #axis arg returns arr with expected shape
    tmp_M = np.array([np.ones(3), np.zeros(3)])
    tmp_stand = standardize(tmp_M, axis=0)
    assert_almost_equal(tmp_stand.shape, tmp_M.shape)

def test_no_demean():
    assert not np.allclose(standardize(S1, demean=False).mean(), 0)

"""
Test funcs_correlate.
"""
