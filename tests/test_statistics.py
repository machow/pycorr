from statistics import isc_within_diff, perm_test
import numpy as np

def test_perm_test_is_permuting():
    np.random.seed(1)
    fun = lambda  A, B: sum(A) - sum(B)
    out = perm_test([0,0,0,0], [1,1,1,1], fun, nreps=10000)
    assert .010 < np.mean(np.array(out) == 4) < .018

def test_isc_within_equiv():
    assert False
