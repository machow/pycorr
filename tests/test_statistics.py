from pieman.tests.gen_corrmat import corr_eig
from pieman.statistics import isc_within_diff, isc_corrmat_within_diff, perm_test
import numpy as np

def test_perm_test_is_permuting():
    np.random.seed(1)
    fun = lambda  A, B: sum(A) - sum(B)
    out = perm_test([0,0,0,0], [1,1,1,1], fun, nreps=10000)
    assert .010 < np.mean(np.array(out) == 4) < .018

C11 = np.ones([3,3]) * .3
C22 = np.ones([3,3]) * .5
C12 = np.ones([3,3]) * 0 
C = np.hstack([np.vstack([C11, C12]), np.vstack([C12, C22])])
C[np.diag_indices_from(C)] = 1

D = corr_eig(None, 6, 20, C)[1].T
a = isc_within_diff(D[0:3], D[3:])
b = isc_corrmat_within_diff(range(3), range(3, 6), C)
#TODO wrap in setup func or class
def test_isc_within_equiv():
    assert np.allclose(a,b)

def test_isc_within_correct():
    assert np.allclose(a, C)

def test_isc_corrmat_within_correct():
    assert np.allclose(b, C)

#TODO PERMUTATION TESTING WITH BOTH FUNCS

