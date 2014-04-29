from pieman.tests.gen_corrmat import corr_eig
from pieman.statistics import isc_within_diff, isc_corrmat_within_diff, perm
import numpy as np

def test_perm_test_is_permuting():
    np.random.seed(1)
    fun = lambda  A, B: sum(A) - sum(B)
    out = perm([0,0,0,0], [1,1,1,1], fun, nreps=1000)
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
def test_isc_within_both_equiv():
    assert np.allclose(a, b)

def test_isc_within_correct():
    assert np.allclose(isc_corrmat_within_diff(range(3), range(3), C), 0)

def test_isc_corrmat_within_correct():
    assert np.allclose(isc_within_diff(D[0:3], D[0:3]), 0)

# this would hold if you were doing actual isc not subject-total corr
#def test_isc_within_correct_diff():
#    assert np.allclose(a, -.2)
#
#def test_isc_corrmat_within_correct_diff():
#    assert np.allclose(b, -.2)
    

#TODO PERMUTATION TESTING WITH BOTH FUNCS
# currently using exact correlations, could use mvnorm data?
# i'm not sure what tests would be useful to have here.. maybe this is fine
def test_perm_isc_corrmat_within():
    A, B = range(3), range(3, 6)
    r_null = perm(B, A, isc_corrmat_within_diff, nreps=10000, C=C)
    r = isc_corrmat_within_diff(B, A, C)
    assert np.mean(np.array(r_null) <= r) == 1   #original grouping yields highest corr


C0 = np.zeros([3,3])
C_ones = np.ones([3,3])
C_null_rho1 = np.hstack([np.vstack([C11, C0]), np.vstack([C0, C11])])
C_null_rho0 = np.hstack([np.vstack([C11, C0]), np.vstack([C0, C11])])
offdiag = np.sqrt(C_ones*C11*C22)
C_rho1 = np.hstack([np.vstack([C11, offdiag]), np.vstack([offdiag, C22])])
for M in [C_null_rho1, C_null_rho0, C_rho1]:
    M[np.diag_indices_from(M)] = 1

def test_isc_corrmat_within_null():
    A, B = range(3), range(3, 6)
    r_null = perm(B, A, isc_corrmat_within_diff, nreps=1000, C=C_null_rho1)
    assert np.all(np.array(r_null) == 0) 
