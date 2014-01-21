# Three ways of calculating ISC
# 1) operating on correlation matrix
# 2) subtract each sub from summed tc and take correlation
# 3) correlate with summed tc, use correction given by Wherry

from nose import with_setup
import numpy as np
from numpy.testing import assert_almost_equal
from funcs_correlate import standardize, corsubs, crosscor, intersubcorr

dims = (2,2, 10)
nsubs = 3

subs = [np.random.random(dims) for ii in range(nsubs)]
for M in subs: M[0,0] = range(dims[-1])
for M in subs: standardize(M)
subs[0][1,1] = np.NAN

C_all = crosscor(subs, standardized=True)
C_all[1,1,0] = np.NAN
isc1 = intersubcorr(C_all)

M_ttl = np.nansum(subs, axis=0)
isc2 = np.array([corsubs(M, M_ttl-M) for M in subs]).transpose([1,2,0])

isc3_list = []
for M in subs:
    r_all = corsubs(M, M_ttl)
    s_all = np.std(M_ttl, axis=-1)
    s_i = np.std(M, axis=-1)
    M_cors = (r_all*s_all - s_i) / \
            np.sqrt(s_i**2 + s_all**2 - 2*s_i*s_all*r_all) #wherry formula
    isc3_list.append(M_cors)
isc3 = np.array(isc3_list).transpose([1,2,0])

def test_intersubcorrXmeantc():
    assert_almost_equal(isc1, isc2)

def test_intersubcorrXwherry():
    assert_almost_equal(isc1, isc3)

def test_meantcXwherry():
    assert_almost_equal(isc2, isc3)



