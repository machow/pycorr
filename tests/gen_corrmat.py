import numpy as np
import sklearn.decomposition as deco
from numpy import random

def corr_eig(lam, nsubs, npoints, C=None):
    """ Generates an npoints x nsubs matrix of data that produces correlation matrix C

    args:
        lam     - off-diagonal correlations
        nsubs   - columns of output mat (number of dvs)
        npoints - rows of output mat (number of observations)
    """
    if not np.any(C):
        # Make correlation matrix
        C = np.identity(nsubs)
        tu = np.triu_indices(nsubs,k=1)  # indices for upper triangle
        C[tu]       = lam
        C[tu[::-1]] = lam

    #Exceptions for correlations of 1

    # Spectral decomposition
    w, V = np.linalg.eigh(C)
    W = np.diag(np.sqrt(w))     # sqrt(eigen values) on diagonal
    U = np.dot(V, W).T          # Transformation matrix (U.T %*% U is C)

    # Generate data
    data = random.randn(npoints, nsubs)
    ortho = deco.PCA(nsubs).fit(data).transform(data)                           # Principal components (so uncorrelated)
    orthonorm = (ortho - np.mean(ortho, axis=0)) / ortho.std(axis=0, ddof=0)    # Z-scored

    return C, np.dot(orthonorm, U)

# nsubs     x corrs x thresh x nan
corrs = np.arange(-.4,1.1,.1)
nsubs = 3
npoints = 10
ttl_pnts = nsubs * len(corrs)* 3 * npoints
out = []
sol = []
for lam in corrs:
    C, data = corr_eig(lam, nsubs, npoints)
    print np.corrcoef(data.T)
    sol.append(C)
    out.append(data)
data

from funcs_correlate import crosscor
np.array(sol)
out = np.array(out).transpose([2,0,1])
out_sol = crosscor(out, standardized=False)
np.allclose(out_sol, np.array(sol))

for ii, sub in enumerate(out): print ii

fourD = out[[(0,1,2)*2*2]].reshape([3,2,2,15,10])
# subject dim x (indx 1 is below thresh) x (indx 1 has two nan subs) x (time)
out_match = fourD + 7000
out_nan = fourD.copy()            #nan for two subjects
out_nan[0:2,...] = np.nan
out_thresh = fourD.copy()         #below thresh

np.allclose(crosscor(out_match), out_sol)
np.allclose(crosscor(fourD, standardized=False), out_sol)
