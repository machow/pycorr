"""
Permutation test for within-group ISC differences. May be used to calculate within-group ISC in parallel.

Outputs:
    thresh_fail -- participants that fall below threshold
    isc_corrmat -- full subject-subject correlation matrix.
    isc_A -- subject-total correlation for group A.
    null -- distribution from perm test (n_reps x spatial_dims)
    r -- within_A - within_B
    p -- proportion of null equal to or more extreme (two-tailed) than r

E.G. permute_isc_within.py -t -x 'all' -o test_out
"""

import os, sys, argparse
basedir = os.path.dirname(__file__)
sys.path.append(os.path.abspath(basedir + '/../..'))
import numpy as np
from pycorr.funcs_correlate import crosscor, intersubcorr
from pycorr.statistics import perm, isc_corrmat_within_diff
from pycorr.subject import Run, Exp
from pycorr.pietools import mkdir_p, parse_section, arr_slice

# ARG PARSING
parser = argparse.ArgumentParser(description= __doc__)
parser.add_argument('-t', help='run test', action='store_true')
parser.add_argument('-a', nargs='*', help='niftis in first group or condname if hdf5')
parser.add_argument('-b', nargs='*', help='niftis in second group or condname if hdf5')
parser.add_argument('-x', type=str, help='slice along row of input arrays to analyze. Can be "all" or slice notation (e.g. ::2)')  
parser.add_argument('-o', '--out', type=str, help='output folder')
parser.add_argument('-m', '--mask', type=str, help='boolean timecourse mask for subsetting')
parser.add_argument('--isc_only', action='store_true', help='just calculate intersubject correlation, instead of running perm test')
parser.add_argument('--hdf5', nargs='?', help='hdf5 pipeline to load niftis from')
parser.add_argument('--thresh', default=6000, help='threshold activation below this level. (not implemented,  hardcoded)')
parser.add_argument('--n_pass', default=.7, help='number of participants above threshold. (not implemented, hardcoded)')
parser.add_argument('--n_reps', type=int, default=1000, help='number of permutations to apply')
args = parser.parse_args()
print args

# TASK ID so script knows where to slice, converts SGE_TASK_ID to 0-indexed
ID = parse_section(args.x) if args.x is not None else int(os.environ['SGE_TASK_ID']) - 1

# OUTPUTS
out = {}

# Load and Slice data ---------------------------------------------------------
mask = np.load(args.mask) if args.mask else slice(None)
if not args.hdf5:
    # LOAD FILES
    if args.t:                     #TESTING FLAG 
        from pycorr.gen_corrmat import fourD
        A_files = fourD + 7000
        B_files = fourD + 7000
    elif args.a and args.b:    
        A_files = [os.path.join(args.a[0], fname) for fname in os.listdir(args.a[0])]  #TODO change back, hack until rondo jobs are fixed
        B_files = [os.path.join(args.b[0], fname) for fname in os.listdir(args.b[0])]
    else: raise BaseException('need either test or specify inputs')

    A = [arr_slice(fname, ID)[...,mask].astype('float') for fname in A_files]
    B = [arr_slice(fname, ID)[...,mask].astype('float') for fname in B_files]
    # Thresholding
    #Hack to get threshold function, which is a class method TODO def move threshold
    import h5py
    Run = Run(h5py.File('dummy.h5'))

    # threshold tcs with low mean
    for dat in A+B: dat[Run.threshold(6000, dat)] = np.nan      
    thresh_pass = [~np.isnan(dat.sum(axis=-1)) for dat in A+B]
    out['thresh_fail'] = Exp.cond_thresh(thresh_pass, mustpassprop=.7)
else:
    E = Exp(args.hdf5)
    A = [run.load(use_subset=mask, standardized=True, threshold=True,  _slice=ID) for run in E.iter_runs(args.a[0])]
    if args.b:  #TODO fix, so hacky.. this script needs structure (want to let arg.b be optional
        B = [run.load(standardized=True, threshold=True, _slice=ID) for run in E.iter_runs(args.b[0])]
    else: B = []
    E.get_cond(args.a[0])
    out['thresh_fail'] = E.get_cond(args.a[0])['threshold'][...]

# Combine group indices for correlation matrix (we will shuffle these) --------
indx_A = range(len(A))
indx_B = range(len(A), len(A + B))
print indx_A
print indx_B

# Cross-Correlation matrix (we will permute rows and columns) -----------------
out['isc_corrmat'] = crosscor(A+B, standardized=False)
out['isc_A'] = intersubcorr(out['isc_corrmat'][..., indx_A, :][..., :, indx_A])

# Permutation Test ------------------------------------------------------------
if not args.isc_only:
    out_shape = (args.n_reps, ) + out['isc_corrmat'].shape[:-2]      #n_reps x spatial_dims
    swap_dims = range(1,len(out_shape)) + [0]                        #list with first and last dims swapped
    out['null'] = perm(indx_A, indx_B, isc_corrmat_within_diff, C = out['isc_corrmat'],
                       nreps=args.n_reps, out=np.zeros(out_shape))
    out['null'] = out['null'].transpose(swap_dims)                            #put corrs on last dim
    out['r'] = isc_corrmat_within_diff(indx_A, indx_B, out['isc_corrmat'])[..., np.newaxis] #since 1 corr, add axis for broadcasting
    out['p'] = np.mean(np.abs(out['r']) <= np.abs(out['null']), axis=-1)

# Output ----------------------------------------------------------------------
outtmp = os.path.join(args.out, "{fold}/{ID}.npy")
for k, v in out.iteritems():
    outfile = outtmp.format(fold=k, ID=args.x or ID)
    mkdir_p(os.path.dirname(outfile))
    np.save(outfile, v)
