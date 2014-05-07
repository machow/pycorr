import os, sys, argparse
basedir = os.path.dirname(__file__)
sys.path.append(os.path.abspath(basedir + '/../..'))
import numpy as np
from pieman.funcs_correlate import crosscor, intersubcorr
from pieman.statistics import perm, isc_corrmat_within_diff
from pieman.subject import Run, Exp
from pieman.pietools import mkdir_p, parse_section, arr_slice

# ARG PARSING
desc = """
Permutation test for within-group ISC differences.
e.g. permute_isc_within.py -t -x 'all' -o test_out
"""
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('-t', help='run test', action='store_true')
parser.add_argument('-a', nargs='*', help='niftis in first group')
parser.add_argument('-b', nargs='*', help='niftis in second group')
parser.add_argument('-x', type=str, default='all', help='slice along row of input arrays to analyze. Can be "all" or slice notation (e.g. ::2)')  
parser.add_argument('-o', '--out', type=str, help='output folder')
parser.add_argument('--thresh', default=6000, help='threshold activation below this level')
parser.add_argument('--n_pass', default=.7, help='number of participants above threshold')
parser.add_argument('--n_reps', type=int, default=1000, help='number of permutations to apply')
args = parser.parse_args()
print args

# TASK ID so script knows where to slice, converts SGE_TASK_ID to 0-indexed
ID = parse_section(args.x) if args.x is not None else int(os.environ['SGE_TASK_ID']) - 1

# OUTPUTS
out = {}

# LOAD FILES
if args.t:                     #TESTING FLAG 
    from pieman.tests.gen_corrmat import fourD
    A_files = fourD + 7000
    B_files = fourD + 7000
elif args.a and args.b:    
    A_files = args.a
    B_files = args.b
else: raise BaseException('need either test or specify inputs')

# Load and Slice data ---------------------------------------------------------
A = [arr_slice(fname, ID) for fname in A_files]
B = [arr_slice(fname, ID) for fname in B_files]

# Thresholding ----------------------------------------------------------------
#Hack to get threshold function, which is a class method TODO def move threshold
import h5py
Run = Run(h5py.File('dummy.h5'))

# threshold tcs with low mean
for dat in A+B: dat[Run.threshold(6000, dat)] = np.nan      
thresh_pass = [~np.isnan(dat.sum(axis=-1)) for dat in A+B]
out['thresh_fail'] = Exp.cond_thresh(thresh_pass, mustpassprop=.7)

# Cross-Correlation matrix (we will permute rows and columns)
out['isc_corrmat'] = crosscor(A+B, standardized=False)

# Combine group indices for correlation matrix (we will shuffle these)
indx_A = range(len(A))
indx_B = range(len(A), len(A + B))
print indx_A
print indx_B

out['null'] = perm(indx_A, indx_B, isc_corrmat_within_diff, nreps=args.n_reps, C = out['isc_corrmat'])
out['r'] = isc_corrmat_within_diff(indx_A, indx_B, out['isc_corrmat'])
out['p'] = np.mean(np.abs(out['r']) <= np.abs(out['null']), axis=-1)
out['isc_A'] = intersubcorr(out['isc_corrmat'][..., indx_A, :][..., :, indx_A])

outtmp = os.path.join(args.out, "{fold}/{ID}.npy")
for k, v in out.iteritems():
    outfile = outtmp.format(fold=k, ID=args.x or os.environ['SGE_TASK_ID'])
    mkdir_p(os.path.dirname(outfile))
    np.save(outfile, v)
