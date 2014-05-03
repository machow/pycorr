import os, sys, gc, argparse
basedir = os.path.dirname(__file__)
sys.path.append(os.path.abspath(basedir + '/../..'))
import numpy as np
import nibabel as nib
from pieman.funcs_correlate import crosscor, intersubcorr
from pieman.statistics import perm, isc_corrmat_within_diff
from pieman.subject import Run, Exp

def slice(nii_file, indx):
    "slices indx from first dim from matrix, for parallelization"
    print 'loading data'
    if type(nii_file) is str:
        nii = nib.load(nii_file)
        dat = nii.get_data()
        del nii
    else: dat = nii_file
    out = dat[indx:indx+1].copy()
    del dat                   #garbage collection
    gc.collect()
    return out

# ARG PARSING
parser = argparse.ArgumentParser()
parser.add_argument('-t', help='run test', action='store_true')
parser.add_argument('-a', nargs='*', help='niftis in first group')
parser.add_argument('-b', nargs='*', help='niftis in second group')
parser.add_argument('-x', type=int, help='slice along row of input arrays to analyze')  
parser.add_argument('-o', '--out', type=str, help='output folder')
parser.add_argument('--thresh', default=6000, help='threshold activation below this level')
parser.add_argument('--n_pass', default=.7, help='number of participants above threshold')
args = parser.parse_args()
print args

# TASK ID so script knows where to slice
ID = args.x if args.x else int(os.environ['SGE_TASK_ID']) - 1

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

#Hack to get threshold function, which is a class method TODO def move threshold
import h5py
Run = Run(h5py.File('dummy.h5'))
#
A = [slice(fname, ID) for fname in A_files]
B = [slice(fname, ID) for fname in B_files]

# threshold tcs with low mean
for dat in A+B: dat[Run.threshold(6000, dat)] = np.nan      
thresh_pass = [~np.isnan(dat.sum(axis=-1)) for dat in A+B]
out['thresh_fail'] = Exp.cond_thresh(thresh_pass, mustpassprop=.7)

#
out['isc_corrmat'] = crosscor(A+B, standardized=False)

indx_A = range(len(A))
indx_B = range(len(A), len(A + B))
print indx_A
print indx_B

out['null'] = perm(indx_A, indx_B, isc_corrmat_within_diff, nreps=100, C = out['isc_corrmat'])
out['r'] = isc_corrmat_within_diff(indx_A, indx_B, out['isc_corrmat'])
out['p'] = np.mean(np.abs(out['r']) <= np.abs(out['null']), axis=-1)
out['isc_A'] = intersubcorr(out['isc_corrmat'][..., indx_A, :][..., :, indx_A])

outtmp = os.path.join(args.out, "{fold}/{ID}.npy")
for k, v in out.iteritems():
    outfile = outtmp.format(fold=k, ID=ID)
    try: os.mkdir(os.path.dirname(outfile))
    except: pass
    np.save(outfile, v)
