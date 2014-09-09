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

import os, argparse

def permute_isc_within(a, b, x, outfile, mask='', meth=None, hdf5=None, thresh=6000, n_pass=.7, n_reps=10000, t=False, kwargs=None):
    if not kwargs: kwargs = {}
    import numpy as np
    from pycorr.funcs_correlate import crosscor, intersubcorr
    from pycorr.subject import Run, Exp
    from pycorr.pietools import mkdir_p, parse_section, arr_slice
    # TASK ID so script knows where to slice, converts SGE_TASK_ID to 0-indexed
    ID = parse_section(x) if x is not None else int(os.environ['SGE_TASK_ID']) - 1

    # OUTPUTS
    out = {}

    # Load and Slice data ---------------------------------------------------------
    mask = np.load(mask) if mask else slice(None)
    if not hdf5:
        # LOAD FILES
        if t:                     #TESTING FLAG 
            from pycorr.gen_corrmat import fourD
            A_files = fourD + 7000
            B_files = fourD + 7000
        elif a and b:    
            A_files = [os.path.join(a[0], fname) for fname in os.listdir(a[0])]  #TODO change back, hack until rondo jobs are fixed
            B_files = [os.path.join(b[0], fname) for fname in os.listdir(b[0])]
        else: raise BaseException('need either test or specify inputs')

        A = [arr_slice(fname, ID)[...,mask].astype('float') for fname in A_files]
        B = [arr_slice(fname, ID)[...,mask].astype('float') for fname in B_files]
        # Thresholding
        #Hack to get threshold function, which is a class method TODO def move threshold
        import h5py
        Run = Run(h5py.File('dummy.h5'))

        # threshold tcs with low mean
        for dat in A+B: dat[Run.threshold(thresh, dat)] = np.nan      
        thresh_pass = [~np.isnan(dat.sum(axis=-1)) for dat in A+B]
        out['thresh_fail'] = Exp.cond_thresh(thresh_pass, mustpassprop=n_pass)
    else:
        E = Exp(hdf5)
        A = [run.load(use_subset=mask, standardized=True, threshold=True,  _slice=ID) for run in E.iter_runs(a[0])]
        if b:  #TODO fix, so hacky.. this script needs structure (want to let arg.b be optional
            B = [run.load(standardized=True, threshold=True, _slice=ID) for run in E.iter_runs(b[0])]
        else: B = []
        out['thresh_fail'] = E.get_cond(a[0])['threshold'][...]

    # Combine group indices for correlation matrix (we will shuffle these) --------
    indx_A = range(len(A))
    indx_B = range(len(A), len(A + B))
    print indx_A
    print indx_B

    # Cross-Correlation matrix (we will permute rows and columns) -----------------
    out['isc_corrmat'] = crosscor(A+B, standardized=False)
    out['isc_A'] = intersubcorr(out['isc_corrmat'][..., indx_A, :][..., :, indx_A])
    out['isc_B'] = intersubcorr(out['isc_corrmat'][..., indx_B, :][..., :, indx_B])

    # Permutation Test ------------------------------------------------------------
    if meth == 'perm':
        from pycorr.stats.perm import run_perm
        res = run_perm(indx_A, indx_B, out['isc_corrmat'], n_reps)
        out.update(res)

    # Bootstrap Test ----------------------------------------------------------
    elif meth == 'boot':
        from pycorr.stats.boot import run_boot_within_isc_diff
        res = run_boot_within_isc_diff(A, B, n_samples=n_reps, **kwargs)
        out.update(res)

    # Output ----------------------------------------------------------------------
    outtmp = os.path.join(outfile, "{fold}/{ID}.npy")
    for k, v in out.iteritems():
        outfile = outtmp.format(fold=k, ID=x or ID)
        mkdir_p(os.path.dirname(outfile))
        np.save(outfile, v)


if __name__ == '__main__':
    import yaml 

    # ARG PARSING
    parser = argparse.ArgumentParser(description= __doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', help='run test', action='store_true')
    parser.add_argument('-a', nargs='*', help='niftis in first group or condname if hdf5')
    parser.add_argument('-b', nargs='*', help='niftis in second group or condname if hdf5')
    parser.add_argument('-x', type=str, help='slice along row of input arrays to analyze. Can be "all" or slice notation (e.g. ::2)')  
    parser.add_argument('-o', '--outfile', type=str, help='output folder')
    parser.add_argument('-m', '--mask', type=str, help='boolean timecourse mask for subsetting')
    parser.add_argument('--meth', default='boot', help='type of test (perm or bootstrap)')
    parser.add_argument('--hdf5', nargs='?', help='hdf5 pipeline to load niftis from')
    parser.add_argument('--thresh', default=6000, help='threshold activation below this level. (not implemented,  hardcoded)')
    parser.add_argument('--n_pass', default=.7, help='number of participants above threshold. (not implemented, hardcoded)')
    parser.add_argument('--n_reps', type=int, default=10000, help='number of permutations to apply')
    parser.add_argument('--kwargs', type=yaml.load, help="""additional arguments to pass to test function in YAML format. (e.g. "{a: 1, b: two}")""")
    args = vars(parser.parse_args())
    print args

    permute_isc_within(**args)
