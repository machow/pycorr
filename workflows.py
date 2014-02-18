import numpy as np
from scipy.stats import nanmean
from funcs_correlate import crosscor, corcomposite

def dset_overwrite(g, name, dat):
    """If dset doesn't exist, create it. Then write dat to dset"""
    dset =  g.require_dataset(name, shape=dat.shape, dtype=dat.dtype, compression='gzip', track_times=True)
    dset[...] = dat


def isc_within(E, condA, method=('inter-subject', 'subject-total'), threshold=True):
    condnameA = condA.name.split('/')[-1]
    g_name = 'correlations/%s'%condnameA
    g_out = condA[g_name] if g_name in condA else condA.create_group(g_name)
    
    iter_runs = E.iter_runs(condnameA)
    dlist = [run.load(standardized=True, threshold=threshold) for run in iter_runs]
    if 'inter-subject' in method:
        #wg isc
        C = crosscor(dlist, standardized=True)
        C = C - np.diag([np.nan]*C.shape[-1])           #not done in place, just in case
        C_mean = np.squeeze(np.apply_over_axes(nanmean, C, [-1,-2]))

        dset_overwrite(g_out, 'isc_mat', C)
        dset_overwrite(g_out, 'inter-subject', C_mean)

    if 'subject-total' in method:
        #WG subject-total
        sub_ttl_corr = [corcomposite(dat, condA['composite']) for dat in dlist]
        sub_ttl_mean = nanmean(sub_ttl_corr, axis=0)
        
        dset_overwrite(g_out, 'subject-total', sub_ttl_mean)

def isc_between(E, condA, condB, method=('inter-subject', 'subject-total', 'total-total'), threshold=True):
    condnameA = condA.name.split('/')[-1]
    condnameB = condB.name.split('/')[-1]
    g_name = 'correlations/%s'%condnameB
    g_out = condA[g_name] if g_name in condA else condA.create_group(g_name)

    iter_runs = E.iter_runs(condnameA)
    A_dlist = [run.load(standardized=True, threshold=threshold) for run in iter_runs]
    if 'inter-subject' in method:
        #wg isc
        B_runs = E.iter_runs(condnameB)
        B_dlist = [run.load(standardized=True, threshold=threshold) for run in B_runs]
        C = crosscor(A_dlist, B_dlist, standardized=True)
        C_mean = np.squeeze(np.apply_over_axes(nanmean, C, [-1,-2]))

        dset_overwrite(g_out, 'isc_mat', C)
        dset_overwrite(g_out, 'inter-subject', C_mean)

    if 'subject-total' in method:
        #WG subject-total
        sub_ttl_corr = [corcomposite(dat, condB['composite']) for dat in A_dlist]
        sub_ttl_mean = nanmean(sub_ttl_corr, axis=0)

        dset_overwrite(g_out, 'subject-total', sub_ttl_mean)
