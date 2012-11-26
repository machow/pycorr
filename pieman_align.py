import numpy as np
import scipy.io as sio
from funcs_correlate import lagcor, shift
from pietools import loadwith
import csv
from os import path

def lagged(Mlist, M_other, lags=None, subsub=False, offset=0):
    """Calculate lagged correlations over an interval"""
    if not lags: lags = int(1/5 * M_other.shape[-1])
    out = []
    newM = M_other
    for M_sub in Mlist:
        if subsub == True: newM = M_other - M_sub
        out.append([float(lagcor(M_sub, newM, h, standardized=False, offset= offset)) for h in lags])
    return out

def csvout(fname, rows):
    f1 = open('lag_' + fname + '.csv', 'w')
    writer = csv.writer(f1)
    writer.writerows(rows)
    f1.close()


func = lambda x: sio.loadmat(x)['subj_roitc'].reshape(300)
pipe_args = dict(subdir = 'pieNDiv/subjects',
                 path = 'analysis/preproc/preproc01.feat/audenv_3mm_gdboth_wordscram_mean_thr0.2_trans_filtered_func_data.mat',
                 func = func)

olga_args = dict(subdir = 'pieNDiv/subsnotpipe',
                 path = 'intact1.feat/audenv_3mm_gdboth_wordscram_mean_thr0.2_trans_filtered_func_data.mat',
                 func = func)

audpath = 'pieNDiv/data/pieman_intact_audenv.mat'
alignfile = 'pieNDiv/intersubj/align_intact01.tsv'
offset=-10
tc_len=280

Mdict = loadwith(**pipe_args)
Mdict.update(loadwith(**olga_args))

if alignfile:
    aligndict = dict([subentry.strip('\n').split('\t') for subentry in open(alignfile).readlines()])
    for key in Mdict.keys():
        h = int(aligndict[key]) if key in aligndict.keys() else 0
        Mdict[key] = shift(Mdict[key], h, tc_len, offset=offset)
    offset = 0

Mlist = Mdict.values()

M_aud = sio.loadmat(audpath)['audenv'].reshape(280)
M_meantc = np.vstack(Mlist).sum(axis=0)
print M_meantc.shape

lagvals = range(-15,15,1)
results = {}
hdr = [['lags'] + Mdict.keys()]
results['Aud_env_cust'] = hdr + zip(lagvals, *lagged(Mlist, M_aud, lagvals, offset))
results['Mean_tc_cust'] = hdr + zip(lagvals, *lagged(Mlist, M_meantc, lagvals, offset))

map(csvout, *zip(*results.viewitems()))
