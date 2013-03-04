import numpy as np
import scipy.io as sio
from funcs_correlate import lagcor, shift
from pietools import loadwith
import csv
import yaml

def lagged(Mlist, M_other, lags=None, subsub=False, offset=0):
    """Return list with lagged correlations for each entry in Mlist

    Parameters:
    Mlist - list of subject time courses
    M_other - time course to compare all in Mlist to
    lags - list of lag values to derive xcorrs for
    subsub - remove subject timecourse from M_other before calculating
    offset - Mlist is offset number of time points past M_other

    """

    if not lags: lags = int(1/5 * M_other.shape[-1])            #default number of lags
    out = []
    newM = M_other
    for M_sub in Mlist:
        if subsub == True: newM = M_other - M_sub
        out.append([float(lagcor(M_sub, newM, h, standardized=False, offset= offset)) for h in lags])
    return out

def csvout(fname, rows):
    f1 = open(fname + '.csv', 'w')
    writer = csv.writer(f1)
    writer.writerows(rows)
    f1.close()


#LOAD CONFIG
cnfg = yaml.load(open('config.yaml'))
condname = ['intact01', 'scram01', 'intact02', 'scram02'][0]                   #select condition
cond = cnfg['cond'][condname]

#INPUTS
subdir = cond['subdir']
nii_aud = cond['filedir'] + cond['nii_aud']             #A1 nifti for each sub
shiftfile = cond['shift_tc'] + condname + '.tsv'        #How much to shift each sub tc
aud_env = cond['aud_env']                               #Auditory envelope tc
outbase = cond['align']                                 #base name for outputs
###

func = lambda x: sio.loadmat(x)['subj_roitc'].reshape(300)
Mdict = loadwith(subdir, nii_aud, func)   #load subject Aud1 timecourses

#SHIFT WITH ALIGN FILE
offset=-10
tc_len=280
if shiftfile:
    aligndict = dict([subentry.strip('\n').split('\t') for subentry in open(shiftfile).readlines()])
    for key in Mdict.keys():
        h = int(aligndict[key]) if key in aligndict.keys() else 0           #don't shift subjects not in file
        Mdict[key] = shift(Mdict[key], h, tc_len, offset=offset)
    offset = 0

#LOAD TIME COURSES FOR XCORR
Mlist = Mdict.values()                                          #subject tc's
M_aud = sio.loadmat(aud_env)['audenv'].reshape(280)     #audio envelope
M_meantc = np.vstack(Mlist).sum(axis=0)                         #mean tc

#LAGGED CORRELATIONS
lagvals = range(-15,15,1)
results = {}
hdr = [['lags'] + Mdict.keys()]

results[outbase + '_aud_env_' + condname] = hdr + zip(lagvals, *lagged(Mlist, M_aud, lagvals, offset = offset))
results[outbase + '_mean_tc_' + condname] = hdr + zip(lagvals, *lagged(Mlist, M_meantc, lagvals, offset = offset))

map(csvout, *zip(*results.viewitems()))
