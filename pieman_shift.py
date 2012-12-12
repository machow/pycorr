import numpy as np
from funcs_correlate import shift, loadData
from os import path


cond = 'scram01'
alignf = 'pieNDiv/intersubj/align_' + cond + '.tsv'
offset=-10
tc_len=280

pipeargs = dict(subdir = 'pieNDiv/subjects',
                mrpaths = ['analysis/preproc/preproc02.feat/trans_filtered_func_data.nii'],
                condnames = [cond],
                behpaths = [])

olgaargs = dict(subdir = 'pieNDiv/subsnotpipe',
                mrpaths = ['scrambled1.feat/trans_filtered_func_data.nii'],
                condnames = [cond],
                behpaths = [])

#load data
subdata = loadData(**pipeargs) + loadData(**olgaargs)

#Load alignment file
aligndict = dict([subentry.strip('\n').split('\t') for subentry in open(alignf).readlines()])
for key in aligndict: aligndict[key] = int(aligndict[key])

#Shift time course and save (appended with _sync)
for subdict in subdata:
    data = subdict['Ni'][cond].get_data()
    h = aligndict[subdict['Sub']] if (subdict['Sub'] in aligndict.keys()) else 0
    data = shift(data, h, tc_len, offset)
    f, ext = path.splitext(subdict['mrpaths'][cond])
    np.save(path.join(subdict['basedir'], f + '_sync'), data)

