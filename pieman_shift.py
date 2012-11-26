import numpy as np
from funcs_correlate import shift, loadData
from os import path


pipeargs = dict(subdir = 'pieNDiv/subjects',
                mrpaths = ['analysis/preproc/preproc01.feat/trans_filtered_func_data.nii'],
                condnames = ['intact01'],
                behpaths = [])

olgaargs = dict(subdir = 'pieNDiv/subsnotpipe',
                mrpaths = ['intact1.feat/trans_filtered_func_data.nii'],
                condnames = ['intact01'],
                behpaths = [])

subdata = loadData(**pipeargs) + loadData(**olgaargs)


aligndict = dict([subentry.strip('\n').split('\t') for subentry in open('pieNDiv/intersubj/align.tsv').readlines()])
for key in aligndict: aligndict[key] = int(aligndict[key])
offset=-10
tc_len=280

for subdict in subdata:
    data = subdict['Ni']['intact01'].get_data()
    h = aligndict[subdict['Sub']] if (subdict['Sub'] in aligndict.keys()) else 0
    data = shift(data, h, tc_len, offset)
    f, ext = path.splitext(subdict['mrpaths']['intact01'])
    np.save(path.join(subdict['basedir'], f + '_sync'), data)

