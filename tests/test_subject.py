import os
if os.path.isfile('test.h5'): os.remove('test.h5')

import h5py, yaml
import numpy as np
from funcs_correlate import lagcor
import workflows
reload(workflows)
from workflows import isc_within, isc_between
import subject
reload(subject)

f = h5py.File('test.h5')
E = subject.Exp(f)

config = yaml.load(open('_newconfig.yaml'))
E.setup(config, create_conds=True)
for condname in config['conds']:
    cond = E.get_cond(condname)
    for run in E.iter_runs(condname):
        run.threshold(cond.attrs['threshold'], save=True)


    #Alignment
    a1 = E.f['rois/audenv_3mm_gdboth_wordscram_mean_thr0.2'][...]
    aud_env = cond['audio_env']

    lags = xrange(-15,16)
    for run in E.iter_runs(condname):
        data = run.load(roi=a1, standardized=True, threshold=True)#slice(None,None,None))
        corrs = [lagcor(data, aud_env, h=h) for h in lags]
        print corrs
        offset = lags[np.argmax(corrs)]
        run.attrs['offset'] = offset
    
    E.gen_composite(condname)

#SECOND RUN for testing
cond1 = E.get_cond('wordscram')
cond2 = E.get_cond('wordscram2')


#test isc
isc_within(E, cond1)
isc_between(E, cond1, cond2, method='inter-subject')
