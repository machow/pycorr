import os
if os.path.isfile('test.h5'): os.remove('test.h5')

import h5py, yaml
import numpy as np
from pieman.funcs_correlate import lagcor
from pieman.workflows import isc_within, isc_between
import pieman.subject as subject

# Setup ---------------------------------------------------

f = h5py.File('test.h5')
E = subject.Exp(f)

config = yaml.load(open('config.yaml'))
E.setup(config, create_conds=True)

# Iterate over conditions ---------------------------------

for condname in config['conds']:
    cond = E.get_cond(condname)
    for run in E.iter_runs(condname):
        run.threshold(cond.attrs['threshold'], save=True)

    # Alignment -------------------------------------------
    a1 = E.f['rois/audenv_3mm_gdboth_wordscram_mean_thr0.2'][...]   #TODO remove hardcoded roi?
    aud_env = cond['audio_env']
    # Calculate lags over range, then set offset
    lags = xrange(-5,5)
    for run in E.iter_runs(condname):
        data = run.load(roi=a1, standardized=True, threshold=True)#slice(None,None,None))
        corrs = [lagcor(data, aud_env, h=h) for h in lags]
        print corrs
        offset = lags[np.argmax(corrs)]
        run.attrs['offset'] = offset

    # With offsets down, can generate composite
    E.gen_composite(condname)

# ISC Calculations ----------------------------------------

cond1 = E.get_cond('wordscram')     # conditions for calculations
cond2 = E.get_cond('wordscram2')

# Run calculations
isc_within(E, cond1)
isc_between(E, cond1, cond2, method='inter-subject')
