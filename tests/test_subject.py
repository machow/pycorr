import os
if os.path.isfile('test.h5'): os.remove('test.h5')

import inspect
import h5py, yaml
import subject
reload(subject)

f = h5py.File('test.h5')
E = subject.Exp(f)

config = yaml.load(open('_newconfig.yaml'))
try:
    E.setup(config)
    for run in E.iter_runs('intact'):
        run.threshold(6000, save=True)
except:
    var = inspect.trace()[-1][0]#[-1][0]
    raise
