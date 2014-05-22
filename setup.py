import os, shutil, sys
from pietools import mkdir_p

basedir = os.path.dirname(os.path.abspath(__file__))
shutil.copy2(os.path.join(basedir, 'config.yaml'), 'config.yaml')

if len(sys.argv) > 1 and sys.argv[1] == 'test':
    try:
        shutil.copytree(os.path.join(basedir, 'tests/data/fullpipe'), 'pipeline')
    except: print "can't find test data..."
    shutil.copy2(os.path.join(basedir, 'tests/run_pipe.py'), 'run.py')
    os.symlink(basedir, 'pieman')
    
    import numpy as np
    sys.path.append(os.path.join(basedir, '..'))
    from pieman.tests.gen_corrmat import fourD, fourD_sol
    mkdir_p('pipeline/subjects')
    for ii, sub in enumerate(fourD):
        np.save('pipeline/subjects/test_sub%s.npy'%ii, sub + 7000)
    np.save('pipeline/subjects/test_sol.npy', fourD_sol)
