import tempfile
import subprocess
import os
import shutil
import sys; sys.path.append('../..')

import h5py, yaml
from pieman import subject
from pieman.workflows import isc_within, isc_between
import numpy as np

#TODO wrap in class, setup method for all copying and pipe setup
orig_path = os.getcwd()
try: 
    tmpdir = tempfile.mkdtemp()
    setup_fname = os.path.abspath('setup.py')
    perm_fname = os.path.abspath('scripts/permute_isc_within.py')
    os.chdir(tmpdir)
    subprocess.call(['python', setup_fname, 'test'])
    print os.listdir('.')
    #subprocess.call('python run.py', shell=True)

    f = h5py.File('test.h5')
    E = subject.Exp(f)

    config = yaml.load(open('config.yaml'))
    E.setup(config, create_conds=True)
    cond = E.get_cond('test_cond')
    for run in E.iter_runs(cond):
        run.threshold(cond.attrs['threshold'], save=True)

    E.gen_composite('test_cond')
    isc_within(E, cond)

    sol = np.load('pipeline/subjects/test_sol.npy')
    C = cond['correlations/test_cond/isc_mat'][...]
    assert np.allclose(C, sol)

    #Permutation test test
    subprocess.call('python %s -t -x all -o test_out --n_reps 10'%perm_fname, shell=True)
    C = np.load('test_out/isc_corrmat/all.npy')
    indx_A = range(3)
    assert np.allclose(C[..., indx_A,:][...,:,indx_A], sol)

finally:
    shutil.rmtree(tmpdir)
    os.chdir(orig_path)
