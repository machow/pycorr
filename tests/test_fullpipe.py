import tempfile
import subprocess
import os
import shutil
import sys; sys.path.append('../..')

import h5py, yaml
from pycorr import subject
from pycorr.workflows import isc_within
import numpy as np

class test_fullpipe(object):
    config = yaml.load(open('examples/test/config.yaml'))
    def setup(self):
        self.orig_path = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        setup_fname = os.path.abspath('setup_example.py')
        self.perm_fname = os.path.abspath('pycorr/cli/permute_isc_within.py')

        os.chdir(self.tmpdir)
        subprocess.call(['python', setup_fname, 'test'])
        print os.listdir('.')

        f = h5py.File('test.h5')
        E = subject.Exp(f)

        print self.config
        E.setup(self.config, create_conds=True)
        cond = E.get_cond('test_cond')
        for run in E.iter_runs(cond):
            run.threshold(cond.attrs['threshold'], save=True)

        E.gen_composite('test_cond')
        isc_within(E, cond)

        E.gen_cond_thresh('test_cond')
        E.f.file.flush()
        self.E = E

    def teardown(self):
        shutil.rmtree(self.tmpdir)
        os.chdir(self.orig_path)

    def test_correct_reference_setting(self):
        arr_type = type(self.E.iter_runs('test_cond').next().data)
        print arr_type
        if self.config['conds']['test_cond'].get('reference'):
            assert arr_type is subject.AttrArray
        else:
            assert arr_type is h5py._hl.dataset.Dataset

    def test_isc_mat(self):
        sol = np.load('pipeline/subjects/test_sol.npy')
        C = self.E.get_cond('test_cond')['correlations/test_cond/isc_mat'][...]
        assert np.allclose(C, sol)

    def test_perm_sim_serial(self):
        perm_fname = self.perm_fname
        sol = np.load('pipeline/subjects/test_sol.npy')

        subprocess.call('python %s -t -x all -o test_out --n_reps 10'%perm_fname, shell=True)
        C = np.load('test_out/isc_corrmat/all.npy')
        indx_A = range(3)
        assert np.allclose(C[..., indx_A,:][...,:,indx_A], sol)

    def test_perm_sim_parallel(self):
        perm_fname = self.perm_fname
        sol = np.load('pipeline/subjects/test_sol.npy')
        for ii, _ in enumerate(sol):
            print ii
            subprocess.call('python %s -t -x %s -o test_out --n_reps 10'%(perm_fname, ii), shell=True) #run analysis
        subprocess.call('python pycorr/cli/splice_dir.py -s test_out/isc_corrmat', shell=True)                        #splice back together
        C = np.load('test_out/isc_corrmat.npy')
        indx_A = range(3)
        assert np.allclose(C[..., indx_A,:][...,:,indx_A], sol)

    def test_perm_serial(self):
        perm_fname = self.perm_fname
        sol = np.load('pipeline/subjects/test_sol.npy')

        subprocess.call('python %s --hdf5 test.h5 -a test_cond -x all -o test_out --isc_only'%perm_fname, shell=True)
        C = np.load('test_out/isc_corrmat/all.npy')
        indx_A = range(3)
        assert np.allclose(C[..., indx_A,:][...,:,indx_A], sol)

    def test_perm_parallel(self):
        perm_fname = self.perm_fname
        sol = np.load('pipeline/subjects/test_sol.npy')
        for ii, _ in enumerate(sol):
            print ii
            subprocess.call('python %s --hdf5 test.h5 -a test_cond -x %s -o test_out --isc_only'%(perm_fname, ii), shell=True) #run analysis
        subprocess.call('python pycorr/cli/splice_dir.py -s test_out/isc_corrmat', shell=True)                        #splice back together
        C = np.load('test_out/isc_corrmat.npy')
        indx_A = range(3)
        assert np.allclose(C[..., indx_A,:][...,:,indx_A], sol)

from copy import deepcopy
class test_fullpipe_reference(test_fullpipe): pass
newcnfg = deepcopy(test_fullpipe_reference.config)
newcnfg['conds']['test_cond']['reference'] = True
test_fullpipe_reference.config = newcnfg

class test_fullpipe_nomaxlen(test_fullpipe): pass
newcnfg = deepcopy(test_fullpipe.config)
newcnfg['conds']['test_cond']['max_len'] = None
test_fullpipe_reference.config = newcnfg
