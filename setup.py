import os, shutil, sys
from pycorr.pietools import mkdir_p
from os.path import join
import os.path

basedir = os.path.dirname(os.path.abspath(__file__))

if len(sys.argv) > 1:
    example_dir = sys.argv[1]
    shutil.copy2(join(basedir, 'examples', example_dir, 'config.yaml'), 'config.yaml')
    shutil.copy2(join(basedir, 'tests/run_pipe.py'), 'run.py')
    os.symlink(join(basedir, 'pycorr'), 'pycorr')

    if example_dir == 'full': 
        print "setting up full example (using nifti files and ROIS)"
        import urllib, zipfile
        print "downloading data"
        datfile = 'pycorr_full_data.zip'
        datpath = os.path.join(basedir, 'examples', example_dir, datfile)
        # Get zipped data if needed, unzip in example folder
        if not os.path.isfile(datpath): 
            urllib.urlretrieve('http://princeton.edu/~machow/pycorr_full_data.zip', datpath)
        zipfile.ZipFile(open(datpath, 'rb')).extractall()

    elif example_dir == 'test':
        import numpy as np
        from pycorr.gen_corrmat import fourD, fourD_sol
        mkdir_p('pipeline/subjects')
        for ii, sub in enumerate(fourD):
            np.save('pipeline/subjects/test_sub%s.npy'%ii, sub + 7000)
        np.save('pipeline/subjects/test_sol.npy', fourD_sol)
