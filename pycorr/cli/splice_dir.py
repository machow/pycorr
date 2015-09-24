import sys, os
import gc
basedir = os.path.dirname(__file__)
sys.path.append(os.path.abspath(basedir + '/../..'))

def splice_dir(dirnames, save, nifti, thresh, check, alldirs, mat, ind_dir='', xyz_shape=None, nan=None):
    from pycorr.pietools import splice_dir, splice_dir_with_ind, copy_nii_hdr, load_nii_or_npy
    import numpy as np

    if alldirs:
        dirs = []
        for dname in dirnames: 
            subdirs = os.walk(dname).next()[1]
            print subdirs
            dirs.extend([os.path.join(dname, fname) for fname in subdirs])
    else:
        dirs = dirnames

    for dirname in dirs:
        if dirname.split('/')[-1] in ['thresh_fail', 'ind'] and alldirs: continue
        print dirname
        if check:
            vols = [int(fname.split('.')[0]) for fname in os.listdir(dirname)]
            print set(range(max(vols))) - set(vols)
            print 'max volume:\t', max(vols)
            continue
        
        # Load and Splice data
        if ind_dir:
            M = splice_dir_with_ind(dirname, ind_dir, xyz_shape)#, mmap='r')
        elif len(os.listdir(dirname)) > 1:
            # splice dir with many (numbered) files
            M = splice_dir(dirname, save=False)#, mmap='r')
        else:
            # just load single file data
            M = load_nii_or_npy(os.path.join(dirname, os.listdir(dirname)[0]))
        print M.shape

        # Threshold
        if thresh: 
            tmp_thresh = load_nii_or_npy(thresh).astype('bool')
            print "thresholding ", tmp_thresh.sum(), " voxels"
            M[tmp_thresh] = np.nan

        # Save in various formats
        if save: np.save(dirname, M)

        if nifti:
            copy_nii_hdr(nifti, M, save=dirname + ".nii.gz", nan=nan)
            #WRITE SAVING LINE HERE

        if mat:
            import scipy.io as sio
            sio.savemat(dirname, {'data': M})
        # Garbage collect, since M can be large
        del M
        gc.collect()






if __name__ == '__main__':
    import sys, os, argparse
    basedir = os.path.dirname(__file__)
    sys.path.append(os.path.abspath(basedir + '/../..'))

    parser = argparse.ArgumentParser()
    parser.add_argument('dirnames', nargs='+', help='directory to reconstruct from parallel analyses')
    parser.add_argument('-s', '--save', action='store_true', help='save output as {directory_name}.npy')
    parser.add_argument('-n', '--nifti', type=str, help='save as nifti, using header from file provided')
    parser.add_argument('-t', '--thresh', type=str, help='threshold mask')
    parser.add_argument('-c', '--check', action='store_true', help='check for missing volumes in folder')
    parser.add_argument('-a', '--alldirs', action='store_true', help='reconstruct all subdirectories')
    parser.add_argument('-m', '--mat', action='store_true', help='save as .mat (in addition to others)')
    parser.add_argument('--ind_dir', type=str)
    parser.add_argument('--xyz_shape', type=int, nargs='*')
    parser.add_argument('--nan', type=float, help='value to substitute for nans in niftis (default leaves them)')
    args = parser.parse_args()

    splice_dir(**args.__dict__)
