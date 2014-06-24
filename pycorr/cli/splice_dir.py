import sys, os
basedir = os.path.dirname(__file__)
sys.path.append(os.path.abspath(basedir + '/../..'))

def splice_dir(dirnames, save, nifti, thresh, check, alldirs, mat):
    from pycorr.pietools import splice_dir, copy_nii_hdr, load_nii_or_npy
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
        if 'thresh_fail' in dirname and alldirs: continue
        print dirname
        if check:
            vols = [int(fname.split('.')[0]) for fname in os.listdir(dirname)]
            print set(range(max(vols))) - set(vols)
            print 'max volume:\t', max(vols)
            continue

        M = splice_dir(dirname, save=False, mmap='r')
        print M.shape
        if thresh: 
            thresh = load_nii_or_npy(thresh)
            M[thresh] = np.nan

        if save: np.save(dirname, M)

        if nifti:
            nii = copy_nii_hdr(nifti, M, save=dirname + ".nii.gz")
            #WRITE SAVING LINE HERE

        if mat:
            import scipy.io as sio
            sio.savemat(dirname, {'data': M})






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
    args = parser.parse_args()

    splice_dir(**args.__dict__)
