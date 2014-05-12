import sys, os
basedir = os.path.dirname(__file__)
sys.path.append(os.path.abspath(basedir + '/../..'))
from pieman.pietools import splice_dir, copy_nii_hdr, load_nii_or_npy
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('dirname', type=str, help='directory to reconstruct from parallel analyses')
parser.add_argument('-s', '--save', action='store_true', help='save output as {directory_name}.npy')
parser.add_argument('-n', '--nifti', type=str, help='save as nifti, using header from file provided')
parser.add_argument('-t', '--thresh', type=str, help='threshold mask')
args = parser.parse_args()

M = splice_dir(args.dirname, save=False)
if args.thresh: 
    thresh = load_nii_or_npy(args.thresh)
    M[thresh] = np.nan

if args.save: np.save(args.dirname, M)

if args.nifti:
    nii = copy_nii_hdr(args.nifti, M, save=args.dirname + ".nii.gz")


