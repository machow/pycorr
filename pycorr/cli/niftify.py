from pycorr.pietools import copy_nii_hdr
import argparse
from os import path

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=copy_nii_hdr.__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('nii')
    parser.add_argument('data')
    parser.add_argument('--save', default='')
    parser.add_argument('--neglog', action='store_true', default=False)
    parser.add_argument('--nan', type=int, default=None)

    args = parser.parse_args()
    if not args.save: args.save = path.splitext(args.data)[0] + '.nii.gz'

    copy_nii_hdr(**vars(args))
