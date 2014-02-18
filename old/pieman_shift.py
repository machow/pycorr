import numpy as np
import nibabel as nib
from funcs_correlate import shift, load_nii_or_npy
from pietools import loadwith
from os import path
import yaml

#CONFIG
condname = ['intact01', 'scram01', 'intact02', 'scram02'][0]
cnfg = yaml.load(open('config.yaml'))
cond = cnfg['cond'][condname]

#INPUTS
subdir = cond['subdir']
alignf = cond['shift_tc'] + condname + '.tsv'
nii_trans = cond['filedir'] + cond['nii_trans']
nii_out = '_sync'

offset=-10
tc_len=280

#LOAD DATA
Mdict = loadwith(subdir, nii_trans, load_nii_or_npy)   #load subject data

#LOAD SHIFT FILE
aligndict = dict([subentry.strip('\n').split('\t') for subentry in open(alignf).readlines()])
for key in aligndict: aligndict[key] = int(aligndict[key])

#SHIFT TIME COURSE AND SAVE (appended with _sync)
for key in Mdict.keys():
    h = int(aligndict[key]) if key in aligndict.keys() else 0           #don't shift subjects not in file
    Mdict[key] = shift(Mdict[key], h, tc_len, offset)
    f, ext = path.splitext(nii_trans)
    print 'writing %s to: \t\t %s'%(key, path.join(f + '_sync'))
    nib.save(nib.Nifti1Image(Mdict[key], affine=None), path.join(subdir, key, f + nii_out))
    #print 'writing %s to: \t\t %s'%(key, path.join(f + '_sync'))
