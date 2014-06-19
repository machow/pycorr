"""
Export ROIS from h5 pipe. Will create directories if needed. Files are named "{roi_name}.csv".

E.G. export_roi.py pipe.h5 outputs/rois -c cond1 cond2 -r RSC
"""

import os
import sys
basedir = os.path.dirname(__file__)
sys.path.append(os.path.abspath(basedir + '/../..'))

from pycorr.subject import Exp
from pycorr.pietools import mkdir_p
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = __doc__)
parser.add_argument('h5file',                help='experiment file (something.h5)')
parser.add_argument('-c', '--cond_name', nargs='*',   help='conditions within experiment')
parser.add_argument('-r', '--roi_name',  nargs='*',   help='names of rois to export')
parser.add_argument('out_dir',   nargs='?', help='output directory')
parser.add_argument('-l', '--list',          help='list conditions', action='store_true')
parser.add_argument('-s', '--subs',          help='list subjects', action='store_true')
args = parser.parse_args()

E = Exp(args.h5file)

if args.list:
    print "=======Conditions========"
    for cond in E.f['conds']: print cond
    print "=======ROIS======="
    for roi in E.f['rois']: print roi
    sys.exit()
if args.subs:
    print "=======Subjects========="
    for sub in E.f['subjects']:
        print
        print sub
        for cond in E.f['subjects'][sub].keys():
            print '\t', cond
    sys.exit()

# select conditions and rois (all if none are specified)
conds = [v for k, v in E.f['conds'].iteritems() if not args.cond_name or k in args.cond_name]
rois  = [v for k, v in E.f['rois'].iteritems()  if not args.roi_name  or k in args.roi_name]
    
for cond in conds:
    cond_name= cond.name.split('/')[-1].encode()
    print cond_name
    for roi in rois:
        roi_name = roi.name.split('/')[-1].encode()
        print '\t', roi_name

        roi = roi[...]
        sub_tc = {run.grp.parent.name.split('/')[1] : run.load(roi=roi, threshold=True, standardized=True) for run in E.iter_runs(cond)}
        data = pd.DataFrame(sub_tc)
        data['roi'] = roi_name
        data['cond'] = cond_name

        outname = os.path.join(args.out_dir or '', cond_name, roi_name + '.csv')
        mkdir_p(os.path.dirname(outname))
        data.to_csv(outname)
