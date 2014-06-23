"""
Export ROIS from h5 pipe. Will create directories if needed. Files are named "{roi_name}.csv".

E.G. export_roi.py pipe.h5 outputs/rois -c cond1 cond2 -r RSC
"""

#import os
#import sys
#basedir = os.path.dirname(__file__)
#sys.path.append(os.path.abspath(basedir + '/../..'))

from pycorr.subject import Exp
from pycorr.pietools import mkdir_p
import pandas as pd
import os, sys
import argparse

def export_roi(h5file, conds, rois, out_dir):
    E = Exp(h5file)

    # select conditions and rois (all if none are specified)
    conds = [v for k, v in E.f['conds'].iteritems() if not conds or k in conds]
    rois  = [v for k, v in E.f['rois'].iteritems()  if not rois  or k in rois]
        
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

            outname = os.path.join(out_dir or '', cond_name, roi_name + '.csv')
            mkdir_p(os.path.dirname(outname))
            data.to_csv(outname)

def summarize(h5file, conds, subs):
    E = Exp(h5file)

    if conds:
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

if __name__ == '__main__':
    import inspect
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('h5file',                help='experiment file (something.h5)')
    parser.add_argument('-c', '--cond_name', nargs='*',   help='conditions within experiment')
    parser.add_argument('-r', '--roi_name',  nargs='*',   help='names of rois to export')
    parser.add_argument('out_dir',   nargs='?', help='output directory')
    parser.add_argument('-l', '--conds',          help='list conditions', action='store_true')
    parser.add_argument('-s', '--subs',          help='list subjects', action='store_true')
    args = parser.parse_args()

    if args.conds or args.subs:
        summarize(args.h5file, args.conds, args.subs)

    sigargs = inspect.getcallargs(export_roi, 'DNE')
    kwargs = {k:v for k,v in args.__dict__.iteritems() if v != 'DNE'}
    export_roi(**kwargs)
