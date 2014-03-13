from subject import Exp
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('h5file',      help='experiment file (something.h5)')
parser.add_argument('cond_name',   help='condition within experiment')
parser.add_argument('roi_name',    help='name of roi to export')
args = parser.parse_args()

E = Exp(args.h5file)
roi = E.f['rois'][args.roi_name][...]
data = pd.DataFrame({run.grp.parent.name.split('/')[1] :
                     run.load(roi=roi, threshold=True, standardized=True) for run in E.iter_runs(args.cond_name)})

data.to_csv(args.cond_name + '_' + args.roi_name + '.csv')
