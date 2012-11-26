#   Create a matrix which holds the correlation matrices for each voxel
#
#   M Chow, Oct 2012
#   machow@princeton.edu
#

from scipy import io as sio
from funcs_correlate import loadData, crosscor

#load data
fname = 'trans_filtered_func_data_sync.npy'
pipeargs = dict(subdir = 'pieNDiv/subjects',
                mrpaths = ['analysis/preproc/preproc01.feat/'+fname],
                condnames = ['intact01'],
                behpaths = [])

olgaargs = dict(subdir = 'pieNDiv/subsnotpipe',
                mrpaths = ['intact1.feat/'+fname],
                condnames = ['intact01'],
                behpaths = [])

subdata = loadData(**pipeargs) + loadData(**olgaargs)       #list of subject dicts

#Get a dictionary with sub : nifti for intact01 condition
ddict = {sub['Sub'] : sub['Ni'].values()[0] for sub in subdata}
sub_indx = sorted(ddict.keys())
dlist = [ddict[sub] for sub in sub_indx]

#Correlation matrix
#This will likely run on the cluster without taking up too much memory, but can be segmented spatially if needed
C = crosscor(dlist)

#Intersubject Correlations (mean off-diagonal crosscorrelation)
N = C.shape[-1]
mean_C = (C.sum(-1) - 1) / (N - 1)      #Adjust for 1 on diagonal

sio.savemat('pieNDiv/intersubj/ccr_intact01.mat', dict(CCR = C, sub_indx=sub_indx))
sio.savemat('pieNDiv/intersubj/ccr_mean_intact01', dict(CCR_mean = mean_C, sub_indx=sub_indx))
