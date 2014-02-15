import h5py
import nibabel as nib
from funcs_correlate import shift
from pietools import load_nii_or_npy

class Run:
    """Class to wrap /Subjects/sub_id/Cond"""

    def __init__(self, h5grp):
        """Load individual scan, given Group object.

        """
        self.grp = h5grp
        self.data    = h5grp['data']        if 'data'       in h5grp else None
        self.thresh  = h5grp['thresh_mask'] if 'thresh_mask'in h5grp else None
        
    def load(self, block=None):
        """
        return data as a numpy array. Select from block, otherwise use offset attr to shift.
        TODO: use slices
        """
        if block: BaseException

        return shift(self.data, h=self.data.attrs['offset'], outlen=self.data.attrs['max_len'])   #shift last dim forward h

    def threshold(self, threshold, data=None, save=False):
        """Boolean mask of values below threshold"""
        if not data:
            data = self.load()
        thresh_mask = data.mean(axis=-1) < threshold
        
        if save:
            if not self.thresh: 
                self.thresh = self.grp.create_dataset('thresh', shape=self.data.shape[:-1], dtype=bool)
            assert self.thresh.shape == thresh_mask.shape           #TODO create thresh if not exist
            self.thresh[...] = thresh_mask
            self.thresh.attrs['threshold'] = threshold

        return thresh_mask

    def create_dataset(self, data, overwrite=False, ref=False, **kwargs):
        """Data can be np.ndarray, nifti, or file name"""
        #TODO should fill attributes be mandatory?
        #TODO ref and overwrite args
        if self.data and not overwrite:
            raise BaseException('data already exists')
        else:
            if type(data) is str: np_arr = load_nii_or_npy(data)                #string
            elif type(data) is nib.nifti1.Nifti1Image: np_arr = data.get_data() #nifti
            else: np_arr = data                                                 #np array
            self.data = self.grp.create_dataset('data', data=np_arr)                        #TODO chunk shapes, what are defaults?
            if kwargs: self.fill_attributes(**kwargs)

    def fill_attributes(self, offset=0, max_len=None, exclude=False, notes="", **kwargs): #TODO should initial arguments be set from Exp setup, so they can be in the yaml?
        """Kwargs are unused (there to soak up uneccesary args)"""
        if not self.data: raise BaseException('data does not exist')
        if not max_len: max_len = self.data.shape[-1]
        self.data.attrs['offset'] = offset
        self.data.attrs['max_len'] = max_len
        self.data.attrs['exclude'] = exclude
        self.data.attrs['notes'] = notes
        #'blocks': pandas,
        #nii_hdr,
        #'date_scanned': unix date?

    def summary(self):
        print 'Data:\t',self.data
        print 'thresh:\t', self.thresh
        width = max(map(len, self.data.attrs))
        for key, val in self.data.attrs.iteritems(): 
            print "{:{}}".format(key, width), ':    ', val


import os
from pietools import query_to_re_groups
from glob import glob
import re

class Exp:
    """

    """

    def __init__(self, f):
        self.f = h5py.File(f) if type(f) is str else f

    def setup(self, config):

        self.f.create_group('subjects')
        self.f.create_group('rois')
        self.f.create_group('conds')
        for condname, cond in config['conds'].iteritems():
            g_cond = self.create_cond(condname, **cond)

            full_query = os.path.join(config['base_dir'], cond['nii_files'])
            #Create subject data
            for m in self.get_subject_files(full_query):
                fname = m.group()
                sub_id = m.groupdict()['sub_id']
                print fname, sub_id
                self.create_subrun(sub_id, condname, fname, **cond)
                #have option to do threshold?
        self.f.flush()

    def create_cond(self, condname, offset=None, max_len=None, **kwargs):
        cond = self.f['conds'].create_group(condname)
        cond.attrs['offset'] = offset
        cond.attrs['max_len'] = max_len
        cond.create_group('blocks')
        #TODO blocks (pandas)
        return cond
    
    @staticmethod
    def get_subject_files(globpath):
        """Iterator. Returns re match object for each matching path.
        
        Full path accessed using m.group(); args from m.groupdict()
        """
        to_glob = globpath.format(sub_id='*')
        to_re = query_to_re_groups(globpath)   #for getting sub_id out of results
        nii_files = glob(to_glob)
        
        for fname in nii_files:
            yield re.match(to_re, fname)

    def create_subrun(self, sub_id, condname, fname_or_arr, ref=False, **kwargs):
        g_sub = self.f['subjects'].create_group('%s/%s'%(sub_id, condname))
        run = Run(g_sub)
        run.create_dataset(data=fname_or_arr, ref=False)
        run.fill_attributes(**kwargs)

        return run

    def get_cond(self, condname):
        return self.f['conds/%s'%condname]

    def iter_runs(self, condname):
        for sname, sub in self.f['subjects'].iteritems():
            for cname, cond in sub.iteritems():
                if cname == condname: yield Run(cond)


    #TODO

    def create_roi(self, roiname):
        roi = self.f['rois'].create_group(roiname)
        return roi

