import h5py
import nibabel as nib
from funcs_correlate import shift, standardize
from pietools import load_nii_or_npy
import numpy as np
from scipy.stats import nanmean


class Run:
    """Class to wrap /Subjects/sub_id/Cond"""

    def __init__(self, h5grp):
        """Load individual scan, given Group object.

        """
        self.grp = h5grp
        self.data    = h5grp['data']    if 'data'   in h5grp else None
        self.thresh  = h5grp['thresh']  if 'thresh' in h5grp else None
        self.attrs = h5grp.attrs
        
    def load(self, fourD=False, block=None, standardized=False, threshold=False, roi=False):
        """
        return data as a numpy array. Select from block, otherwise use offset attr to shift.
        TODO: use slices, change roi so that it properly takes mean along xyz axes (instead of assuming flat)
        Roi forces standardize
        """
        if block: BaseException
        #ROI (could also set roi equal to ellipsis or have shift return a slice
        if np.any(roi):
            if threshold: roi = roi & ~self.thresh[...]
            M_roi = self.data[...][roi]    #Have to load into mem :(
            shifted = shift(M_roi, h=self.data.attrs['offset'], outlen=self.data.attrs['max_len'])   #shift last dim forward h
            if standardized: standardize(shifted, inplace=True)
            return standardize(shifted.mean(axis=0))
        else: 
            shifted = shift(self.data, h=self.data.attrs['offset'], outlen=self.data.attrs['max_len'])   #shift last dim forward h

            #other args
            if standardized: standardize(shifted, inplace=True)
            if threshold: shifted[self.thresh[...]] = np.nan
            if fourD: shifted = shifted.reshape(self.data.attrs['shape'])
        
            return shifted

    def threshold(self, threshold, data=None, save=False):
        """Boolean mask of values below threshold.
        
        May fail if data contains nans
        """
        data = data if np.any(data) else self.load()
        thresh_mask = data.mean(axis=-1) < threshold
        
        if save:
            self.thresh = self.grp.require_dataset('thresh', shape=self.data.shape[:-1], dtype=bool)
            assert self.thresh.shape == thresh_mask.shape           #TODO create thresh if not exist
            self.thresh[...] = thresh_mask
            self.thresh.attrs['threshold'] = threshold
        self.grp.file.flush()
        return thresh_mask

    def create_dataset(self, data, overwrite=False, ref=False, compression='gzip', chunks=(10,10,10), **kwargs):
        """Data can be np.ndarray, nifti, or file name"""
        #TODO should fill attributes be mandatory?
        #TODO ref and overwrite args
        if self.data and not overwrite:
            raise BaseException('data already exists')

        if type(data) is str: np_arr = load_nii_or_npy(data)                #string
        elif type(data) is nib.nifti1.Nifti1Image: np_arr = data.get_data() #nifti
        else: np_arr = data                                                 #np array
        #rdy_arr = np_arr.reshape([-1, np_arr.shape[-1]])
        chunks += np_arr.shape[-1:]
        self.data = self.grp.create_dataset('data', data=np_arr, chunks=chunks, compression=compression)
        if kwargs: self.fill_attributes(**kwargs)
        self.grp.file.flush()
        #self.data.attrs['shape'] = np_arr.shape

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
        self.grp.file.flush()

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
from funcs_correlate import sum_tc

class Exp:
    """

    """

    def __init__(self, f):
        self.f = h5py.File(f) if type(f) is str else f

    def setup(self, config, create_conds=False):
        #if 'data_storage' in config:
        #    for k,v in config['data_storage'].iteritems():
        #        self.f.attrs[k] = v

        self.f.create_group('subjects')
        self.f.create_group('rois')
        self.f.create_group('conds')
        self.f.attrs['sub_folder'] = config['sub_folder']
        try: os.mkdir(config['sub_folder'])
        except: pass #TODO add exception expected

        #Load roi files in from config
        if 'roi_files' in config:
            for m_roi in self.get_subject_files(config['roi_files']):
                fname = m_roi.group()
                roi_id= m_roi.groupdict()['roi_id']
                self.create_roi(roi_id, fname)
        #create new conditions
        if create_conds:
            for condname, cond in config['conds'].iteritems():
                self.create_cond(condname, **cond)
        self.f.flush()

    def create_cond(self, condname, offset=0, max_len=None, threshold=0, audio_env=None, 
                    base_dir="", nii_files=None, **kwargs):
        cond = self.f['conds'].create_group(condname)
        cond.attrs['offset'] = offset
        cond.attrs['max_len'] = max_len
        cond.attrs['threshold'] = threshold
        cond.create_group('blocks')
        cond.create_group('correlations')
        cond.create_group('analyses')
        
        #Load conditions in from config
        
        if audio_env:
            aud_path = os.path.join(base_dir, audio_env)
            print aud_path
            cond.create_dataset('audio_env', data=load_nii_or_npy(aud_path))

        full_query = os.path.join(base_dir, nii_files)
        #Create subject data
        for m in self.get_subject_files(full_query):
            fname = m.group()
            sub_id = m.groupdict()['sub_id']
            print fname, sub_id
            self.create_subrun(sub_id, condname, fname, **cond)
            #have option to do threshold?
            self.f.flush()
        #TODO blocks (pandas)
        return cond
    
    def create_subrun(self, sub_id, condname, fname_or_arr, ref=False, **kwargs):
        path = '%s/%s'%(sub_id, condname)
        #Remote link if sub_folder is specified
        if self.f.attrs['sub_folder'] and not sub_id in self.f['subjects']:
            fname = os.path.join(self.f.attrs['sub_folder'], sub_id+'.h5')
            with h5py.File(fname) as f_new:
                f_new.create_group(sub_id)
            self.f['subjects/%s'%sub_id] = h5py.ExternalLink(fname, sub_id)

        #add condition to subject group
        g_sub = self.f['subjects'].create_group('%s/%s'%(sub_id, condname))
        run = Run(g_sub)
        run.create_dataset(data=fname_or_arr, ref=False)
        run.fill_attributes(**kwargs)

        return run

    def get_cond(self, condname):
        #TODO change to getitem method?
        return self.f['conds/%s'%condname]

    def iter_runs(self, condname):
        for sname, sub in self.f['subjects'].iteritems():
            for cname, cond in sub.iteritems():
                if cname == condname: yield Run(cond)

    def N_runs(self, condname):
        return len(list(self.iter_runs(condname)))

    def gen_composite(self, condname, overwrite=False):
        #TODO overwrite
        data = (run.load(standardized=True, threshold=True) for run in self.iter_runs(condname))
        shape = self.iter_runs(condname).next().load().shape
        composite = sum_tc(data, shape=shape, standardize_out=True)
        self.get_cond(condname).create_dataset('composite', data=composite)

    def gen_cond_thresh(self, condname, overwrite=False):
        cond = self.get_cond(condname)

        dlist = (~run.thresh[...] for run in self.iter_runs(condname))              #all that pass threshold
        mustpassprop = cond.attrs['prop_pass_thresh']
        n_must_pass = (mustpassprop) * len(dlist)
        above_thresh = np.sum(dlist, axis=0)
        thresh_fail = above_thresh < n_must_pass 			#sum num of failed subjects per voxel

        if 'threshold' not in cond:
            cond.create_dataset('threshold', data=thresh_fail, dtype=bool)
        else: 
            cond['threshold'][...] = thresh_fail

    def summarize(self, condname):
        cond = self.get_cond(condname)
        def printvalue(name):
            print name, ':\t', cond[name]
        cond.visit(printvalue)

    @staticmethod
    def get_subject_files(globpath):
        """Iterator. Returns re match object for each matching path.
        
        Full path accessed using m.group(); args from m.groupdict()
        """
        #to_glob = globpath.format(sub_id='*')
        to_glob, to_re = query_to_re_groups(globpath)   #for getting sub_id out of results
        nii_files = glob(to_glob)
        
        for fname in nii_files:
            yield re.match(to_re, fname)

    def create_roi(self, roiname, fname):
        roi_data = load_nii_or_npy(fname)
        roi = self.f['rois'].create_dataset(roiname, data=roi_data.astype(bool))
        return roi

class Cond:
    """
    """
    def __init__(self):
        pass


