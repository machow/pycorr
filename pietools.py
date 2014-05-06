from os import path, walk
import os
from glob import glob
import fnmatch
try: import dicom
except: "lib 'dicom' unavailable"
from collections import OrderedDict
import re
import nibabel as nib
import numpy as np

def load_nii_or_npy(fname):
    """convenience function to load nib or np filetypes"""
    if os.path.splitext(fname)[1] == '.npy': return np.load(fname)
    else: return nib.load(fname).get_data()

def copy_nii_hdr(nii, data, save=""):
    """Returns a nifti object using another nifti's header
    
    Dimensions are changed to correspond with data shape

    parameters:
        nii - nifti file to take header from
        data - data file to convert to Nifti1Image
        save - file name to save to, no save if blank
    """
    if type(nii) is str: nii = nib.load(nii)
    hdr = nii.get_header()
    affine = nii.get_affine()
    new_hdr = hdr.copy()
    new_hdr.set_data_shape(data.shape)
    new_nii = nib.Nifti1Image(data, affine, new_hdr)
    if save: nib.save(new_nii, save)
    return new_nii


def query_to_re_groups(query):
    """For a string with with glob syntax, allow for named groups for wildcards
    
    E.g. "*/{sub_id}/file.txt" is converted to a regular expression where {sub_id} is a wildcard,
         and it's match is a group named sub_id.
        
    Supports multiple {bracketed_args}
    """
    var_matches = re.finditer('{(?P<var>.*?)}', query)
    glob_exp = query
    re_exp = fnmatch.translate(query)
    for m in var_matches: 
        glob_exp = glob_exp.format(**{m.groupdict()['var']: '*'})
        re_group_exp = '(?P<%s>.*?)'%m.groupdict()['var']
        re_exp = re_exp.replace(re.escape(m.group()), re_group_exp)
    return glob_exp, re_exp


def subs_getfile(dirname, matchme, verbose = False):
    """Recursively find all matching files, return as dictionary with entry for each subdir

    Parameters:
    dirname -- directory to search
    matchme -- string to match.  May use wildcards (glob used)

    """
    matches = {}
    subs = os.listdir(dirname)
    for sub in subs:
        matches[sub] = searchr(path.join(dirname, sub), matchme, verbose)
    return matches

def splice_dir(dirname, save=False):
    data = np.vstack([np.load(os.path.join(dirname, fname)) for fname in os.listdir(dirname)])
    if save: np.save(dirname + '.npy', data)
    return data

import errno

def mkdir_p(path):
        try:
            os.makedirs(path)
        except OSError as exc: # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else: raise

# Script functions ------------------------------------------------------------
def parse_section(snapshot):
    if snapshot == 'all':
        return slice(None)
    try:
        section = int(snapshot)
    except ValueError:
        section = [int(s) if s else None for s in snapshot.split(':')]
        if len(section) > 3:
            raise ValueError('snapshots input incorrect')
        section = slice(*section)
    return section

import gc
def arr_slice(nii_file, _slice):
    "slices indx from first dim from matrix, for parallelization"
    print 'loading data'
    if type(nii_file) is str:
        nii = nib.load(nii_file)
        dat = nii.get_data()
    else: dat = nii_file
    out = dat[_slice].copy()
    del dat, nii                   #garbage collection
    gc.collect()
    return out

######################

def searchr(dirname, matchme, verbose = True):
    """Walk dir, return list of all matches.  Use glob to match."""

    matches = []
    for root, subfold, fnames in walk(dirname):
        for f in glob(path.join(root, matchme)):
            matches.append(f)
            if verbose: print """{subfold}: \t\t{files}""".format(subfold = root, files = path.basename(f))
    return matches



def print_Blocks(dirname):
    """Given subject dir, print number of TRs for each block"""

    data = [fname.split('_')[-2:] for fname in os.listdir(dirname)]
    if len([entry for entry in data if len(entry) == 1]) > 0:
        print "Can't parse some file names.  Are there other files in the directory?  Correcting.."
        data = [entry for entry in data if len(entry) == 2]
    blocks = set([entry[0] for entry in data])
    for b in sorted(blocks, key=lambda x: x[-3:]):
        print 'Block %s'%b, '\t', max([entry[1] for entry in data if entry[0] == b])


def loadwith(subdir, path, func, **kwargs):
    """For all subs in basedir, apply func to path, return list.

    Parameteres:
    subdir --	base directory to iterate through. If none, just glob path.
    path --	path to search for in each dir in subdir
    func --	function to apply to each path found

    """
    out = OrderedDict()
    if subdir:
        for sub in os.walk(subdir).next()[1]:
            fname = os.path.join(subdir, sub, path)
            out[sub] = func(fname, **kwargs)
    else:
        fnames = glob(path)
        for name in fnames:
            subj = os.path.split(name)[-1][:9]
            out[subj] = func(name, **kwargs)
    return out
#########
###
###

def print_dicom(dirname):
    lookup = ['SeriesNumber', 'SeriesDescription']
    entries = []
    tallies = {}
    count = 1
    fnames = os.listdir(dirname)
    print fnames[:50]
    sorted(fnames, reverse = True)

    for fname in fnames:
        data = dicom.ReadFile(os.path.join(dirname, fname), force = True)
        hdr = [str(getattr(data, entry, "NO %s"%entry)) for entry in lookup]
        if hdr in entries:
            count +=1
        else:
            entries.append(hdr)
            tallies[str(hdr)] = str(count)
            print '\t'.join(hdr + [str(count)])
            count = 1
    return tallies

def test(dirname):
    print [ii for ii in os.walk(dirname)][:3]
