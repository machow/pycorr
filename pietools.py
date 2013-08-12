from os import path, walk
import os
from glob import glob
import argparse
import fnmatch
try: import dicom
except: "lib 'dicom' unavailable"
from collections import Counter, OrderedDict

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
