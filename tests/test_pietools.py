#import tempfile, os
#from helper import temp_chdir 
from pycorr.pietools import sort_int_folder

def test_sort_int_folder():
    files = ['1.npy', '10.npy', '2.npy']
    sort_int_folder(files)
    assert files == ['1.npy', '2.npy', '10.npy']
