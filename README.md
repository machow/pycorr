pycorr
======

Running Full Example 
--------------------

Create a new directory. Then, from this directory run `python setup.py test`.
This will copy the necessary data and scripts, including `run.py`. Enter 
`python run.py` to setup the hdf5 datastore and analyses.

TODO
----

Quick:

1. update of funcs_correlate
2. unit tests load_roi and roi_mask (should rethink these functions)
3. assign runs to group when creating condition (done?)

Extensive:

3. update Align (yield?), Shift (yield?), ISC (rewrite ROI script)
5. command-line interface

Down the Road:

1. ROI event plots
2. Class for unscrambling, subsetting?
3. Bootstrap
4. GUI
5. SEM
