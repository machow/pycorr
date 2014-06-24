.. pycorr documentation master file, created by
   sphinx-quickstart on Wed Jun 18 22:54:39 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pycorr's documentation!
==================================
Pycorr allows you to produce quick, replicable pipelines for inter-subject correlation analyses.
Oftentimes, analysis scripts become a tangled mess of experiment specific parameters and file paths.
Sometimes you're left with a dozen different folders with the data saved in slightly different ways, and names like "sub1_cropleft2_cropright4.nii".
We've all been there. No judgement. pycorr can help alleviate your research woes by doing a few things:

1. Seperating parameter configuration from analysis scripts
2. Modular structure, with functions for loading, analyzing, and saving data
3. Support for HDF5 using `H5py`_. This allows you to store parameters in the data, rather than in it's filename :).
4. Command-line tools for running analyses
5. Unit-tests and documentation

_`H5py`: http://www.h5py.org/

General Workflow
================
The bulk of most pycorr analyses consist of two parts:

1. Configuration file
2. Analysis script

The configuration file has basic settings, such as filepaths to the data.
The analysis scripts takes advantage of the two classes, ``Exp`` and ``Run``.
``Exp`` loads the data using the config file, and ``Run`` wraps each functional run.

Contents:

.. toctree::
   :maxdepth: 2

   setup
   tutorial
   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

