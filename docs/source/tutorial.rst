Tutorial
========
This tutorial will take you through how to setup and run a full pipeline.

Setting up Pipeline Directory
-----------------------------
Make a new directory for your pipeline, then run ``python $PATH_TO_PYCORR/setup_example.py full``.
This will install the following in the current directory:

1. config.yaml
2. run.py
3. pycorr_full_data

Example Data Structure
----------------------
The example data is contained in ``pycorr_full_data``. 
Looking at how this directory is set up will give a sense for how the config file works.
The directory should look like this::

    ├── rois
    │   ├── BHipp_3mm_thr30.nii.gz
    │   ├── ...
    │   └── precun_3mm_thr50.nii.gz
    └── subjects
        ├── OG_062312
        │   ├── intact
        │   │   └── trans_filtered_func_data.nii.gz
        │   ├── mini_wordscram
        │   │   └── trans_filtered_func_data.nii.gz
        │   ├── short_wordscram
        │   │   └── trans_filtered_func_data.nii.gz
        │   └── wordscram
        │       └── trans_filtered_func_data.nii.gz
        ├── ...
        │   ├── ...
        │   │   └── ...
        ├── aud_env_wordscram.npy
        ├── pieman_intact_audenv.mat
        └── pieman_wordscram_audenv.mat

With a little bit of squinting, it's clear the structure is::

    ├── rois
    │   └── {roi_id}
    └── subjects
        ├── {sub_id}
        │   ├── runA
        │   │   └── trans_filtered_func_data.nii.gz
        │   └── runB
        │       └── trans_filtered_func_data.nii.gz
        │
        └── OTHER THINGS

This is very config-able!

Config file
-----------
The configuration file uses `Yaml`_, and consists of experiment-wide and condition specific settings.
Here, necessary parameters are in **bold**.

Experiment-wide parameters are:
    * **roi_files**: path to rois to-be-loaded into experiment
    * **sub_folder**: where to store data from each Run. A single hdf5 file will be created per participant.

Condition-wide parameters are:
    * base_dir: used as prefix for all paths in condition
    * **nii_files**: path to niftis (or npy files) to-be-loaded
    * audio_env: path to audio envelope for stimulus
    * offset: default 0. Shift each timecourse forward by offset.
    * max_len: default None (all). Maximum length of timecourse.
    * threshold: default ???. Drop voxels with mean activation below threshold.
    * prop_pass_thresh: default ???. How many participants need to be above threshold to keep voxel at group level. 

**Using ``{var}`` magic**:
Setting ``roi_files`` and ``nii_files`` can be done by inserting ``{roi_file}`` or ``{nii_file}`` into their respective paths.
For example, in the data we have in our folder, the path to rois is ``rois/{roi_id}.nii*``.
In searching for the path, ``pycorr`` will use both ``{roi_id}`` and ``*`` as wildcards.
However, it will enter each matched nifti with the name of whatever matches ``{roi_id}``.

This takes advantage of :py:func:`pycorr.subject.get_subject_files`.

**Repeating default condition parameters**.
Yaml lets you set parameters from one condition as defaults for another.
An example of this is in the ``config.yaml``.
The following snippet shows a default condition in action::

    default: &default
        max_offset: 10

    conditions:
        condA:
            <<: *default

Running Pipeline
----------------
To start the pipeline, enter ``python run.py`` in the console.
The starter script, ``run.py``, iterates over conditions and performs
alignment, thresholding, makes a composite, and performs inter-subject correlation.
 

.. _Yaml: http://www.yaml.org/start.html

