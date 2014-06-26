Setup
=====

Installing Dependencies
-----------------------
All dependencies are listed in ``requirements.txt``.
It's a good idea to install everything into a `virtual environment`_, but feel free to be young and reckless :O.
Alternatively, the `Anaconda python distribution`_ contains all the necessary packages, except for `nibabel`_, which can be installed using ``pip install nibabel``.

Installing Pycorr
-----------------
cd to the directory with pycorr, then ``python setup.py install``.
Since the package is still somewhat of a moving target, you might want to symlink it by using ``python setup.py develop``.

Using pycorr on Rondo
--------------------------

pycorr is already installed in a virtual environment in the Hasson lab fileshare.
Activate it by editing your PATH in ``~/.local``::
  setenv PATH=/jukebox/hasson/michael/venv/anaconda:$PATH

Then, you can load the python environment by entering ``source activate trw-interl`` in the console.

If you have a different shell, such as bash, enabled by default--then it may be better to edit ``~/.bash_login``.
Replace ``setenv`` with ``export``.


.. _virtual environment: http://docs.python-guide.org/en/latest/dev/virtualenvs/
.. _Anaconda python distribution: https://store.continuum.io/cshop/anaconda/
.. _nibabel: http://nipy.org/nibabel/
