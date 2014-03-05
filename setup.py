import os, shutil, sys

basedir = os.path.dirname(os.path.realpath(__file__))
shutil.copy2(os.path.join(basedir, 'config.yaml'), 'config.yaml')


if len(sys.argv) > 1 and sys.argv[1] == 'test':
    shutil.copytree(os.path.join(basedir, 'tests/data/fullpipe'), 'pipeline')
    shutil.copy2(os.path.join(basedir, 'tests/test_subject.py'), 'run.py')
    os.symlink(os.path.abspath('../pieman'), 'pieman')
