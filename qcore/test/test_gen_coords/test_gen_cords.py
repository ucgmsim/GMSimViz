""" Command to run this test: 'python -m pytest -v -s test_gen_cords.py' 
If the test passed it will delete the files in the output folder. 
Otherwise it would not delete files 
 
 Instructions: Sample1 folder contains a sample output taken from hypocentre. Its path is noted in the readme file. In that path you will find the 
 params_vel.py along with other 5 output files. Use them as the benchmark files.If you want another sample to be tested, 
 create a similar folder like sample1 and store the relevant files there (e.g:sample2). While running the test change sample1 to sample2

Just to run : py.test (or) python -m pytest -s -v test_gen_cords.py
To know the code coverage : py.test --cov=test_gen_cords.py
To know the test coverage :python -m pytest --cov ../../gen_cords.py test_gen_cords.py
"""

from qcore import shared
import os
import shutil
import getpass
from datetime import datetime
import sys

PATH_TO_SAMPLE_DIR = os.path.join(os.getcwd(),"sample1")
PATH_TO_SAMPLE_OUTDIR = os.path.join(PATH_TO_SAMPLE_DIR, "output")
# print "PATH_TO_SAMPLE_OUTDIR: ",PATH_TO_SAMPLE_OUTDIR

DIR_NAME = (os.path.join("/home/",getpass.getuser(),("tmp_" + os.path.basename(__file__)[:-3] + '_' + ''.join(str(datetime.now()).split())).replace('.', '_')).replace(
            ':', '_'))
print "PATH_TO_NEW_OUTDIR: ", DIR_NAME
PATH_FOR_PRG_TOBE_TESTED = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), "gen_cords.py"))
print "PATH_FOR_PRG_TOBE_TESTED: ",PATH_FOR_PRG_TOBE_TESTED

def setup_module(scope="module"):
    """ create a symbolic link for params_vel.py"""
    print "---------setup_module------------"
    try:
        os.mkdir(DIR_NAME)
    except Exception as (e):
        sys.exit(e)
    sample_path = os.path.join(PATH_TO_SAMPLE_DIR, "params_vel.py")
    os.symlink(sample_path,os.path.join(os.getcwd(), "params_vel.py"))

def test_gencords():
    """ test qcore/gen_coords.py """
    print "---------test_gencords------------"
    shared.exe("python "+PATH_FOR_PRG_TOBE_TESTED+" "+DIR_NAME)
    out,err=shared.exe("diff -qr "+DIR_NAME+" "+PATH_TO_SAMPLE_OUTDIR)
    assert out == "" and err == ""
    shutil.rmtree(DIR_NAME)

def teardown_module():
    """ delete the symbolic link for params_vel.py"""
    print "---------teardown_module------------"
    if (os.path.isfile(os.path.join(os.getcwd(), "params_vel.py"))):
        os.remove(os.path.join(os.getcwd(), "params_vel.py"))