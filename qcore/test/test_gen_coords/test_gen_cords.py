from qcore import shared
import pytest
import os
""" Command to run this test: 'python -m pytest -v -s test_gen_cords.py' 
If the test passed it will delete the files in the output folder. 
Otherwise it would not delete files """
""" Instructions: Sample1 folder contains a sample output taken from hypocentre. Its path is noted in the readme file. In that path you will find the 
 params_vel.py along with other 5 output files. Use them as the benchmark files.If you want another sample to be tested, 
 create a similar folder like sample1 and store the relevant files there (e.g:sample2). While running the test change sample1 to sample2"""

PATH_TO_SAMPLE_DIR = os.path.join(os.getcwd(),"sample1")
PATH_TO_SAMPLE_OUTDIR = os.path.join(PATH_TO_SAMPLE_DIR, "output")
# print "PATH_TO_SAMPLE_OUTDIR: ",PATH_TO_SAMPLE_OUTDIR
PATH_TO_NEW_OUTDIR = os.path.join(os.getcwd(), "output")
# print "PATH_TO_NEW_OUTDIR: ", PATH_TO_NEW_OUTDIR
PATH_FOR_PRG_TOBE_TESTED = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), "gen_coords.py"))
# print "PATH_FOR_PRG_TOBE_TESTED: ",PATH_FOR_PRG_TOBE_TESTED

def setup_module():
    """ create a symbolic link for params_vel.py"""
    print "---------setup_module------------"
    sample_path = os.path.join(PATH_TO_SAMPLE_DIR, "params_vel.py")
    os.symlink(sample_path,os.path.join(os.getcwd(), "params_vel.py"))

def test_gencords():
    """ test qcore/gen_coords.py """
    print "---------test_gencords------------"
    shared.exe("python "+PATH_FOR_PRG_TOBE_TESTED+" "+PATH_TO_NEW_OUTDIR)
    out,err=shared.exe("diff -qr "+PATH_TO_NEW_OUTDIR+" "+PATH_TO_SAMPLE_OUTDIR)
    assert out == "" and err == ""
    remove_files()

def remove_files():
    """ delete the output files generated if test_gencords passed"""
    print "---------remove_files------------"
    filelist = [f for f in os.listdir(PATH_TO_NEW_OUTDIR) ]
    for f in filelist:
        os.remove(os.path.join(PATH_TO_NEW_OUTDIR, f))

def teardown_module():
    """ delete the symbolic link for params_vel.py"""
    print "---------teardown_module------------"
    if (os.path.isfile(os.path.join(os.getcwd(), "params_vel.py"))):
        os.remove(os.path.join(os.getcwd(), "params_vel.py"))