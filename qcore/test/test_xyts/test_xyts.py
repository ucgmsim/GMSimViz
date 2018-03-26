import os
from qcore import xyts
import pytest
from qcore.test.tool import utils
import numpy as np
import getpass
from datetime import datetime
import sys
import shutil
import errno

XYTS_ACTUAL_PATH ="/nesi/projects/nesi00213/sample_data/test_qcore/xyts.e3d"
SAMPLE_OUT_DIR_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)),"sample1/output")
XYTS_FILE = os.path.join(os.path.abspath(os.path.dirname(__file__)),"sample1/input/xyts.e3d")
os.symlink(XYTS_ACTUAL_PATH, XYTS_FILE)
OBJ_XYTS = xyts.XYTSFile(XYTS_FILE)
SAMPLE_PGV = os.path.join(SAMPLE_OUT_DIR_PATH,"sample_pgvout")
SAMPLE_MMI = os.path.join(SAMPLE_OUT_DIR_PATH,"sample_mmiout")
TMP_DIR_NAME = (os.path.join("/home/",getpass.getuser(),("tmp_" + os.path.basename(__file__)[:-3] + '_' + ''.join(str(datetime.now()).split())).replace('.', '_')).replace(
            ':', '_'))

def setup_module(scope="module"):
    """ create a tmp directory for storing output from test"""
    print "----------setup_module----------"
    try:
        os.mkdir(TMP_DIR_NAME)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def teardown_module():
    """ delete the symbolic link 
    delete the tmp directory if it is empty"""
    print "---------teardown_module------------"
    if (os.path.isfile(XYTS_FILE)):
        os.remove(XYTS_FILE)
    if (len(os.listdir(TMP_DIR_NAME)) == 0):
        try:
            shutil.rmtree(TMP_DIR_NAME)
        except (IOError, OSError) as (e):
            sys.exit(e)


@pytest.mark.parametrize("gmt_format, expected_corners",[(False,[[170.91449, -43.2916], [172.26257, -44.21773], [171.01305, -45.15298], [169.66354, -44.21213]]) ,\
(True,([[170.91449, -43.2916], [172.26257, -44.21773], [171.01305, -45.15298], [169.66354, -44.21213]], '170.91449 -43.2916\n172.26257 -44.21773\n171.01305 -45.15298\n169.66354 -44.21213')
)])
def test_corners(gmt_format, expected_corners):
   assert OBJ_XYTS.corners(gmt_format=gmt_format)  == expected_corners


@pytest.mark.parametrize("corners, expected_region",[(None,(169.66354000000001, 172.26257000000001, -45.152979999999999, -43.291600000000003)) ,\
([[170.91449, -43.2916], [172.26257, -44.21773], [171.01305, -45.15298], [169.66354, -44.21213]],(169.66354000000001, 172.26257000000001, -45.152979999999999, -43.291600000000003))])
def test_region(corners, expected_region):
    assert OBJ_XYTS.region(corners) == expected_region


@pytest.mark.parametrize("mmi, pgvout, mmiout, sample_pgv, sample_mmi",[(True, None, None, SAMPLE_PGV, SAMPLE_MMI),(False, None, None, SAMPLE_PGV, SAMPLE_MMI), \
                                                                        (True, "test_pgv_path1", "test_mmi_path1", SAMPLE_PGV, SAMPLE_MMI),(False, "test_pgv_path2", None, SAMPLE_PGV, SAMPLE_MMI)])
def test_pgv(mmi, pgvout, mmiout, sample_pgv, sample_mmi):
    files_to_del = []
    if pgvout:
        pgvout = os.path.join(TMP_DIR_NAME,pgvout)
        files_to_del.append(pgvout)
    if mmiout:
        mmiout = os.path.join(TMP_DIR_NAME,mmiout)
        files_to_del.append(mmiout)

    xyts_test_output_array = OBJ_XYTS.pgv(mmi=mmi, pgvout=pgvout, mmiout=mmiout)

    if pgvout:
        sample_pgv_array = np.fromfile(sample_pgv, dtype='3<f4')
        test_pgvout_array = np.fromfile(pgvout, dtype='3<f4')
        utils.compare_np_array(sample_pgv_array,test_pgvout_array)
        if mmiout:
            sample_mmi_array = np.fromfile(sample_mmi, dtype='3<f4')
            test_mmiout_array = np.fromfile(mmiout, dtype='3<f4')
            utils.compare_np_array(sample_mmi_array, test_mmiout_array)

    else:
        if not mmi:
            sample_pgv_array = np.fromfile(sample_pgv, dtype='3<f4')
            utils.compare_np_array(sample_pgv_array, xyts_test_output_array)
        elif mmiout == None:
            pgv,mmi = xyts_test_output_array
            sample_pgv_array = np.fromfile(sample_pgv, dtype='3<f4')
            utils.compare_np_array(sample_pgv_array, pgv)
            sample_mmi_array = np.fromfile(sample_mmi, dtype='3<f4')
            utils.compare_np_array(sample_mmi_array, mmi)


    for f in files_to_del:
        utils.remove_file(f)


@pytest.mark.parametrize("step, comp, test_outfile, sample_outfile",[(10, -1, None, "out_tslice-1"),\
                        (10, -1,"test_tslice-1", "out_tslice-1"),(10, 0, "test_tslice0", "out_tslice2"), \
                        (10, 1, "test_tslice1", "out_tslice1"),(10, 2, "test_tslice2", "out_tslice0")])
def test_tslice_get(step, comp, test_outfile, sample_outfile):
    files_to_del = []
    if test_outfile:
        test_outfile = os.path.join(TMP_DIR_NAME, test_outfile)
        sample_outfile = os.path.join(SAMPLE_OUT_DIR_PATH, sample_outfile)
        result = OBJ_XYTS.tslice_get(step, comp, test_outfile)
        if result:
            sample_array = np.fromfile(sample_outfile, dtype='3<f4')
            utils.compare_np_array(sample_array, result)
        else:
            utils.compare_binary_files(sample_outfile, test_outfile)
            files_to_del.append(test_outfile)
    for f in files_to_del:
        utils.remove_file(f)





