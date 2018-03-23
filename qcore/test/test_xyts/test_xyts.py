import os
from qcore import xyts
import pytest

sample_out_dir_path = "/nesi/projects/nesi00213/sample_data/test_qcore"
xyts_file = os.path.join(sample_out_dir_path,"xyts.e3d")
obj_xyts = xyts.XYTSFile(xyts_file)


@pytest.mark.parametrize("gmt_format, expected_corners",[(False,[[170.91449, -43.2916], [172.26257, -44.21773], [171.01305, -45.15298], [169.66354, -44.21213]]) ,\
(True,([[170.91449, -43.2916], [172.26257, -44.21773], [171.01305, -45.15298], [169.66354, -44.21213]], '170.91449 -43.2916\n172.26257 -44.21773\n171.01305 -45.15298\n169.66354 -44.21213')
)])
def test_corners(gmt_format, expected_corners):
   assert obj_xyts.corners(gmt_format=gmt_format)  == expected_corners


@pytest.mark.parametrize("corners, expected_region",[(None,(169.66354000000001, 172.26257000000001, -45.152979999999999, -43.291600000000003)) ,\
([[170.91449, -43.2916], [172.26257, -44.21773], [171.01305, -45.15298], [169.66354, -44.21213]],(169.66354000000001, 172.26257000000001, -45.152979999999999, -43.291600000000003))])
def test_region(corners, expected_region):
    assert obj_xyts.region(corners) == expected_region


@pytest.mark.parametrize("mmi, pgvout, mmiout",[(True, None, None),(False, None, None),(True, "pgvout_path1", "mmiout_path1"),(False, "pgvout_path2", None)])
def test_pgv(mmi, pgvout, mmiout):
    if pgvout:
        pgvout = os.path.join(sample_out_dir_path,pgvout)
    if mmiout:
        mmiout = os.path.join(sample_out_dir_path,mmiout)
    obj_xyts.pgv(mmi=mmi, pgvout=pgvout, mmiout=mmiout)

    # assert outpgv ==


