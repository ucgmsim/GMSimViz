from qcore import srf, shared
import pytest

import os
""" Command to run this test: 'python -m pytest -v -s test_srf.py'  """

SRF_1_PATH = os.path.join(os.getcwd(), "sample1/Hossack_HYP01-01_S1244.srf")
SRF_2_PATH = os.path.join(os.getcwd(), "sample2/Tuakana13_HYP01-01_S1244.srf")
SRF_1_CNR_PATH = os.path.join(os.getcwd(), "sample1/output/cnrs.txt")
SRF_2_CNR_PATH = os.path.join(os.getcwd(), "sample2/output/cnrs.txt")
HEADERS = ['centre', 'nstrike', 'ndip', 'length', 'width', 'strike', 'dip', 'dtop', \
           'shyp', 'dhyp']
SRF_1_PLANES = srf.read_header(SRF_1_PATH, True)
SRF_2_PLANES = srf.read_header(SRF_2_PATH, True)





def test_nplane1():
    assert len(SRF_1_PLANES) == 1

def test_nplane2():
    assert len(SRF_2_PLANES) == 2

@pytest.mark.parametrize("plane, expected_values",[( SRF_1_PLANES[0], [[176.2354,-38.3404], 34, 92, 3.44, 9.24, 230, 60, 0.00, 0.00, 5.54]),
                                                (SRF_2_PLANES[0],[[176.8003, -37.0990], 46, 104, 4.57, 10.44, 21, 50, 0.00, 0.00, 6.27]),
                                                (SRF_2_PLANES[1], [[176.8263, -37.0622], 49, 104, 4.89, 10.44, 37, 50, 0.00, -999.90, -999.90])])
def test_plane(plane, expected_values):
    """ Tests for the header lines  """
    for i in xrange(len(HEADERS)):
        assert(plane[HEADERS[i]] == expected_values[i])

@pytest.mark.parametrize("test_dt, expected_dt",[(SRF_1_PATH, 2.50000e-02),
                                                (SRF_2_PATH,2.50000e-02),])
def test_dt(test_dt, expected_dt):
    assert srf.srf_dt(test_dt) == expected_dt

@pytest.mark.parametrize("test_dxy, expected_dxy",[(SRF_1_PATH, (0.10,0.10)),
                                                (SRF_2_PATH,(0.1,0.10)),])
def test_dxy(test_dxy, expected_dxy):
    assert srf.srf_dxy(test_dxy) == expected_dxy


@pytest.mark.parametrize("test_srf,seg,depth",[(SRF_1_PATH, -1, True),(SRF_2_PATH, -1, True), \
                                               (SRF_1_PATH, -1, False),(SRF_2_PATH, -1, False)])
def test_get_bounds(test_srf, seg, depth):
    print srf.get_bounds(test_srf, seg = seg, depth = depth)
    assert True

@pytest.mark.parametrize("test_srf,filename,sample_cnr_file_path",[(SRF_1_PATH,'cnrs1.txt',SRF_1_CNR_PATH), (SRF_2_PATH,'cnrs2.txt',SRF_2_CNR_PATH)])
def test_srf2corners(test_srf,filename,sample_cnr_file_path):
    srf.srf2corners(test_srf,cnrs=filename)
    out, err = shared.exe("diff -qr " + sample_cnr_file_path + " " + filename)
    assert out == "" and err == ""

@pytest.mark.parametrize("test_srf,expected_latlondepth",[(SRF_1_PATH, {'lat': -38.3354, 'depth': 0.0431, 'lon': 176.2414}),\
                                                          (SRF_2_PATH, {'lat': -37.1105, 'depth': 0.0381, 'lon': 176.7958}
)])
def test_read_latlondepth(test_srf,expected_latlondepth): #give you so many lat,lon,depth points

    points = srf.read_latlondepth(test_srf)
    assert points[9] == expected_latlondepth # 10th point in the srf file

@pytest.mark.parametrize("test_srf,lonlat,depth",[(SRF_1_PATH, True, True),(SRF_2_PATH, True, True),(SRF_1_PATH, True, False), \
                                                  (SRF_2_PATH, True, False),(SRF_1_PATH, False, True),(SRF_2_PATH, False, True), \
                                                  (SRF_1_PATH, False, False),(SRF_2_PATH, False, False)])
def test_get_hypo(test_srf,lonlat, depth):
    print srf.get_hypo(test_srf,lonlat = lonlat, depth = depth)