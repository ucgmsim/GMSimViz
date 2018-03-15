from qcore import srf, shared
import pytest

import os
""" Command to run this test: 'python -m pytest -v -s test_srf.py'  """

SRF_1_PATH = os.path.join(os.getcwd(), "sample1/Hossack_HYP01-01_S1244.srf")
SRF_2_PATH = os.path.join(os.getcwd(), "sample2/Tuakana13_HYP01-01_S1244.srf")
SRF_1_CNR_PATH = os.path.join(os.getcwd(), "sample1/output/cnrs.txt")
SRF_2_CNR_PATH = os.path.join(os.getcwd(), "sample2/output/cnrs.txt")
SRF_1_PLANES = srf.read_header(SRF_1_PATH, True)
SRF_2_PLANES = srf.read_header(SRF_2_PATH, True)
HEADERS = ['centre', 'nstrike', 'ndip', 'length', 'width', 'strike', 'dip', 'dtop', \
           'shyp', 'dhyp']

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

# TODO: check if sample cnrs files are OK to be used as sample output
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

@pytest.mark.parametrize("test_srf,seg,depth,expected_bounds",[(SRF_1_PATH, -1, True,[[(176.2493, -38.3301, 0.0431), (176.2202, -38.3495, 0.0431), (176.1814, -38.3221, 7.886), (176.2105, -38.3027, 7.886)]]
),(SRF_2_PATH, -1, True,[[(176.7922, -37.118, 0.0381), (176.8101, -37.0806, 0.0381), (176.876, -37.1089, 7.8931), (176.8581, -37.1464, 7.8931)], [(176.8107, -37.0798, 0.038), (176.8433, -37.0455, 0.038), (176.9092, -37.0739, 7.8672), (176.8765, -37.1082, 7.8672)]]), \
                                               (SRF_1_PATH, -1, False,[[(176.2493, -38.3301), (176.2202, -38.3495), (176.1814, -38.3221), (176.2105, -38.3027)]]
),(SRF_2_PATH, -1, False,[[(176.7922, -37.118), (176.8101, -37.0806), (176.876, -37.1089), (176.8581, -37.1464)], [(176.8107, -37.0798), (176.8433, -37.0455), (176.9092, -37.0739), (176.8765, -37.1082)]]
)])
def test_get_bounds(test_srf, seg, depth, expected_bounds):
    assert srf.get_bounds(test_srf, seg = seg, depth = depth) == expected_bounds

@pytest.mark.parametrize("test_srf, expected_nseg",[(SRF_1_PATH, 1),(SRF_2_PATH,2)])
def test_get_nseg(test_srf, expected_nseg):
    assert srf.get_nseg(test_srf) == expected_nseg

@pytest.mark.parametrize("test_srf, expected_result",[(SRF_1_PATH, True),(SRF_2_PATH,True)])
def test_is_ff(test_srf, expected_result):
    assert srf.is_ff(test_srf) == expected_result

@pytest.mark.parametrize("test_srf_planes, expected_result",[(SRF_1_PLANES, 1),(SRF_2_PLANES,2)])
def test_nplane1(test_srf_planes, expected_result):
    assert len(test_srf_planes) == expected_result

# # TODO:TO BE DECIDED LATER
# @pytest.mark.parametrize("test_srf,value,depth",[(SRF_1_PATH, None, True),(SRF_2_PATH, None, True),\
#                                                   (SRF_1_PATH, None, False),(SRF_2_PATH, None, False)])
# def test_get_lonlat(test_srf, value , depth):
#     print srf.get_lonlat(test_srf,value = value, depth  = depth)

# @pytest.mark.parametrize("test_srf,lonlat,depth",[(SRF_1_PATH, True, True),(SRF_2_PATH, True, True),(SRF_1_PATH, True, False), \
#                                                   (SRF_2_PATH, True, False),(SRF_1_PATH, False, True),(SRF_2_PATH, False, True), \
#                                                   (SRF_1_PATH, False, False),(SRF_2_PATH, False, False)])
# def test_get_hypo(test_srf,lonlat, depth):
#     print srf.get_hypo(test_srf,lonlat = lonlat, depth = depth)

# def test_pf_params():
#     with pytest.raises(AssertionError) as e:
#         print srf.ps_params(SRF_2_PATH)
#     # assert ("Error parsing documentation",) == e.value.args