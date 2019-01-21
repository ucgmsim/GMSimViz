import numpy as np
import pytest

from gmsimviz import geo


@pytest.fixture
def rot_mat1():
    return geo.gen_mat(0, 170, -40)[0]


@pytest.fixture
def rot_mat2():
    return geo.gen_mat(-10, -170, 40)[0]


@pytest.mark.parametrize(
    "xy, expected",
    [
        (np.array([[0, 0]]), np.array([[170, -40]])),
        (
            np.array([[10, 0], [0, 10]]),
            np.array([[170.117267, -39.999942], [170.000000, -40.089832]]),
        ),
    ],
)
def test_xy2ll_mat1(rot_mat1, xy, expected):
    assert np.allclose(geo.xy2ll(xy, rot_mat1), expected, atol=1e-5)


@pytest.mark.parametrize(
    "xy, expected",
    [
        (np.array([[-20, -40]]), np.array([[-170.313905, 40.322247]])),
        (np.array([[-5000, 5000]]), np.array([[162.899713, -5.841725]])),
    ],
)
def test_xy2ll_mat2(rot_mat2, xy, expected):
    assert np.allclose(geo.xy2ll(xy, rot_mat2), expected, atol=1e-5)


@pytest.mark.parametrize(
    "gp, nx, ny, hh, expected",
    [
        (np.array([[0, 0]]), 100, 100, 0.1, np.array([[-4.95, -4.95]])),
        (np.array([[99, 99]]), 100, 100, 0.1, np.array([[4.95, 4.95]])),
        (np.array([[250, 9]]), 500, 20, 0.1, np.array([[0.05, -0.05]])),
    ],
)
def test_gp2xy(gp, nx, ny, hh, expected):
    assert np.allclose(geo.gp2xy(gp, nx, ny, hh), expected, atol=1e-5)


@pytest.mark.parametrize(
    "lon1, lat1, lon2, lat2, expected",
    [
        (170, -40, 180, -40, 852.309109),
        (170, -40, -180, -40, 852.309109),
        (1, 0.4, -1, -0.4, 239.787918),
    ],
)
def test_ll_dist(lon1, lat1, lon2, lat2, expected):
    assert geo.ll_dist(lon1, lat1, lon2, lat2) == pytest.approx(expected)


@pytest.mark.parametrize(
    "lon1, lat1, lon2, lat2, midpoint, expected",
    [
        (172, -43, 172, -40, False, 0.0),
        (-10, 5, -10, 1, False, 180.0),
        (180, 5, -180, -5, True, 180.0),
        (180, 5, -180, -5, False, 180.0),
        (170, -40, 180, -30, True, 39.074631),
        (10, 20, -10, -30, False, 201.598154),
    ],
)
def test_ll_bearing(lon1, lat1, lon2, lat2, midpoint, expected):
    assert geo.ll_bearing(lon1, lat1, lon2, lat2, midpoint=midpoint) == pytest.approx(
        expected
    )


@pytest.mark.parametrize(
    "b1, b2, expected",
    [
        (0, 45, 45),
        (0, 200, -160),
        (170, 190, 20),
        (190, 170, -20),
        (-10, 10, 20),
        (10, -1, -11),
        (270, 250, -20),
        (360, 359, -1),
    ],
)
def test_ll_bearing(b1, b2, expected):
    assert geo.angle_diff(b1, b2) == pytest.approx(expected)


@pytest.mark.parametrize("angles, expected", [
    ([[1, 100], [2, 100]], 1.5),
    ([[0, 100], [2, 50]], 0.666637),
    ([[0, 10000], [350, 50], [10, 50]], 0),
    ([[-10, 0.001], [10, 0.001]], 0),
    ([[-10, 100], [190, 50]], 332.122013),
])
def test_avg_wbearing(angles, expected):
    assert geo.avg_wbearing(angles) % 360 == pytest.approx(expected)
