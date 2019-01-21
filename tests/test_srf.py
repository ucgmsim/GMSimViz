import os

import pytest

from gmsimviz import srf


@pytest.fixture
def srf_file():
    return os.path.join(
        os.path.dirname(__file__), os.pardir, "sample_data", "fault.srf"
    )


def test_file(srf_file):
    assert os.path.exists(srf_file)


@pytest.fixture
def srf_header(srf_file):
    return srf.read_header(srf_file, idx=True)


@pytest.mark.parametrize(
    "attr, expected",
    [
        ("centre", (170.9208, -44.2494)),
        ("nstrike", 348),
        ("ndip", 173),
        ("length", 34.83),
        ("width", 17.32),
        ("strike", 336),
        ("dip", 60),
        ("shyp", -10),
        ("dhyp", 10.39),
        ("dtop", 0),
    ],
)
def test_read_header(srf_header, attr, expected):
    assert srf_header[0][attr] == pytest.approx(expected)


def test_get_bounds(srf_file):
    assert srf.get_bounds(srf_file)[0] == [
        (171.0117, -44.3912),
        (170.8322, -44.107),
        (170.9298, -44.0743),
        (171.1097, -44.3583),
    ]


@pytest.mark.parametrize(
    "lonlat, depth, expected",
    [
        (True, False, (171.0319, -44.311)),
        (True, True, (171.0319, -44.311, 8.9734)),
        (False, True, (0, 7.415, 10.39)),
        (False, False, (0, 7.415)),
    ],
)
def test_get_hypo(srf_file, lonlat, depth, expected):
    assert srf.get_hypo(srf_file, lonlat=lonlat, depth=depth) == pytest.approx(expected)


def test_srf_dxy(srf_file):
    assert srf.srf_dxy(srf_file) == (0.1, 0.1)


def test_srf_dt(srf_file):
    assert srf.srf_dt(srf_file) == 0.025
