import os

import numpy as np
import pytest

from gmsimviz import xyts


@pytest.fixture
def xyts_instance():
    return xyts.XYTSFile(
        os.path.join(os.path.dirname(__file__), os.pardir, "sample_data", "xyts.e3d")
    )


@pytest.mark.parametrize(
    "attr, expected",
    [
        ("x0", 0),
        ("y0", 0),
        ("z0", 1),
        ("t0", 0),
        ("nx", 75),
        ("ny", 72),
        ("nz", 1),
        ("nt", 694),
        ("dx", 2),
        ("dy", 2),
        ("hh", 0.4),
        ("dt", 0.1),
        ("mrot", 44),
        ("mlat", -44.222310),
        ("mlon", 170.962890),
        ("dxts", 5),
        ("dyts", 5),
        ("nx_sim", 375),
        ("ny_sim", 360),
        ("dip", 0),
        ("comps", {"X": 2.338741, "Y": 0.767945, "Z": 1.570796}),
    ],
)
def test_attribute(xyts_instance, attr, expected):
    assert getattr(xyts_instance, attr) == pytest.approx(expected)


@pytest.mark.parametrize(
    "gmt_format, expected",
    [
        (
            False,
            [
                [170.914383, -43.291595],
                [172.262466, -44.217743],
                [171.012955, -45.153000],
                [169.663437, -44.212143],
            ],
        ),
        (
            True,
            "170.914382935 -43.291595459\n172.262466431 -44.2177429199\n171.012954712 -45.1529998779\n169.66343689 -44.2121429443",
        ),
    ],
)
def test_corners(xyts_instance, gmt_format, expected):
    if gmt_format:
        assert xyts_instance.corners(gmt_format)[1] == expected
    else:
        assert np.allclose(xyts_instance.corners(gmt_format), expected, atol=1e-5)


@pytest.mark.parametrize(
    "expected",
    [(169.66343688964844, 172.26246643066406, -45.15299987792969, -43.291595458984375)],
)
def test_region(xyts_instance, expected):
    assert xyts_instance.region() == expected
