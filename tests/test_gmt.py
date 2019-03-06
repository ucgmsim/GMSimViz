#!/usr/bin/env python2

import os
from shutil import rmtree

import pytest

try:
    from imageio import imread
except ImportError:
    from scipy.misc import imread

from gmsimviz import gmt

TEMP_DIR = os.path.abspath("gmt_output")
if os.path.exists(TEMP_DIR):
    rmtree(TEMP_DIR)
os.makedirs(TEMP_DIR)


class Figure:
    """
    Creates an object with a savefig function for pytest-mpl.
    """

    # to make matplotlib.pyplot.close(Figure) work
    int = None

    def __init__(self, gmt_plot, dpi=100, clip=True):
        self.p = gmt_plot
        self.dpi = dpi
        self.clip = clip

    def savefig(self, filename):
        self.p.finalise()
        self.p.png(out_name=os.path.splitext(filename)[0], dpi=self.dpi, clip=self.clip)
        rmtree(self.p.wd)


def gmt_plot_factory(test_name):
    """
    Creates a gmt plot object with a unique name for the test.
    """
    ps_name = os.path.join(TEMP_DIR, test_name, "{}.ps".format(test_name))
    if not os.path.isdir(os.path.dirname(ps_name)):
        os.makedirs(os.path.dirname(ps_name))
    return gmt.GMTPlot(ps_name)


@pytest.fixture
def wd(request):
    temp_dir = os.path.join(TEMP_DIR, request.function.__name__)
    if not os.path.isdir(temp_dir):
        os.path.makedirs(temp_dir)
    return temp_dir


@pytest.fixture
def p(request):
    return gmt_plot_factory(request.function.__name__)


@pytest.mark.mpl_image_compare(tolerance=20)
def test_coastlines(p):
    p.spacial("M", (170, 175, -44, -38), sizing=3)
    p.coastlines(width="0.4p")
    return Figure(p, dpi=200)


@pytest.mark.mpl_image_compare(tolerance=20)
def test_land(p):
    p.spacial("M", (170.1, 179.91, -37, -34), sizing=7)
    p.land(fill="darkred")
    return Figure(p)


@pytest.mark.mpl_image_compare(tolerance=20)
def test_ticks(p):
    p.spacial("M", (160.992, 174.9122, -44, -34.01), sizing=2)
    p.ticks(major="1d", minor="20m", sides="ew")
    return Figure(p)


@pytest.mark.mpl_image_compare(tolerance=20)
def test_ticks2(p):
    p.spacial(
        "T",
        (160.992, 174.9122, -44, -34.01),
        sizing=5,
        x_shift=0.5,
        y_shift=0.5,
        lon0=174.9122,
    )
    p.ticks(major="1d", minor="20m", sides="ew")
    return Figure(p)


@pytest.mark.mpl_image_compare(tolerance=20)
def test_cpt(p, wd):
    cptf = os.path.join(wd, "cpt.cpt")
    gmt.makecpt("hot", cptf, 0, 120, inc=0.1, invert=True, wd=wd)
    p.spacial("X", (0, 15, 0, 4), sizing="15/4")
    p.cpt_scale(
        6.1, 2.05, cptf, 20, 5, label="test_scale", length=3.05, thickness="0.3i"
    )
    return Figure(p, dpi=320)


@pytest.mark.mpl_image_compare(tolerance=20)
def test_cpt2(p, wd):
    cptf = os.path.join(wd, "cpt.cpt")
    gmt.makecpt(
        "polar",
        cptf,
        -1.5,
        1.5,
        inc=0.25,
        invert=False,
        wd=wd,
        bg="0/0/80",
        fg="80/0/0",
    )
    p.spacial("X", (0, 4, 0, 2), sizing="4/2", x_shift=1, y_shift=1)
    p.cpt_scale(
        0,
        0,
        cptf,
        0.5,
        0.25,
        cross_tick=0.5,
        align="LB",
        length=3,
        thickness="0.3i",
        arrow_f=True,
        arrow_b=True,
    )
    return Figure(p, dpi=222)


@pytest.mark.mpl_image_compare(tolerance=20)
def test_fill(p, wd):
    gmt.gmt_defaults(wd=wd, ps_media="A5")
    p.background(1.5, 1)
    p.background(3, 2, x_margin=1.5, colour="blue")
    p.background(1.5, 1, y_margin=1, colour="red")
    p.background(3, 2, x_margin=4.5, colour="firebrick")
    return Figure(p, clip=False)


def test_autotick():
    major, minor = gmt.auto_tick(170, 180, 10)
    assert major == 1 and minor == 0.1


def test_autotick2():
    major, minor = gmt.auto_tick(170, 170.04, 4)
    assert major == 0.01 and minor == 0.001


def test_autotick3():
    major, minor = gmt.auto_tick(-178, 178, 5)
    assert major == 100 and minor == 10
