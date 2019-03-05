#!/usr/bin/env python2
"""
MPI pattern based on
github.com/jbornschein/mpi4py-examples/blob/master/09-task-pull.py
"""

from hashlib import sha1
import os
from shutil import rmtree

import pytest
try:
    from imageio import imread
except ImportError:
    from scipy.misc import imread

from gmsimviz import gmt

TEMP_DIR = os.path.abspath('GMT_TESTING')
if os.path.exists(TEMP_DIR):
    rmtree(TEMP_DIR)
os.makedirs(TEMP_DIR)

###
### INTERNAL TESTING FUNCTIONS
###

def __test_coastlines(pf):
    p = gmt.GMTPlot(pf)
    p.spacial('M', (170, 175, -44, -38), sizing = 3)
    p.coastlines(width = '0.4p')
    p.finalise()
    p.png(dpi = 200, clip = True)

def __test_land(pf):
    p = gmt.GMTPlot(pf)
    p.spacial('M', (170.1, 179.91, -37, -34), sizing = 7)
    p.land(fill = 'darkred')
    p.finalise()
    p.png(dpi = 100, clip = True)

def __test_ticks(pf):
    p = gmt.GMTPlot(pf)
    p.spacial('M', (160.992, 174.9122, -44, -34.01), sizing = 2)
    p.ticks(major = '1d', minor = '20m', sides = 'ew')
    p.finalise()
    p.png(dpi = 100, clip = True)

def __test_ticks2(pf):
    p = gmt.GMTPlot(pf)
    p.spacial('T', (160.992, 174.9122, -44, -34.01), sizing = 5, \
            x_shift = 0.5, y_shift = 0.5, lon0 = 174.9122)
    p.ticks(major = '1d', minor = '20m', sides = 'ew')
    p.finalise()
    p.png(dpi = 100, clip = True)

def __test_cpt(pf):
    cptf = '%s.cpt' % (os.path.splitext(pf)[0])
    gmt.makecpt('hot', cptf, 0, 120, inc = 0.1, invert = True, \
            wd = os.path.dirname(pf))
    p = gmt.GMTPlot(pf)
    p.spacial('X', (0, 15, 0, 4), sizing = '15/4')
    p.cpt_scale(6.1, 2.05, cptf, 20, 5, label = 'test_scale', \
            length = 3.05, thickness = '0.3i')
    p.finalise()
    p.png(dpi = 320, clip = True)

def __test_cpt2(pf):
    cptf = '%s.cpt' % (os.path.splitext(pf)[0])
    gmt.makecpt('polar', cptf, -1.5, 1.5, inc = 0.25, invert = False, \
            wd = os.path.dirname(pf), bg = '0/0/80', fg = '80/0/0')
    p = gmt.GMTPlot(pf)
    p.spacial('X', (0, 4, 0, 2), sizing = '4/2', x_shift = 1, y_shift = 1)
    p.cpt_scale(0, 0, cptf, 0.5, 0.25, cross_tick = 0.5, align = 'LB', \
            length = 3, thickness = '0.3i', arrow_f = True, arrow_b = True)
    p.finalise()
    p.png(dpi = 222, clip = True)

def __test_fill(pf):
    p = gmt.GMTPlot(pf)
    gmt.gmt_defaults(wd = os.path.dirname(pf), ps_media = 'A5')
    p.background(1.5, 1)
    p.background(3, 2, x_margin = 1.5, colour = 'blue')
    p.background(1.5, 1, y_margin = 1, colour = 'red')
    p.background(3, 2, x_margin = 4.5, colour = 'firebrick')
    p.finalise()
    p.png(dpi = 100, clip = False)

def __test_autotick(wd):
    major, minor = gmt.auto_tick(170, 180, 10)
    if major != 1 or minor != 0.1:
        return False, '1 0.1 vs actual %s %s' % (major, minor)
    return True, None

def __test_autotick2(wd):
    major, minor = gmt.auto_tick(170, 170.04, 4)
    if major != 0.01 or minor != 0.001:
        return False, '0.01 0.001 vs actual %s %s' % (major, minor)
    return True, None

def __test_autotick3(wd):
    major, minor = gmt.auto_tick(-178, 178, 5)
    if major != 100 or minor != 10:
        return False, '1 0.1 vs actual %s %s' % (major, minor)
    return True, None


###
### ALL TESTS ARE RUN FROM HERE
###

@pytest.mark.parametrize("test, expected_hash, min_version", [
    (__test_coastlines, '61efdfbe4cd9bfd5a90ad86c67013c8d9494abc6', 5.0),
    (__test_land, '8c31ef6345aaad068cc4eb3e1d4b3f78cd6fc9a3', 5.0),
    (__test_ticks, 'd625bfb10464664d38f60537686999945eafafe7', 5.0),
    (__test_ticks2, '903830d2ba5d19a415a6d4e15fa62175bd7b7242', 5.0),
    (__test_cpt, '960032a10ef5fef8f013aab8892082fc4a606c9c', 5.2),
    (__test_cpt2, 'd83a3f811e8e514a6e88dacf38fe24bc4dea3280', 5.2),
    (__test_fill, 'eaa069d015eefb20f9826061a6de09d17a420a91', 5.0),
    (__test_autotick, None, 5.0),
    (__test_autotick2, None, 5.0),
    (__test_autotick3, None, 5.0),
])
def test_all(test, expected_hash, min_version):
    wd = os.path.join(TEMP_DIR, test.__name__)
    if not os.path.exists(wd):
        os.makedirs(wd)

    if float('.'.join(gmt.GMT_VERSION.split('.')[:2])) < min_version:
        return

    if expected_hash != None:
        ps = '%s/%s.ps' % (wd, test.__name__)
        png = '%s/%s.png' % (wd, test.__name__)
        test(ps)
        os.symlink(png, os.path.join(TEMP_DIR, os.path.basename(png)))
        assert os.path.isfile(png)
        png_hash = sha1(imread(png)).hexdigest()
        assert png_hash == expected_hash
    else:
        assert test(wd)[0]
