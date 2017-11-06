#!/usr/bin/env python2
"""
Plot a single seismogram animation.
"""

import multiprocessing as mp
import os
import sys
from shutil import copy, copyfile

import numpy as np

import qcore.gmt as gmt
import qcore.timeseries as ts

try:
    velfile = os.path.abspath(sys.argv[1])
    basename = os.path.splitext(velfile)[0]
    gmt_temp = '%s/gmt_wd' % (os.path.dirname(velfile))
    if not os.path.exists(gmt_temp):
        os.makedirs(gmt_temp)
    assert(os.path.exists(velfile))
except IndexError:
    print('First parameter must be the velocity file.')
    exit(1)
except AssertionError:
    print('Could not find velocity file %s.' % (velfile))
    exit(1)

out_dir = os.path.abspath('Png')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

series = ts.read_ascii(velfile, t0 = True)
max_v = max(abs(series))
max_t = len(series) * 0.005
time = np.arange(0, max_t, 0.005)
xy = str(np.dstack((time, series)).tolist()) \
        .replace('], [', '\n').rstrip(']').lstrip('[').split('\n')

t = gmt.GMTPlot('%s/template.ps' % (gmt_temp))
gmt.gmt_defaults(wd = gmt_temp, font_annot_primary = '28,white', \
        extra = ['MAP_TICK_PEN_PRIMARY', 'thick,white', \
                'MAP_FRAME_PEN', 'thick,white', 'COLOR_BACKGROUND', 'white'])
t.background(11.2, 2.7, colour = 'black', y_margin = 0.5, x_margin = 0.1)
t.spacial('X', (0, max_t, -80, 80), \
        sizing = '10/2', x_shift = 1, y_shift = 1)
t.ticks(major = 40, minor = 1, sides = 'sw')
t.text(0, 80, 'cm/s', align = 'LT', dx = 0.1, dy = -0.08, size = 28, colour = 'white')
t.text(250, -80, 'sec', align = 'RB', dx = -0.1, dy = 0.08, size = 28, colour = 'white')
t.leave()

def render(i):
    wd = '%s/%.4d' % (gmt_temp, i)
    if not os.path.exists(wd):
        os.makedirs(wd)
    copyfile(t.pspath, '%s/%.4d.ps' % (wd, i))
    copy('%s/gmt.conf' % (os.path.dirname(t.pspath)), wd)
    copy('%s/gmt.history' % (os.path.dirname(t.pspath)), wd)
    p = gmt.GMTPlot('%s/%.4d.ps' % (wd, i), append = True)
    p.path('\n'.join(xy[:20 * i]), is_file = False, colour = 'white', \
            width = 1.2)
    p.finalise()
    p.png(out_dir = out_dir)

pool = mp.Pool(40)
pool.map(render, xrange(len(series) // 20))
