#!/usr/bin/env python2
"""
Plots deagg data.

Requires:
numpy >= 1.13 for numpy.unique(axis)

TODO:
box style legend or gmt style legend
automation in discretisation
"""

from argparse import ArgumentParser
from io import BytesIO
import math
import os
from subprocess import call

import numpy as np

from qcore import gmt
# DEBUG
import gmt

# TODO: automate defaults
dx = 10
dy = 0.25
x_len = 4.5
y_len = 4.0
z_len = 2.5

###
### LOAD DATA
###
parser = ArgumentParser()
parser.add_argument('deagg_file', help = 'deagg file to plot')
args = parser.parse_args()
assert(os.path.exists(args.deagg_file))
rrup_mag_e_c = np.loadtxt(args.deagg_file, skiprows = 4, usecols = (2, 1, 5, 4))

x_min = 0
x_max = int(math.ceil(max(rrup_mag_e_c[:, 0] / float(dx)))) * dx
x_inc = 20

y_min = 5
y_max = 9
y_inc = 0.5

# bins to put data in
bins_x = (np.arange(int(x_max / dx)) + 1) * dx
bins_y = (np.arange(int((y_max - y_min) / dy)) + 1) * dy + y_min
bins_z = np.array([-2, -1, -0.5, 0, 0.5, 1, \
                   max(2, np.max(rrup_mag_e_c[:, 2]) + 1)])
# convert data into bin indexes
rrup_mag_e_c[:, 0] = np.digitize(rrup_mag_e_c[:, 0], bins_x)
rrup_mag_e_c[:, 1] = np.digitize(rrup_mag_e_c[:, 1], bins_y)
rrup_mag_e_c[:, 2] = np.digitize(rrup_mag_e_c[:, 2], bins_z)
# combine duplicate bins
blocks = np.zeros(tuple(map(int, np.append(np.max(rrup_mag_e_c[:, :3], axis = 0) + 1, 2))))
unique = np.unique(rrup_mag_e_c[:, :3], axis = 0)
for rrup, mag, e in unique[unique[:, 2].argsort()]:
    # get base
    blocks[int(rrup), int(mag), int(e), 1] = \
            np.max(blocks[int(rrup), int(mag), :, 0])
    # current top = value + base
    value = np.add.reduce(rrup_mag_e_c[ \
            np.ix_(np.minimum.reduce( \
            rrup_mag_e_c[:, :3] == (rrup, mag, e), axis = 1), (3,))])
    if value:
        blocks[int(rrup), int(mag), int(e), 0] = \
                value + blocks[int(rrup), int(mag), int(e), 1]
# move indexes into array
# top base
blocks = blocks.reshape((-1, 2))

# z axis depends on max contribution tower
z_inc = int(math.ceil(np.max(np.add.reduce(blocks, axis = 2)) / 5.0))
z_max = z_inc * 5

###
### PLOT AXES
###
p = gmt.GMTPlot('deagg.ps')
os.remove('gmt.conf')
# setup axes
p.spacial('X', (x_min, x_max, y_min, y_max, 0, z_max), \
        sizing = '%si/%si' % (x_len, y_len), z = 'Z%si' % (z_len), \
        p = '150/30', x_shift = '2', y_shift = 2)
p.ticks(axis = 'x', major = x_inc, minor = None, \
        label = 'Rupture Distance (km)', sides = 's')
p.ticks(axis = 'y', major = y_inc, minor = None, \
        label = 'Magnitude', sides = 'e')
p.ticks(axis = 'z', major = z_inc, minor = None, gridline = z_inc, \
        label = '%Contribution', sides = 'z')
# GMT won't plot gridlines without box, manually add gridlines
gridlines = []
for z in xrange(z_inc, z_max + z_inc, z_inc):
    gridlines.append('%s %s %s\n%s %s %s\n%s %s %s' \
                     % (x_min, y_min, z, x_min, y_max, z, x_max, y_max, z))
gridlines.append('%s %s 0\n%s %s %s' % (x_min, y_max, x_min, y_max, z_max))
gridlines.append('%s %s 0\n%s %s %s' % (x_max, y_max, x_max, y_max, z_max))
p.path('\n>\n'.join(gridlines), is_file = False, width = '0.5p', z = True)

###
### PLOT CONTENTS
###
symbol = 'o%si/%sib' % (float(x_len) / len(bins_x) * 0.9, \
                        float(y_len) / len(bins_x) * 0.9)

x = np.arange(blocks.shape[0]) * dx + 0.5 * dx
y = np.arange(blocks.shape[1]) * dy + 0.5 * dy + y_min
elayer = np.column_stack((np.repeat(x, len(y)), np.tile(y, len(x)), blocks[:, :, -2].flatten(), np.zeros(len(x) * len(y))))
gmt_in = BytesIO()
np.savetxt(gmt_in, elayer, fmt = '%.6f')
print gmt_in.getvalue()
# X Y TOP CPT BASE
#p.points(gmt_in.getvalue(), is_file = False, shape = symbol, size = None, z = True, fill = 'blue', line = None)
p.points('85 8.125 33 -10000 11\n85 8.125 11 0 0', is_file = False, shape = symbol, z = True, size = None, line = None, cpt = '/home/nesi00213/PlottingData/cpt/nz_topo_grey1.cpt')
print symbol

###
### SAVE
###
p.finalise()
p.png(portrait = True, background = 'white', dpi = 300)
