#!/usr/bin/env python2
"""
Plots deagg data.

Requires:
numpy >= 1.13 for numpy.unique(axis)

TODO:
automation in discretisation
"""

from argparse import ArgumentParser
from io import BytesIO
import math
import os
from tempfile import mkstemp

import numpy as np
# requires fairly new version of numpy for axis parameter in np.unique
npv = map(int, np.__version__.split('.'))
if npv[0] < 1 or (npv[0] == 1 and npv[1] < 13):
    print('requires numpy >= 1.13')
    exit(1)

from qcore import gmt

X_LEN = 4.5
Y_LEN = 4.0
Z_LEN = 2.5
ROT = 30
TILT = 60
LEGEND_SPACE = 0.7
EPSILON_COLOURS = ['215/38/3', '252/94/62', '252/180/158', '254/220/210', \
                   '217/217/255', '151/151/255', '0/0/255', '0/0/170']

###
### LOAD DATA
###
parser = ArgumentParser()
parser.add_argument('deagg_file', help = 'deagg file to plot')
args = parser.parse_args()
assert(os.path.exists(args.deagg_file))
rrup_mag_e_c = np.loadtxt(args.deagg_file, skiprows = 4, usecols = (2, 1, 5, 4))


###
### PROCESS DATA
###
# TODO: automate, logic from given original code missing
dx = 10
dy = 0.25

# x axis
x_max = int(math.ceil(max(rrup_mag_e_c[:, 0] / float(dx)))) * dx
if x_max < 115:
    x_inc = 10
elif x_max < 225:
    x_inc = 20
elif x_max < 335:
    x_inc = 30
elif x_max < 445:
    x_inc = 40
else:
    x_inc = 50

# y axis
y_min = 5
y_max = 9
if y_max - y_min < 5:
    y_inc = 0.5
else:
    y_inc = 1.0

# bins to put data in
bins_x = (np.arange(int(x_max / dx)) + 1) * dx
bins_y = (np.arange(int((y_max - y_min) / dy)) + 1) * dy + y_min
bins_z = np.array([-2, -1, -0.5, 0, 0.5, 1, 2, \
                   max(3, np.max(rrup_mag_e_c[:, 2]) + 1)])

# convert data into bin indexes
rrup_mag_e_c[:, 0] = np.digitize(rrup_mag_e_c[:, 0], bins_x)
rrup_mag_e_c[:, 1] = np.digitize(rrup_mag_e_c[:, 1], bins_y)
rrup_mag_e_c[:, 2] = np.digitize(rrup_mag_e_c[:, 2], bins_z)

# combine duplicate bins
# TODO: change datatype to int?
blocks = np.zeros(tuple(map(int, \
                  np.append(np.max(rrup_mag_e_c[:, :3], axis = 0) + 1, 2))))
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
del rrup_mag_e_c, unique

# move indexes into array
top, base = blocks.reshape((-1, 2)).T
cpt = np.tile(np.arange(blocks.shape[2]), np.prod(blocks.shape[0:2]))
y = np.tile(np.repeat(np.arange(blocks.shape[1]) * dy + 0.5 * dy + y_min, \
                      blocks.shape[2]), blocks.shape[0])
x = np.repeat(np.arange(blocks.shape[0]) * dx + 0.5 * dx, \
              np.prod(blocks.shape[1:3]))
gmt_rows = np.column_stack((x, y, top, cpt, base))
del x, y, top, cpt, base
# don't plot if top == 0
gmt_rows = np.delete(gmt_rows, \
                     np.argwhere(gmt_rows[:, 2] == 0).flatten(), axis = 0)

# z axis depends on max contribution tower
z_inc = int(math.ceil(np.max(np.add.reduce(blocks, axis = 2)) / 5.0))
z_max = z_inc * 5
del blocks

###
### PLOT AXES
###
p = gmt.GMTPlot('deagg.ps')
os.remove('gmt.conf')
# setup axes
p.spacial('X', (0, x_max, y_min, y_max, 0, z_max), \
        sizing = '%si/%si' % (X_LEN, Y_LEN), z = 'Z%si' % (Z_LEN), \
        p = '%s/%s' % (180 - ROT, 90 - TILT), x_shift = '2', y_shift = 2)
p.ticks(axis = 'x', major = x_inc, minor = None, \
        label = 'Rupture Distance (km)', sides = 's')
p.ticks(axis = 'y', major = y_inc, minor = None, \
        label = 'Magnitude', sides = 'e')
p.ticks(axis = 'z', major = z_inc, minor = None, gridline = z_inc, \
        label = '%Contribution', sides = 'z')
# GMT won't plot gridlines without box, manually add gridlines
gridlines = []
for z in xrange(z_inc, z_max + z_inc, z_inc):
    gridlines.append('0 %s %s\n0 %s %s\n%s %s %s' \
                     % (y_min, z, y_max, z, x_max, y_max, z))
gridlines.append('0 %s 0\n0 %s %s' % (y_max, y_max, z_max))
gridlines.append('%s %s 0\n%s %s %s' % (x_max, y_max, x_max, y_max, z_max))
p.path('\n>\n'.join(gridlines), is_file = False, width = '0.5p', z = True)

###
### PLOT CONTENTS
###
cpt = mkstemp(suffix = '.cpt')[1]
gmt.makecpt(','.join(EPSILON_COLOURS), cpt, \
            0, len(EPSILON_COLOURS), inc = 1)
gmt_in = BytesIO()
np.savetxt(gmt_in, gmt_rows, fmt = '%.6f')
p.points(gmt_in.getvalue(), is_file = False, z = True, line = 'black', \
        shape = 'o', size = '%si/%sib' % (float(X_LEN) / len(bins_x) - 0.05, \
                                          float(Y_LEN) / len(bins_x) - 0.05), \
        line_thickness = '0.5p', cpt = cpt)

###
### PLOT LEGEND
###
# x y diffs from start to end
angle = math.radians(ROT)
x_end = (X_LEN + math.cos(angle) * math.sin(angle) \
                 * (Y_LEN - math.tan(angle) * X_LEN)) / X_LEN
y_end = math.tan(angle) * x_end * X_LEN * (y_max - y_min) / Y_LEN
x_end *= x_max
# x y diffs at start
dip = (LEGEND_SPACE) / math.cos(math.radians(TILT)) + math.sin(angle) * X_LEN
x0 = dip * math.sin(angle) * (x_max / X_LEN)
y0 = y_min - dip * math.cos(angle) * ((y_max - y_min) / Y_LEN)
# legend cube definitions
legend_boxes = []
for i, x in enumerate(np.arange(0, 1.01, 1.0 / (len(EPSILON_COLOURS) - 1.0))):
    legend_boxes.append('%s %s %s %s' % (x0 + x * x_end, y0 + x * y_end, \
                                         z_inc / 2.0, i))
# cubes of legend
p.points('\n'.join(legend_boxes), is_file = False, z = True, line = 'black', \
        shape = 'o', size = '%si/%sib0' % (Z_LEN / 10.0, Z_LEN / 10.0), \
        line_thickness = '0.5p', cpt = cpt, clip = False)
os.remove(cpt)

###
### SAVE
###
p.finalise()
p.png(portrait = True, background = 'white', dpi = 300)
