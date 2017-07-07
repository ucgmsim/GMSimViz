#!/usr/bin/env python2
"""
1.
Creates binary files with LONG, LAT, VALUE for each segment.
Same format as a timeslice to be plotted in GMT.
Creates corresponding path files for masking area and GMT mask files.

2.
Creates a standard corners file to plot fault planes.

3.
Creates a GMT plot of the SRF file.

USAGE:
with internally set defaults:
plot_srf.py
standard usage:
plot_srf.py <srf_file> <dpi_override> <'active_faults' to plot faults>
plot_srf.py mysrf.srf 300 active_faults
"""

from math import log10
import os
from shutil import rmtree
import sys

import numpy as np

import qcore_path
from gmt import *
from srf import *
from geo import *

# illumination file should be in the same directory
topo = '/nesi/projects/nesi00213/PlottingData/Topo/srtm_all_filt_nz.grd'
topo_low = '/nesi/projects/nesi00213/PlottingData/Topo/nztopo.grd'
faults = '/nesi/projects/nesi00213/PlottingData/Paths/faults/FAULTS_20161219.ll'
cpt = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
        'cpt', 'slip.cpt')

# can specify here or pass as command line argument
srf = 'default.srf'
dpi = 300
plot_faults = False
if len(sys.argv) > 1:
    srf = sys.argv[1]
if len(sys.argv) > 2:
    dpi = int(sys.argv[2])
if len(sys.argv) > 3 and sys.argv[3] == 'active_faults':
    plot_faults = True

srf = os.path.abspath(srf)
if not os.path.exists(srf):
    print('SRF file %s not found.' % (srf))
    exit(1)
srf_dir = os.path.dirname(srf)
if srf_dir == '':
    srf_dir = '.'

dx, dy = srf_dxy(srf)
text_dx = '%s km' % (dx)
text_dy = '%s km' % (dy)
# plot at greater resolution to increase smoothness
# also considering rotation, grid will not be exactly matching
plot_dx = '%sk' % (dx * 0.6)
plot_dy = '%sk' % (dy * 0.6)

# output directory for srf resources
out_dir = os.path.abspath('GMT_WD_SRF')
if os.path.exists(out_dir):
    rmtree(out_dir)
os.makedirs(out_dir)
# output for plane data
os.makedirs(os.path.join(out_dir, 'PLANES'))

###
### OUTPUT 1: binary file for GMT grid plotting
###
print('Loading SRF file data from %s...' % (srf))
# get all corners
bounds = get_bounds(srf)
# get all tinit values, set a sane countour interval
# contour interval should probably also depend on area
tinit = srf2llv_py(srf, value = 'tinit')
tinit_max = max([np.max(tinit[p][:, 2]) for p in xrange(len(bounds))])
contour_int = 2
if tinit_max < 10:
    contour_int = 1
if tinit_max < 2:
    contour_int = 0.3
# gridding is computationally expensive
# each segment should include minimal region
seg_regions = []
# for percentile to automatically calculate colour palette range
values = np.array([], dtype = np.float32)
# find total extremes
np_bounds = np.array(bounds)
x_min, y_min = np.min(np.min(np_bounds, axis = 0), axis = 0)
x_max, y_max = np.max(np.max(np_bounds, axis = 0), axis = 0)
plot_region = (x_min - 0.1, x_max + 0.1, y_min - 0.1, y_max + 0.1)
# for plotting region on NZ-wide map
plot_bounds = '%f %f\n%f %f\n%f %f\n%f %f\n' % \
        (plot_region[0], plot_region[2], plot_region[1], plot_region[2], \
        plot_region[1], plot_region[3], plot_region[0], plot_region[3])
# read all max slip values (all at once is much faster)
seg_llslips = srf2llv_py(srf, value = 'slip')
for seg in xrange(len(bounds)):
    # create binary llv file for GMT overlay
    seg_llslips[seg].astype(np.float32).tofile( \
            '%s/PLANES/slip_map_%d.bin' % (out_dir, seg))
    values = np.append(values, seg_llslips[seg][:, -1])
    # also store tinit values retrieved previously
    tinit[seg].astype(np.float32).tofile( \
            '%s/PLANES/tinit_map_%d.bin' % (out_dir, seg))
    # create a mask path for GMT overlay
    path_from_corners(corners = bounds[seg], min_edge_points = 100, \
            output = '%s/PLANES/plane_%d.bounds' % (out_dir, seg))
    # create mask from path
    x_min, y_min = np.min(np_bounds[seg], axis = 0)
    x_max, y_max = np.max(np_bounds[seg], axis = 0)
    seg_regions.append((x_min, x_max, y_min, y_max))
    grd_mask('%s/PLANES/plane_%d.bounds' % (out_dir, seg), \
            '%s/PLANES/plane_%d.mask' % (out_dir, seg), \
            dx = plot_dx, dy = plot_dy, \
            region = seg_regions[seg])
percentile = np.percentile(values, 95)
maximum = np.max(values)
average = np.average(values)
subfaults = len(values)
# round percentile significant digits for colour pallete
if percentile < 1000:
    # 1 sf
    cpt_max = round(percentile, -int(log10(abs(percentile))))
else:
    # 2 sf
    cpt_max = round(percentile, -int(log10(abs(percentile))) + 1)
print('Loading complete.')


###
### OUTPUT 2: corners file for fault plane and hypocentre plot
###
print('Creating corners file...')
# find hypocentre, use bounds from previous step
hypocentre = get_hypo(srf)
# standard corners file format
with open('%s/corners.txt' % (out_dir), 'w') as cf:
    cf.write('>hypocentre\n')
    cf.write('%f %f\n' % (hypocentre[0], hypocentre[1]))
    cf.write('>corners for each segment\n')
    for c, corners in enumerate(bounds):
        cf.write('>segment %d\n' % (c))
        for i in xrange(5):
            cf.write('%f %f\n' % tuple(corners[i % 4]))
print('Corners written to %s.' % ('%s/corners.txt' % (out_dir)))


###
### OUTPUT 3: GMT MAP
###
print('Plotting SRF on map...')
makecpt(cpt, '%s/slip.cpt' % (out_dir), 0, cpt_max, 1)
# also make a cpt to properly stretch topo colours
topo_cpt = '%s/topo.cpt' % (out_dir)
makecpt('gray', topo_cpt, -10000, 3000, inc = 10)
gmt_defaults(wd = out_dir)
# gap on left of maps
gap = 1
# width of NZ map, if changed, other things also need updating
# including tick font size and or tick increment for that map
full_width = 4

### PART A: zoomed in map
p = GMTPlot('%s/%s_map.ps' % (out_dir, \
        os.path.splitext(os.path.basename(srf))[0]))
# cover a larger space, will get cropped at the end
p.background(20, 15)
# this is how high the NZ map will end up being
full_height = mapproject(nz_region[0], nz_region[3], region = nz_region, \
        projection = 'M%s' % (full_width), wd = out_dir)[1]
# match height of zoomed in map with full size map
zoom_width, zoom_height = \
        map_width('M', full_height, plot_region, wd = out_dir)
p.spacial('M', plot_region, sizing = zoom_width, \
        x_shift = gap, y_shift = 2.5)
p.land()
p.topo(topo, cpt = topo_cpt)
p.water()
if plot_faults:
    p.path(faults, is_file = True, close = False, width = '0.4p', colour = 'red')
for seg in xrange(len(bounds)):
    p.overlay('%s/PLANES/slip_map_%d.bin' % (out_dir, seg), \
            '%s/slip.cpt' % (out_dir), dx = plot_dx, dy = plot_dy, \
            climit = 2, crop_grd = '%s/PLANES/plane_%s.mask' % (out_dir, seg), \
            land_crop = False, custom_region = seg_regions[seg], \
            transparency = 30)
    p.overlay('%s/PLANES/tinit_map_%d.bin' % (out_dir, seg), None, \
            dx = plot_dx, dy = plot_dy, climit = 2, \
            crop_grd = '%s/PLANES/plane_%s.mask' % (out_dir, seg), \
            land_crop = False, custom_region = seg_regions[seg], \
            transparency = 30, contours = contour_int)
p.fault(srf, is_srf = True, hyp_colour = 'red')

p.coastlines()
p.sites(sites.keys())
# work out an ideal tick increment (ticks per inch)
# x axis is more constrainig
major_tick = 0.05
while ((plot_region[1] - plot_region[0]) / major_tick) / zoom_width > 1.4:
    major_tick += 0.05
p.ticks(major = '%sd' % (major_tick), \
        minor = '%sd' % (major_tick / 2.0), sides = 'ws')

### PART B: NZ map
# draw NZ wide map to show rough location in country
p.spacial('M', nz_region, sizing = full_width, \
    x_shift = zoom_width + gap)
# height of NZ map
full_height = mapproject(nz_region[0], nz_region[3], wd = out_dir)[1]
p.land()
p.topo(topo_low, cpt = topo_cpt)
p.water()
if plot_faults:
    p.path(faults, is_file = True, close = False, width = '0.1p', colour = 'red')
p.path(plot_bounds, is_file = False, close = True, colour = 'black')
# get displacement of box to draw zoom lines later
window_bottom = mapproject(plot_region[1], plot_region[2], wd = out_dir)
window_top = mapproject(plot_region[1], plot_region[3], wd = out_dir)
# also draw fault planes / hypocentre
p.fault(srf, is_srf = True, plane_width = '0.4p', top_width = '0.6p', \
        hyp_colour = 'red')
p.coastlines(width = '0.1p')
p.ticks(major = '2d', minor = '30m', sides = 'ws')

### PART C: zoom lines
# draw zoom lines that extend from view box to original plot
p.spacial('X', \
        (0, window_bottom[0] + gap, 0, max(zoom_height, full_height)), \
        x_shift = -gap, sizing = '%s/%s' % \
        (window_top[0] + gap, max(zoom_height, full_height)))
p.path('%f %f\n%f %f\n' % \
        (0, 0, window_bottom[0] + gap, window_bottom[1]), \
        width = '0.6p', is_file = False, split = '-', \
        straight = True, colour = 'black')
p.path('%f %f\n%f %f\n' % \
        (0, zoom_height, window_top[0] + gap, window_top[1]), \
        width = '0.6p', is_file = False, split = '-', \
        straight = True, colour = 'black')

### PART D: surrounding info
# add text and colour palette
# position to enclose both plots
total_width = zoom_width + gap + full_width
total_height = max(zoom_height, full_height)
p.spacial('X', (0, total_width, 0, total_height + 2), \
        sizing = '%s/%s' % (total_width, total_height + 2), \
        x_shift = -zoom_width)
# SRF filename
p.text(total_width / 2.0, total_height, os.path.basename(srf), \
        align = 'CB', size = '20p', dy = 0.8)
# max slip
p.text(zoom_width / 2.0, total_height, 'Maximum slip: ', \
        align = 'RB', size = '14p', dy = 0.5)
p.text(zoom_width / 2.0 + 0.1, total_height, '%.1f cm' % (maximum), \
        align = 'LB', size = '14p', dy = 0.5)
# 95th percentile
p.text(zoom_width / 2.0, total_height, '95th percentile: ', \
        align = 'RB', size = '14p', dy = 0.3)
p.text(zoom_width / 2.0 + 0.1, total_height, '%.1f cm' % (percentile), \
        align = 'LB', size = '14p', dy = 0.3)
# average slip
p.text(zoom_width / 2.0, total_height, 'Average slip: ', \
        align = 'RB', size = '14p', dy = 0.1)
p.text(zoom_width / 2.0 + 0.1, total_height, '%.1f cm' % (average), \
        align = 'LB', size = '14p', dy = 0.1)
# planes
p.text(total_width - 4 / 2.0, total_height, 'Planes: ', \
        align = 'RB', size = '14p', dy = 0.5)
p.text(total_width - 4 / 2.0 + 0.1, total_height, len(bounds), \
        align = 'LB', size = '14p', dy = 0.5)
# dx and dy
p.text(total_width - 4 / 2.0, total_height, 'dX, dY: ', \
        align = 'RB', size = '14p', dy = 0.3)
p.text(total_width - 4 / 2.0 + 0.1, total_height, \
        '%s, %s' % (text_dx, text_dy), \
        align = 'LB', size = '14p', dy = 0.3)
# subfaults
p.text(total_width - 4 / 2.0, total_height, 'Subfaults: ', \
        align = 'RB', size = '14p', dy = 0.1)
p.text(total_width - 4 / 2.0 + 0.1, total_height, subfaults, \
        align = 'LB', size = '14p', dy = 0.1)
# scale
p.cpt_scale(zoom_width / 2.0, -0.5, '%s/slip.cpt' % (out_dir), \
        cpt_max / 4.0, cpt_max / 8.0, \
        label = 'Slip (cm)', length = zoom_width)

p.finalise()
p.png(dpi = dpi, out_dir = srf_dir)
rmtree(out_dir)
print('SRF plot at %s.' % ('%s/srf_map.png' % (out_dir)))
