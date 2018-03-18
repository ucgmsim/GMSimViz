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

from math import floor, log10
import os
from shutil import rmtree
import sys
from tempfile import mkdtemp

import numpy as np

import qcore.geo as geo
import qcore.gmt as gmt
import qcore.srf as srf

faults = '/nesi/projects/nesi00213/PlottingData/Paths/faults/FAULTS_20161219.ll'
cpt = gmt.CPTS['slip']

# can specify here or pass as command line argument
srf_file = 'default.srf'
dpi = 300
plot_faults = False
if len(sys.argv) > 1:
    srf_file = sys.argv[1]
if len(sys.argv) > 2:
    dpi = int(sys.argv[2])
if len(sys.argv) > 3 and sys.argv[3] == 'active_faults':
    plot_faults = True

srf_file = os.path.abspath(srf_file)
if not os.path.exists(srf_file):
    print('SRF file %s not found.' % (srf_file))
    exit(1)
srf_dir = os.path.dirname(srf_file)
if srf_dir == '':
    srf_dir = '.'
# whether we are plotting a finite fault or point source
finite_fault = srf.is_ff(srf_file)

if finite_fault:
    dx, dy = srf.srf_dxy(srf_file)
    text_dx = '%s km' % (dx)
    text_dy = '%s km' % (dy)
    # plot at greater resolution to increase smoothness
    # also considering rotation, grid will not be exactly matching
    plot_dx = '%sk' % (dx * 0.6)
    plot_dy = '%sk' % (dy * 0.6)
else:
    text_dx = 'N/A'
    text_dy = 'N/A'

# output directory for srf resources
out_dir = os.path.abspath(mkdtemp(prefix = '_GMT_WD_SRF_', dir = '.'))
if finite_fault:
    # output for plane data
    os.makedirs(os.path.join(out_dir, 'PLANES'))

###
### OUTPUT 1: binary file for GMT grid plotting
###
if finite_fault:
    print('Loading SRF file data from %s...' % (srf_file))
    # get all corners
    bounds = srf.get_bounds(srf_file)
    # get all tinit values, set a sane countour interval
    # contour interval should probably also depend on area
    tinit = srf.srf2llv_py(srf_file, value = 'tinit')
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
    # read all max slip values (all at once is much faster)
    seg_llslips = srf.srf2llv_py(srf_file, value = 'slip')
    for seg in xrange(len(bounds)):
        # create binary llv file for GMT overlay
        seg_llslips[seg].astype(np.float32).tofile( \
                '%s/PLANES/slip_map_%d.bin' % (out_dir, seg))
        values = np.append(values, seg_llslips[seg][:, -1])
        # also store tinit values retrieved previously
        tinit[seg].astype(np.float32).tofile( \
                '%s/PLANES/tinit_map_%d.bin' % (out_dir, seg))
        # create a mask path for GMT overlay
        geo.path_from_corners(corners = bounds[seg], min_edge_points = 100, \
                output = '%s/PLANES/plane_%d.bounds' % (out_dir, seg))
        # create mask from path
        x_min, y_min = np.min(np_bounds[seg], axis = 0)
        x_max, y_max = np.max(np_bounds[seg], axis = 0)
        seg_regions.append((x_min, x_max, y_min, y_max))
        gmt.grd_mask('%s/PLANES/plane_%d.bounds' % (out_dir, seg), \
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
        cpt_max = round(percentile, - int(floor(log10(abs(percentile)))))
    else:
        # 2 sf
        cpt_max = round(percentile, 1 - int(floor(log10(abs(percentile)))))
    print('Loading complete.')
else:
    bounds = []


###
### OUTPUT 2: corners file for fault plane and hypocentre plot
###
if finite_fault:
    print('Creating corners file...')
    # find hypocentre, use bounds from previous step
    hypocentre = srf.get_hypo(srf_file)
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
else:
    hypocentre = srf.get_hypo(srf_file, depth = True)
    plot_region = (hypocentre[0] - 0.2, hypocentre[0] + 0.2, \
            hypocentre[1] - 0.1, hypocentre[1] + 0.1)
    subfaults = 1
    maximum = srf.srf2llv_py(srf_file, value = 'slip')[0][0][-1]
    percentile = maximum
    average = maximum
    # arbitrary, only changes size of beachball which is relative anyway
    mw = 8
    strike, dip, rake = srf.ps_params(srf_file)
# for plotting region on NZ-wide map
plot_bounds = '%f %f\n%f %f\n%f %f\n%f %f\n' % \
        (plot_region[0], plot_region[2], plot_region[1], plot_region[2], \
        plot_region[1], plot_region[3], plot_region[0], plot_region[3])


###
### OUTPUT 3: GMT MAP
###
print('Plotting SRF on map...')
nz_region = gmt.nz_region
if finite_fault:
    gmt.makecpt(cpt, '%s/slip.cpt' % (out_dir), 0, cpt_max, 1)
gmt.gmt_defaults(wd = out_dir)
# gap on left of maps
gap = 1
# width of NZ map, if changed, other things also need updating
# including tick font size and or tick increment for that map
full_width = 4

### PART A: zoomed in map
p = gmt.GMTPlot('%s/%s_map.ps' % (out_dir, \
        os.path.splitext(os.path.basename(srf_file))[0]))
# cover a larger space, will get cropped at the end
p.background(20, 15)
# this is how high the NZ map will end up being
full_height = gmt.mapproject(nz_region[0], nz_region[3], region = nz_region, \
        projection = 'M%s' % (full_width), wd = out_dir)[1]
# match height of zoomed in map with full size map
zoom_width, zoom_height = \
        gmt.map_width('M', full_height, plot_region, wd = out_dir)
p.spacial('M', plot_region, sizing = zoom_width, \
        x_shift = gap, y_shift = 2.5)
p.basemap(land = 'lightgray', topo_cpt = 'grey1')
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
if finite_fault:
    p.fault(srf_file, is_srf = True, hyp_colour = 'red')
else:
    p.beachballs('%s %s %s %s %s %s %s %s %s\n' \
            % (hypocentre[0], hypocentre[1], hypocentre[2], \
            strike, dip, rake, mw, hypocentre[0], hypocentre[1]), \
            is_file = False, fmt = 'a', scale = 0.4)

p.sites(gmt.sites.keys())
major_tick, minor_tick = gmt.auto_tick(plot_region[0], plot_region[1], \
        zoom_width)
p.ticks(major = '%sd' % (major_tick), \
        minor = '%sd' % (minor_tick), sides = 'ws')

### PART B: NZ map
# draw NZ wide map to show rough location in country
p.spacial('M', nz_region, sizing = full_width, \
    x_shift = zoom_width + gap)
# height of NZ map
full_height = gmt.mapproject(nz_region[0], nz_region[3], wd = out_dir)[1]
p.basemap(land = 'lightgray', topo = gmt.TOPO_LOW, topo_cpt = 'grey1', road = None)
if plot_faults:
    p.path(faults, is_file = True, close = False, width = '0.1p', colour = 'red')
p.path(plot_bounds, is_file = False, close = True, colour = 'black')
# get displacement of box to draw zoom lines later
window_bottom = gmt.mapproject(plot_region[1], plot_region[2], wd = out_dir)
window_top = gmt.mapproject(plot_region[1], plot_region[3], wd = out_dir)
if finite_fault:
    # also draw fault planes / hypocentre
    p.fault(srf_file, is_srf = True, plane_width = '0.4p', top_width = '0.6p', \
            hyp_colour = 'red')
else:
    p.beachballs('%s %s %s %s %s %s %s %s %s\n' \
            % (hypocentre[0], hypocentre[1], hypocentre[2], \
            strike, dip, rake, mw, hypocentre[0], hypocentre[1]), \
            is_file = False, fmt = 'a', scale = 0.05)
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
p.text(total_width / 2.0, total_height, os.path.basename(srf_file), \
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
if finite_fault:
    # scale
    p.cpt_scale(zoom_width / 2.0, -0.5, '%s/slip.cpt' % (out_dir), \
            cpt_max / 4.0, cpt_max / 8.0, \
            label = 'Slip (cm)', length = zoom_width)

p.finalise()
p.png(dpi = dpi, out_dir = srf_dir)
rmtree(out_dir)
print('SRF plot at %s.' % ('%s/srf_map.png' % (srf_dir)))
