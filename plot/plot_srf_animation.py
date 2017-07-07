#!/usr/bin/env python2
"""
Plot animation of SRF file showing risetime and finally slip.
Only for plotting finite faults.
"""

from math import ceil, floor, log10
import multiprocessing as mp
import os
from shutil import copyfile, rmtree
import sys

import numpy as np

import qcore_path
from geo import path_from_corners
import gmt
import srf

# for repository resources
script_dir = os.path.dirname(os.path.abspath(__file__))

try:
    srf_file = sys.argv[1]
    assert(os.path.exists(srf_file))
except IndexError:
    print('Missing SRF filepath parameter.')
    exit(1)
except AssertionError:
    print('SRF file not found: %s' % (srf_file))
    exit(1)
srf_file = os.path.abspath(srf_file)

if len(sys.argv) > 2:
    if sys.argv[2] == 'publish':
        draft = False
    elif sys.argv[2] == 'draft':
        draft = True
    else:
        print('Second parameter invalid. \'publish\' or \'draft\'')
        exit()
else:
    draft = True

# frame output directory
out = os.path.abspath('png_srf_animation')
if not os.path.exists(out):
    os.makedirs(out)
# global base working directory
gwd = os.path.abspath('tmp_srf_animation')
if not os.path.exists(gwd):
    os.makedirs(gwd)

###
### STAGE 0: Metadata
###
# dpi is very important for keeping zoomed aspect ratios within same pixel
# 80 -> 720p, 120 -> 1080p
dpi = 120
# don't change following to keep 16:9 movie ratio
page_width = 16
page_height = 9
# space around map for titles, tick labels and scales etc
margin_top = 0.6
margin_bottom = 0.4
margin_left = 1.0
margin_right = 1.7
map_width = page_width - margin_left - margin_right
map_height = page_height - margin_top - margin_bottom
# number of frames it takes to zoom in
zoom_frames = 50
# srf dt decimation
srf_ddt = 4
# subfault dx, dy assuming flat plane
srf_dx, srf_dy = srf.srf_dxy(srf_file)
plot_dx = '%sk' % (srf_dx * 0.6)
plot_dy = '%sk' % (srf_dy * 0.6)
# bounds for all planes
srf_bounds = srf.get_bounds(srf_file)
srf_x_min, srf_y_min = np.min(np.min(srf_bounds, axis = 0), axis = 0)
srf_x_max, srf_y_max = np.max(np.max(srf_bounds, axis = 0), axis = 0)
# prevent re-calculating SRF corners (speed)
# yet use the file as initial source (accuracy)
srf_corners = '%s/srf_corners.txt' % (gwd)
srf.srf2corners(srf_file, cnrs = srf_corners)
# resources for plane processing
plane_regions = []
for plane in xrange(len(srf_bounds)):
    # create a mask path for GMT overlay
    path_from_corners(corners = srf_bounds[plane], min_edge_points = 100, \
            output = '%s/plane_%d.bounds' % (gwd, plane))
    # create mask from path
    x_min, y_min = np.min(srf_bounds[plane], axis = 0)
    x_max, y_max = np.max(srf_bounds[plane], axis = 0)
    plane_regions.append((x_min, x_max, y_min, y_max))
    gmt.grd_mask('%s/plane_%d.bounds' % (gwd, plane), \
            '%s/plane_%d.mask' % (gwd, plane), \
            dx = plot_dx, dy = plot_dy, region = plane_regions[plane])

# total length of rupture
slip_end = srf.srf2llv_py(srf_file, value = 'ttotal')
rup_time = max([max(slip_end[p][:, 2]) for p in xrange(len(slip_end))])
srf_dt = srf.srf_dt(srf_file)
ftime = srf_ddt * srf_dt
srf_frames = int(ceil(rup_time / ftime))

###
### STAGE 1: Calculate Region Sizing
###
region_nz = gmt.nz_region
region_srf = (srf_x_min, srf_x_max, srf_y_min, srf_y_max)
map_width, map_height, region_srf = \
        gmt.fill_space(map_width, map_height, region_srf, \
                proj = 'M', dpi = dpi, wd = gwd)

###
### STAGE 2: Zoom from NZ Region
###
gmt.gmt_defaults(wd = gwd, ps_media = 'Custom_%dix%dix' % (page_width, page_height))
def zoom_sequence(frame):
    # working directory for this process
    pwd = '%s/zs%.4d' % (gwd, frame)
    if not os.path.exists(pwd):
        os.makedirs(pwd)
    ps_file = '%s/seq_%.4d.ps' % (pwd, frame)
    copyfile('%s/gmt.conf' % (gwd), '%s/gmt.conf' % (pwd))
    z = gmt.GMTPlot(ps_file, reset = False)
    plot_region, plot_width, x_margin, y_margin = gmt.region_transition( \
            'M', region_nz, region_srf, map_width, map_height, \
            dpi, frame, zoom_frames, wd = pwd)

    # extend map to cover margins
    map_width_a, map_height_a, borderless_region = gmt.fill_margins( \
            plot_region, plot_width, dpi, left = margin_left + x_margin, \
            right = margin_right + x_margin, top = margin_top + y_margin, \
            bottom = margin_bottom + y_margin, wd = pwd)
    z.spacial('M', borderless_region, sizing = map_width_a)
    # topo, water, overlay cpt scale
    z.basemap()
    # map margins are semi-transparent
    z.background(map_width_a, map_height_a, \
            colour = 'white@25', spacial = True, \
            window = (margin_left + x_margin, margin_right + x_margin, \
                    margin_top + y_margin, margin_bottom + y_margin))
    # leave space for left tickmarks and bottom colour scale
    z.spacial('M', plot_region, sizing = plot_width, \
            x_shift = margin_left + x_margin, \
            y_shift = margin_bottom + y_margin)
    # title is SRF name
    z.text(sum(plot_region[:2]) / 2., plot_region[3], \
            os.path.basename(srf_file), dy = 0.2, align = 'CB', size = '20p')

    # for scaling fault line thickness
    zfactor = (region_srf[1] - region_srf[0]) \
            / (plot_region[1] - plot_region[0])

    z.fault(srf_corners, is_srf = False, top_width = zfactor * 2.4, \
            plane_width = zfactor * 1.2, hyp_width = zfactor * 1.2)
    tick_major, tick_minor = \
            gmt.auto_tick(plot_region[0], plot_region[1], plot_width)
    z.ticks(major = tick_major, minor = tick_minor)
    z.finalise()
    z.png(dpi = dpi, out_dir = out, clip = False)
    print('Opening zoom sequence %.3d/%.3d complete.' \
            % (frame + 1, zoom_frames))

if not draft:
    #for i in xrange(zoom_frames):
    #    zoom_sequence(i)
    pool = mp.Pool(40)
    pool.map(zoom_sequence, xrange(zoom_frames))

###
### STAGE 3: Slip Animation
###
print('Time series for all subfaults...')
spec_sr = 'sliprate-%s-%s' % (ftime, srf_frames * ftime)
slip_pos, slip_rate = srf.srf2llv_py(srf_file, value = spec_sr)
# produce the max sliprate array for colour palette and stats
sliprate_max = np.array([], dtype = np.float32)
for plane in xrange(len(slip_rate)):
    sliprate_max = np.append(sliprate_max, \
            np.nanmax(slip_rate[plane], axis = 1))
# maximum range to 2 sf
cpt_max = np.nanpercentile(sliprate_max, 50)
cpt_max = round(cpt_max, 1 - int(floor(log10(abs(cpt_max)))))
print('Time series processed.')
# colour palette for slip
cpt_sliprate = '%s/sliprate.cpt' % (gwd)
gmt.makecpt('%s/cpt/slip.cpt' % (script_dir), cpt_sliprate, 0, cpt_max, \
        inc = 2, invert = False)

# create basemap
b = gmt.GMTPlot('%s/basemap.ps' % (gwd), reset = False)
if not draft:
    # extend map to cover margins
    map_width_a, map_height_a, borderless_region = gmt.fill_margins( \
            region_srf, map_width, dpi, left = margin_left, \
            right = margin_right, top = margin_top, bottom = margin_bottom)
    b.spacial('M', borderless_region, sizing = map_width_a)
    # topo, water, overlay cpt scale
    b.basemap()
    # map margins are semi-transparent
    b.background(map_width_a, map_height_a, \
            colour = 'white@25', spacial = True, \
            window = (margin_left, margin_right, margin_top, margin_bottom))
else:
    # background can be larger as whitespace is later cropped
    b.background(page_width, page_height)
# leave space for left tickmarks and bottom colour scale
b.spacial('M', region_srf, sizing = map_width, \
        x_shift = margin_left, y_shift = margin_bottom)
if draft:
    # topo, water, overlay cpt scale
    b.basemap(topo = None, highway = None, road = None, res = 'f')
# title is SRF name
b.text(sum(region_srf[:2]) / 2., region_srf[3], os.path.basename(srf_file), \
        dy = 0.2, align = 'CB', size = '20p')
b.cpt_scale('R', 'B', cpt_sliprate, cpt_max / 10., cpt_max / 50., \
        cross_tick = cpt_max / 10., arrow_f = True, horiz = False, \
        length = map_height - 0.2, thickness = 0.2, \
        dx = 0.2, pos = 'rel_out', align = 'LB', label = 'Slip Rate (cm/s)')
b.leave()
tick_major, tick_minor = \
        gmt.auto_tick(region_srf[0], region_srf[1], map_width)

def slip_sequence(frame):
    # working directory for this process
    pwd = '%s/ss%.4d' % (gwd, frame)
    if not os.path.exists(pwd):
        os.makedirs(pwd)
    ps_file = '%s/seq_%.4d.ps' \
            % (pwd, zoom_frames - zoom_frames * draft + frame)
    copyfile('%s/basemap.ps' % (gwd), ps_file)
    copyfile('%s/gmt.conf' % (gwd), '%s/gmt.conf' % (pwd))
    copyfile('%s/gmt.history' % (gwd), '%s/gmt.history' % (pwd))

    # finish plot
    s = gmt.GMTPlot(ps_file, append = True, reset = False)
    # current time
    s.text(region_srf[1], region_srf[3], 't = %ss' % (frame * ftime), \
            dx = 0.5, dy = 0.2, align = 'LB', size = '16p')
    # prepare overlays
    for plane in xrange(len(srf_bounds)):
        plane_data = np.empty((len(slip_pos[plane]), 3))
        plane_data[:, :2] = slip_pos[plane]
        plane_data[:, 2] = slip_rate[plane][:, frame]
        # do not plot if all values are nan (no data at this time in plane)
        if np.min(np.isnan(plane_data[:, 2])):
            continue
        xyv_file = '%s/plane_%d.bin' % (pwd, plane)
        plane_data.astype(np.float32).tofile(xyv_file)
        gmt.table2grd(xyv_file, '%s.grd' % (xyv_file), \
                grd_type = 'nearneighbor', sectors = 1, min_sectors = 1, \
                region = plane_regions[plane], dx = plot_dx, dy = plot_dy, \
                wd = pwd, search = '%sk' % (srf_dx))
        s.overlay('%s.grd' % (xyv_file), cpt_sliprate, \
                dx = plot_dx, dy = plot_dy, \
                climit = 0.5, crop_grd = '%s/plane_%d.mask' % (gwd, plane), \
                land_crop = False, custom_region = plane_regions[plane], \
                transparency = 0)
    if draft:
        s.coastlines(width = 0.2, res = 'f')
    else:
        s.coastlines(width = 0.2)
    # TODO: fault thickness must be related to zfactor used during zoom
    s.fault(srf_corners, is_srf = False, \
            plane_width = '1', hyp_width = '1', top_width = '2')
    s.ticks(major = tick_major, minor = tick_minor)
    s.finalise()
    s.png(dpi = dpi, out_dir = out, clip = False)
    print('Slip sequence %.3d/%.3d complete.' % (frame + 1, srf_frames))

#for i in xrange(srf_frames):
#    slip_sequence(i)
pool = mp.Pool(26)
pool.map(slip_sequence, xrange(srf_frames))
#pool.map(slip_sequence, xrange(240))

###
### STAGE 4: Show Maximum Sliprate
###
#TODO: future

###
### STAGE 5: Movie and Cleaning
###
gmt.make_movie(os.path.join(out, 'seq_%04d.png'), \
        '%s_animation.mov' % os.path.splitext(srf_file)[0], fps = 20)
rmtree(gwd)
rmtree(out)
