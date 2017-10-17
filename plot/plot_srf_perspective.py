#!/usr/bin/env python2

import math
import os
from shutil import rmtree
import sys
from tempfile import mkdtemp

import geo
import gmt
import srf

page_width = 16
page_height = 9
# 120 for 1920x1080, 16ix9i
dpi = 600

# retrieve srf file parameter
if len(sys.argv) > 1:
    srf_file = os.path.abspath(sys.argv[1])
else:
    print('First parameter is the path of the SRF file.')
    sys.exit(1)
# check file exists
try:
    assert(os.path.isfile(srf_file))
except AssertionError:
    print('Could not find SRF: %s' % (srf_file))
    sys.exit(1)
# load plane data
try:
    planes = srf.read_header(srf_file, idx = True)
except (ValueError, IndexError):
    print('Failed to read SRF: %s' % (srf_file))
    sys.exit(1)
# information from plane data
avg_strike = geo.avg_wbearing([(p['strike'], p['length']) for p in planes])
avg_dip = planes[0]['dip']
s_azimuth = avg_strike + 90
map_tilt = max(90 - avg_dip, 20)
# plane domains
bounds = srf.get_bounds(srf_file, depth = True)
hlon, hlat, hdepth = srf.get_hypo(srf_file, depth = True)
top_left = bounds[0][0]
top_right = bounds[-1][1]
top_mid = geo.ll_mid(top_left[0], top_left[1], top_right[0], top_right[1])
gmt_bottom = '\n>\n'.join(['\n'.join([' '.join(map(str, b)) \
        for b in p]) for p in bounds])
gmt_top = '\n>\n'.join(['\n'.join([' '.join(map(str, b)) \
        for b in p[:2]]) for p in bounds])
map_region = (top_mid[0] - 2, top_mid[0] + 2, top_mid[1] - 0.5, top_mid[1] + 0.4, -8, 0)
# smallest size to fill page
gmt_x_size, gmt_y_size, sx, by = gmt.perspective_fill(page_width, page_height, \
        view = s_azimuth, tilt = map_tilt, zlevel = 0)


# TODO: use oblique mercator 'O'
new_y_size, map_region = gmt.adjust_latitude('M', gmt_x_size, gmt_y_size, map_region, \
        wd = '.', abs_diff = True, accuracy = 0.4 * 1. / dpi, reference = 'left', \
        top = True, bottom = True)

# plotting resources
gmt_temp = mkdtemp(prefix = 'GMT_WD_PERSPECTIVE_', \
        dir = os.path.dirname(srf_file))
gmt_ps = os.path.join(gmt_temp, '%s_perspective.ps' \
        % (os.path.splitext(os.path.basename(srf_file))[0]))
p = gmt.GMTPlot(gmt_ps)
# use custom page size
gmt.gmt_defaults(wd = gmt_temp, \
        ps_media = 'Custom_%six%si' % (page_width, page_height))
# page background
p.spacial('X', (0, 1, 0, 1), sizing = '%s/%s' % (page_width, page_height))
p.path('0 0\n0 1\n1 1\n1 0', is_file = False, close = True, \
        fill = 'white', width = None)
# geographic projection
z_scale = -0.1
p.spacial('M', map_region, z = 'z%s' % (z_scale), sizing = gmt_x_size, \
        p = '%s/%s/0' % (s_azimuth, map_tilt), x_shift = - sx, y_shift = - by)
# load srf plane data
srf_data = gmt.srf2map(srf_file, gmt_temp, prefix = 'plane', value = 'slip', \
        cpt_percentile = 95, wd = gmt_temp, z = True, xy = True, \
        pz = z_scale * math.cos(math.radians(map_tilt)), dpu = dpi * 0.25)
p.basemap(topo = None)
p.ticks(major = 0.1, minor = 0.01)

# slip distribution has been reprojected onto x, y of page area
p.spacial('X', (0, page_width + sx, 0, page_height + by), \
        sizing = '%s/%s' % (page_width + sx, page_height + by))
for s in xrange(len(srf_data[3])):
    p.overlay('%s/plane_%d_slip_xy.grd' % (gmt_temp, s), \
            '%s/plane.cpt' % (gmt_temp), transparency = 40, \
            crop_grd = '%s/plane_%d_mask_xy.grd' % (gmt_temp, s))
p.background(page_width, page_height, spacial = False, \
        window = (0.5, 0.5, 0.8, 0.3), x_margin = sx, y_margin = by, \
        colour = 'white@50')
p.text(sx + page_width / 2.0, by + page_height - 0.3, \
        os.path.basename(srf_file), align = 'CT', size = 26)

# place outlines on top of slip distribution
p.spacial('M', map_region, z = 'z%s' % (z_scale), sizing = gmt_x_size, \
        p = '%s/%s/0' % (s_azimuth, map_tilt))
# srf outline: underground, surface, hypocentre
p.path(gmt_bottom, is_file = False, colour = 'blue', width = '1p', \
        split = '-', close = True, z = True)
p.path(gmt_top, is_file = False, colour = 'blue', width = '2p', z = True)
p.points('%s %s %s\n' % (hlon, hlat, hdepth), is_file = False, shape = 'a', \
        size = 0.35, line = 'blue', line_thickness = '1p', \
        z = True, clip = False)
p.sites(gmt.sites_major)

# finish, clean up
p.finalise()
p.png(dpi = dpi, clip = False, out_dir = os.path.dirname(srf_file))
rmtree(gmt_temp)
