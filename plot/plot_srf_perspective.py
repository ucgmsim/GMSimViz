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
map_tilt = max(90 - avg_dip, 10)
# plane domains
bounds = srf.get_bounds(srf_file, depth = True)
top_left = bounds[0][0]
top_right = bounds[-1][1]
top_mid = geo.ll_mid(top_left[0], top_left[1], top_right[0], top_right[1])
gmt_bottom = '\n>\n'.join(['\n'.join([' '.join(map(str, b)) \
        for b in p]) for p in bounds])
gmt_top = '\n>\n'.join(['\n'.join([' '.join(map(str, b)) \
        for b in p[:2]]) for p in bounds])
map_region = (top_mid[0] - 0.6, top_mid[0] + 0.6, top_mid[1] - 0.5, top_mid[1] + 0.4, 0, 1)

# will be function parameters
width = page_width
height = page_height
tilt = map_tilt
view = s_azimuth
# part of the map view (outside) edge is bs, bl, ss, sl
#        /\
#    bl /  \ bs        s|\
#  ___ /____\____      l| \             /|\
#  ss /|GMT |\ sl      y|  \ sl        / | \
# ___/ |PAGE| \___      |___\      bl /  |b \ bs
# sl \ |AREA| / ss      |sx\/sxss    /   |y  \
#  ___\|____|/___      s|  / ss     /____|____\
#      \    /          s| /          blx   bsx
#    bs \  / bl        y|/
#        \/
# repeated values
s_tilt = math.sin(math.radians(map_tilt))
c_tilt = math.cos(math.radians(map_tilt))
s_view = abs(math.sin(math.radians(view)))
c_view = abs(math.cos(math.radians(view)))
gmt_x_angle = math.atan2(s_view * s_tilt, c_view)
gmt_y_angle = math.atan2(s_view, c_view * s_tilt)
# bottom and top edge segments
bs = width * s_view
bsx = bs * s_view
by = math.sqrt(bs ** 2 - bsx ** 2) * s_tilt
bs = math.sqrt(bsx ** 2 + by ** 2)
bl = math.sqrt((width - bsx) ** 2 + by ** 2)
# side segments
ss = height * s_view
sx = ss * c_view
ssy = math.sqrt(ss ** 2 - sx ** 2)
try:
    sx = ssy / math.tan(math.atan((ssy * s_tilt) / sx))
except ZeroDivisionError:
    sx = 0
ss = math.sqrt(ssy ** 2 + sx ** 2)
sl = math.sqrt((height - ssy) ** 2 + sx ** 2)
# result sizes
page_x_size = abs(bl) + abs(ss)
page_y_size = abs(bs) + abs(sl)
# GMT lifts map upwards slightly when map is tilted back
# 'by' is still as before, can only be used for offsets from now
by += c_tilt / 10.
# gmt_x_size and gmt_y_size are pre-tilt dimensions
# with tilt applied, they will be equivalent to page_x_size and page_y_size
gmt_x_size = math.sqrt((page_x_size * math.cos(gmt_x_angle)) ** 2 \
        + (page_x_size * math.sin(gmt_x_angle) / s_tilt) ** 2)
gmt_y_size = math.sqrt((page_y_size * math.sin(gmt_y_angle)) ** 2 \
        + (page_y_size * math.cos(gmt_y_angle) / s_tilt) ** 2)


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

p.spacial('M', map_region, z = 'z-0.1i', sizing = gmt_x_size, \
        p = '%s/%s/0' % (s_azimuth, map_tilt), x_shift = - sx, y_shift = - by)
#p.basemap(topo = None)
p.land(fill = 'darkgreen@80')
p.water(colour = 'black@80')
p.ticks(major = 0.1, minor = 0.01)
p.path(gmt_bottom, is_file = False, colour = 'blue', width = '1p', split = '-', close = True, z = True)
p.path(gmt_top, is_file = False, colour = 'blue', width = '2p', z = True)
#p.sites(gmt.sites_major)

# finish, clean up
p.finalise()
p.png(dpi = dpi, clip = False, out_dir = os.path.dirname(srf_file))
rmtree(gmt_temp)
