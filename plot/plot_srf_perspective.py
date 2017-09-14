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
dpi = 120

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
s_azimuth = 190
map_tilt = 50
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


# part of the map view (outside) edge is bs, bl, ss, sl
#        /\
#    bl /  \ bs        s|\
#  ___ /____\____      l| \             /|\
#  ss /|GMT |\ sl      y|  \ sl        / | \
# ___/ |PAGE| \___      |___\      bl /  |b \ bs
# sl \ |AREA| / ss      |sx /        /   |y  \
#  ___\|____|/___      s|  / ss     /____|____\
#      \    /          s| /          blx   bly
#    bs \  / bl        y|/
#        \/


bs = (page_width - 6) / math.sin(math.radians(90)) * math.sin(math.radians(s_azimuth))
bsx = bs / math.sin(math.radians(90)) * abs(math.sin(math.radians(s_azimuth)))
blx = (page_width - 6) - bsx
# adjusted for tilt and re-calculate bs
by = math.sqrt(abs(bs ** 2 - bsx ** 2)) * math.sin(math.radians(map_tilt))
bs = math.sqrt(bsx ** 2 + by ** 2)
bl = math.sqrt(blx ** 2 + by ** 2)
ss = (page_height - 6) / math.sin(math.radians(90)) * math.sin(math.radians(s_azimuth))
sl = math.sqrt((page_height - 6) ** 2 - ss ** 2)
ssy = ss / math.sin(math.radians(90)) * abs(math.sin(math.radians(s_azimuth)))
sx = math.sqrt(abs(ss ** 2 - ssy ** 2))
ssy *= math.sin(math.radians(map_tilt))
sly = (page_height - 6) - ssy
ss = math.sqrt(ssy ** 2 + sx ** 2)
sl = math.sqrt(sly ** 2 + sx ** 2)

x_size = abs(bl) + abs(ss)
y_size = abs(bs) + abs(sl)
y = abs(bs * math.cos(math.radians(s_azimuth)))
yx = bs * math.sin(math.radians(s_azimuth))
x = abs(ss * math.cos(math.radians(s_azimuth)))
xy = ss * math.sin(math.radians(s_azimuth))
print x_size, y_size
sdiff = 1. / math.sin(math.radians(map_tilt)) - 1
print sdiff
x_size += x_size * sdiff * abs(math.sin(math.radians(s_azimuth % 180))) * - math.cos(math.radians(map_tilt))
y_size += y_size * sdiff * abs(math.cos(math.radians(s_azimuth % 180))) * math.sin(math.radians(map_tilt))
# GMT lifts map upwards slightly when map is tilted back
y += (ssy + sly) * math.cos(math.radians(map_tilt)) * math.sin(math.radians(20)) * 0.1

new_y_size, map_region = gmt.adjust_latitude('M', x_size, y_size, map_region, \
        wd = '.', abs_diff = True, accuracy = 0.4 * 1. / dpi, reference = 'left', \
        top = True, bottom = True)
print x_size, y_size, new_y_size
#exit()

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

p.spacial('X', (0, 1, 0, 1), sizing = '%s/%s' % (page_width - 6, page_height - 6), \
        x_shift = 3, y_shift = 3)
p.path('0 0\n0 1\n1 1\n1 0', is_file = False, close = True, colour = 'red')

p.spacial('M', map_region, z = 'z-0.1i', sizing = x_size, \
        p = '%s/%s/0' % (s_azimuth, map_tilt), x_shift = - x, y_shift = - y)
#p.basemap(topo = None)
p.land(fill = 'darkgreen@80')
p.water(colour = 'lightblue@80')
p.ticks(major = 0.1, minor = 0.01)
p.path(gmt_bottom, is_file = False, colour = 'blue', width = '1p', split = '-', close = True, z = True)
p.path(gmt_top, is_file = False, colour = 'blue', width = '2p', z = True)
#p.sites(gmt.sites_major)
print gmt.mapproject(map_region[1], map_region[3], wd = gmt_temp, projection = None, region = None, \
    inverse = False, unit = None, z = '-Jz1', p = True)
print gmt.mapproject(map_region[0], map_region[2], wd = gmt_temp, projection = None, region = None, \
    inverse = False, unit = None, z = '-Jz1', p = True)
print gmt.mapproject(map_region[0], map_region[3], wd = gmt_temp, projection = None, region = None, \
    inverse = False, unit = None, z = '-Jz1', p = True)
print gmt.mapproject(map_region[1], map_region[2], wd = gmt_temp, projection = None, region = None, \
    inverse = False, unit = None, z = '-Jz1', p = True)

p.spacial('X', (0, page_width, 0, page_height), sizing = '%s/%s' % (page_width, page_height), \
        x_shift = x - 3, y_shift = y - 3)
#p.path('%s %s\n3 3\n%s %s\n%s %s' % (3 - x, 3 + (9 - 6) - xy, 3 + yx, 3 - y, 3 + 16 - 6, 3), is_file = False, close = False)

# finish, clean up
p.finalise()
p.png(dpi = dpi, clip = False, out_dir = os.path.dirname(srf_file))
rmtree(gmt_temp)
