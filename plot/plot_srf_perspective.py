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
s_azimuth = 170
map_tilt = 90
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

bs = (page_width - 6) / math.sin(math.radians(90)) * math.sin(math.radians(s_azimuth))
bl = math.sqrt((page_width - 6) ** 2 - bs ** 2)
#print bs, bl
#bsx = bs / math.sin(math.radians(90)) * abs(math.cos(math.radians(s_azimuth)))
#blx = (page_width - 6) - bsx
#print bsx, blx
## adjusted for tilt and re-calculate bs
#by = math.sqrt(abs(bs ** 2 - bsx ** 2)) * math.sin(math.radians(map_tilt))
#bs = math.sqrt(bsx ** 2 + by ** 2)
#bl = math.sqrt(blx ** 2 + by ** 2)
#print bs, bl
ss = (page_height - 6) / math.sin(math.radians(90)) * math.sin(math.radians(s_azimuth))
sl = math.sqrt((page_height - 6) ** 2 - ss ** 2)
x_size = abs(bl) + abs(ss)
y_size = abs(bs) + abs(sl)
y = abs(bs * math.cos(math.radians(s_azimuth)))
yx = bs * math.sin(math.radians(s_azimuth))
x = abs(ss * math.cos(math.radians(s_azimuth)))
xy = ss * math.sin(math.radians(s_azimuth))

sdiff = 1. / math.sin(math.radians(map_tilt)) - 1
#x_size += x_size * sdiff * math.sin(math.radians(s_azimuth % 180))
#y_size += y_size * sdiff * math.cos(math.radians(s_azimuth % 180))

new_y_size, map_region = gmt.adjust_latitude('M', x_size, y_size, map_region, \
        wd = '.', abs_diff = True, accuracy = 0.4 * 1. / dpi, reference = 'left', \
        top = True, bottom = True)
print x_size, y_size, new_y_size

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
p.path('0 0\n0 1\n1 1\n1 0', is_file = False, close = True)

p.spacial('M', map_region, z = 'z-0.1i', sizing = x_size, \
        p = '%s/%s/0' % (s_azimuth, map_tilt), x_shift = - x, y_shift = - y)
#p.basemap(topo = None)
p.land(fill = 'darkgreen@50')
p.water(colour = 'lightblue@50')
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
p.path('%s %s\n3 3\n%s %s\n%s %s' % (3 - x, 3 + (9 - 6) - xy, 3 + yx, 3 - y, 3 + 16 - 6, 3), is_file = False, close = False)

# finish, clean up
p.finalise()
p.png(dpi = dpi, clip = False, out_dir = os.path.dirname(srf_file))
rmtree(gmt_temp)
