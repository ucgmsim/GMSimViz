#!/usr/bin/env python2

import math
import os
from shutil import copy, rmtree
import sys
from tempfile import mkdtemp

from mpi4py import MPI

import geo
import gmt
import srf

page_width = 16
page_height = 9
# 120 for 1920x1080, 16ix9i
dpi = 600

###
### MASTER
###
if len(sys.argv) > 1:
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('srf_file', help = 'srf file to plot')
    parser.add_argument('-a', '--animate', help = 'create animation', \
        action = 'store_true')
    parser.add_argument('-n', '--nproc', help = 'number of processes to run', \
            type = int, default = int(os.sysconf('SC_NPROCESSORS_ONLN'))
    args = parser.parse_args()
    srf_file = os.path.abspath(sys.argv[1])
    # check file exists
    try:
        assert(os.path.isfile(srf_file))
    except AssertionError:
        print('Could not find SRF: %s' % (srf_file))
        sys.exit(2)

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

    # working directory
    gmt_temp = mkdtemp(prefix = '_GMT_WD_PERSPECTIVE_', \
            dir = os.path.dirname(srf_file))

    

    # cleanup
    rmtree(gmt_temp)
    # stop mpi
    comm.Disconnect()

###
### SLAVE
###
else:
    # connect to parent
    try:
        comm = MPI.Comm.Get_parent()
        rank = comm.Get_rank()
    except:
        print('First parameter must be input file. Parameter not found.')
        print('Alternatively MPI cannot connect to parent.')
        exit(1)

    # ask for work until stop sentinel
    logbook = []
    for task in iter(lambda: comm.sendrecv(None, dest = MASTER), StopIteration):

        t0 = time()
        # no jobs available yet
        if task == None:
            sleep(1)
            logbook.append(('sleep', time() - t0))
        elif task[0] is column_overlay:
            column_overlay(task[1], task[2], task[3], task[4])
            logbook.append(('column overlay', task[1], time() - t0))
        elif task[0] is render:
            render(task[1], task[2], task[3])
            logbook.append(('ps2png', task[3], time() - t0))
        elif task[0] is template_2:
            template_2(task[1])
            logbook.append(('basemap', time() - t0))
        elif task[0] is fault_prep:
            fault_prep(srf_file = task[1], cnrs_file = task[2])
            logbook.append(('fault processing', time() - t0))
        else:
            print('Slave recieved unknown task to complete: %s.' % (task))

    # reports to master
    comm.gather(sendobj = logbook, root = MASTER)
    # shutdown
    comm.Disconnect()


def timeslice(i, n, meta):

    # working directory for current image
    swd = mkdtemp(prefix = '_%.4d_' % (i), dir = meta['wd'])
    copy(os.path.join(meta['wd'], gmt.GMT_CONF), \
            os.path.join(swd, gmt.GMT_CONF))
    copy(os.path.join(meta['wd'], gmt.GMT_HISTORY), \
            os.path.join(swd, gmt.GMT_HISTORY))
    gmt_ps = os.path.join(swd, '%s_perspective%s.ps' \
        % (os.path.splitext(os.path.basename(meta['srf_file']))[0], \
        '_%.4d' % (i) * (n > 1)))

    # smallest size to fill page
    gmt_x_size, gmt_y_size, sx, by = gmt.perspective_fill( \
            page_width, page_height, view = s_azimuth, tilt = map_tilt)
    # TODO: use oblique mercator 'O'
    new_y_size, map_region = gmt.adjust_latitude('M', gmt_x_size, gmt_y_size, \
            map_region, wd = '.', abs_diff = True, accuracy = 0.4 * 1. / dpi, \
            reference = 'left', top = True, bottom = True)

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
    srf_data = gmt.srf2map(srf_file, swd, prefix = 'plane', value = 'slip', \
            cpt_percentile = 95, wd = swd, z = True, xy = True, \
            pz = z_scale * math.cos(math.radians(map_tilt)), dpu = dpi * 0.25)
    p.basemap(topo = None)

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
