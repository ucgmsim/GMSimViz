#!/usr/bin/env python2

import math
import os
from shutil import copy, move, rmtree
import sys
from tempfile import mkdtemp
from time import time

from mpi4py import MPI

import geo
import gmt
import srf

MASTER = 0
TILT_MAX = 20
PAGE_WIDTH = 16
PAGE_HEIGHT = 9
# 120 for 1920x1080, 16ix9i
DPI = 120

def timeslice(i, n, meta):

    # working directory for current image
    swd = mkdtemp(prefix = '_%.4d_' % (i), dir = meta['wd'])
    #copy(os.path.join(meta['wd'], gmt.GMT_CONF), \
    #        os.path.join(swd, gmt.GMT_CONF))
    #copy(os.path.join(meta['wd'], gmt.GMT_HISTORY), \
    #        os.path.join(swd, gmt.GMT_HISTORY))
    gmt_ps = os.path.join(swd, '%s_perspective%s.ps' \
        % (os.path.splitext(os.path.basename(meta['srf_file']))[0], \
        '_%.4d' % (i) * (n > 1)))

    # position in animation
    if meta['s_azimuth'] <= 180:
        azimuth = 180 - (i / float(n)) * (180 - meta['s_azimuth'])
    else:
        azimuth = 180 + (i / float(n)) * (meta['s_azimuth'] - 180)
    tilt = 90 - (i / float(n)) * (90 - meta['map_tilt'])
    # borders on page
    window_t = 0.8
    window_b = 0.3
    scale_p_final = 0.7
    if meta['animate']:
        scale_p = min(i, meta['t_frames']) / meta['t_frames'] * scale_p_final
    else:
        scale_p = scale_p_final

    # smallest size to fill page
    gmt_x_size, gmt_y_size, sx, by = gmt.perspective_fill( \
            PAGE_WIDTH, PAGE_HEIGHT, view = azimuth, tilt = tilt)
    # oblique mercator 'O' would be ideal but doesn't seem to work with 3D
    new_y_size, map_region = gmt.adjust_latitude('M', gmt_x_size, gmt_y_size, \
            meta['map_region'], wd = swd, abs_diff = True, \
            accuracy = 0.4 * 1. / DPI, reference = 'left', top = True, \
            bottom = True)

    p = gmt.GMTPlot(gmt_ps)
    # use custom page size
    gmt.gmt_defaults(wd = swd, \
            ps_media = 'Custom_%six%si' % (PAGE_WIDTH, PAGE_HEIGHT))
    # page background
    p.spacial('X', (0, 1, 0, 1), sizing = '%s/%s' % (PAGE_WIDTH, PAGE_HEIGHT))
    p.path('0 0\n0 1\n1 1\n1 0', is_file = False, close = True, \
            fill = 'white', width = None)
    # geographic projection
    z_scale = -0.1
    p.spacial('M', map_region, z = 'z%s' % (z_scale), sizing = gmt_x_size, \
            p = '%s/%s/0' % (azimuth, tilt), x_shift = - sx, y_shift = - by)
    # load srf plane data
    srf_data = gmt.srf2map(meta['srf_file'], swd, prefix = 'plane', \
            value = 'slip', cpt_percentile = 95, wd = swd, \
            xy = True, pz = z_scale * math.cos(math.radians(tilt)), \
            dpu = DPI)
    p.basemap()
    p.rose('C', 'M', '2i', pos = 'rel', \
            dxp = PAGE_WIDTH / 2 - 2, dyp = PAGE_HEIGHT / 2 - 1.8, \
            fill = 'white@80', clearance = '-0.4i/0.2i', pen = 'thick,red')

    # slip distribution has been reprojected onto x, y of page area
    p.spacial('X', (0, PAGE_WIDTH + sx, 0, PAGE_HEIGHT + by), \
            sizing = '%s/%s' % (PAGE_WIDTH + sx, PAGE_HEIGHT + by))
    for s in xrange(len(srf_data[1])):
        if not os.path.exists('%s/plane_%d_slip_xy.grd' % (swd, s)):
            continue
        p.overlay('%s/plane_%d_slip_xy.grd' % (swd, s), \
                '%s/plane.cpt' % (swd), transparency = 40, \
                crop_grd = '%s/plane_%d_mask_xy.grd' % (swd, s), \
                custom_region = srf_data[1][s])
    p.background(PAGE_WIDTH, PAGE_HEIGHT, spacial = False, \
            window = (0.5, 0.5, window_t, max(window_b, scale_p)), \
            x_margin = sx, y_margin = by, colour = 'white@50')
    p.cpt_scale(PAGE_WIDTH / 2.0 + sx, scale_p + by, '%s/plane.cpt' % (swd), \
            length = PAGE_WIDTH / 1.618, align = 'CT', \
            dy = 0.1, thickness = 0.25, major = srf_data[2][2] / 5., \
            minor = srf_data[2][2] / 20., cross_tick = srf_data[2][2] / 20.)
    # cpt label
    p.text((PAGE_WIDTH - PAGE_WIDTH / 1.618) / 2.0 + sx, scale_p + by, \
            'Slip (cm)', align = 'RM', dx = - 0.2, dy = scale_p_final / -2.0, \
            size = 16)
    # title
    p.text(sx + PAGE_WIDTH / 2.0, by + PAGE_HEIGHT, \
            os.path.basename(meta['srf_file']), align = 'RM', size = 26, \
            dy = window_t / -2.0, dx = - 0.2)
    # slip max, 95%, avg, 25%
    p.text(sx + PAGE_WIDTH / 2.0, by + PAGE_HEIGHT, 'slip max: %s cm' % \
            int(round(srf_data[2][0])), align = 'LB', size = 14, \
            dy = - window_t + 0.45, dx = 0.2)
    p.text(sx + PAGE_WIDTH / 2.0, by + PAGE_HEIGHT, 'slip 95%%: %s cm' % \
            int(round(srf_data[2][1])), align = 'LB', size = 14, \
            dy = - window_t + 0.2, dx = 0.2)

    # place outlines on top of slip distribution
    p.spacial('M', map_region, z = 'z%s' % (z_scale), sizing = gmt_x_size, \
            p = '%s/%s/0' % (azimuth, tilt))
    # srf outline: underground, surface, hypocentre
    p.path(meta['gmt_bottom'], is_file = False, colour = 'blue', width = '1p', \
            split = '-', close = True, z = True)
    p.path(meta['gmt_top'], is_file = False, colour = 'blue', width = '2p', \
            z = True)
    p.points('%s %s %s\n' % (meta['hlon'], meta['hlat'], meta['hdepth']), \
            is_file = False, shape = 'a', size = 0.35, line = 'blue', \
            line_thickness = '1p', z = True, clip = False)
    p.sites(gmt.sites_major)

    # finish, clean up
    p.finalise()
    p.png(dpi = DPI, clip = False, out_dir = meta['wd'])

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
            type = int, default = int(os.sysconf('SC_NPROCESSORS_ONLN')))
    parser.add_argument('-f', '--framerate', help = 'animation framerate', \
            type = int, default = 25)
    parser.add_argument('-t', '--time', help = 'animation transition time (s)', \
            type = float, default = 6.0)
    parser.add_argument('-m', '--mtime', help = 'minor animation transition time (s)', \
            type = float, default = 0.5)
    parser.add_argument('-d', '--delay', help = 'animation start delay (s)', \
            type = float, default = 1.5)
    parser.add_argument('-e', '--end', help = 'animation end delay (s)', \
            type = float, default = 3.0)
    args = parser.parse_args()
    if args.mtime > args.time:
        print('Failed constraints (mtime <= time).')
        sys.exit(1)
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
    map_tilt = max(90 - avg_dip, TILT_MAX)
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

    meta = {'wd':gmt_temp, 'srf_file':srf_file, 's_azimuth':s_azimuth, \
            'map_tilt':map_tilt, 'hlon':hlon, 'hlat':hlat, 'hdepth':hdepth, \
            'gmt_bottom':gmt_bottom, 'gmt_top':gmt_top, 'animate':args.animate, \
            'map_region':map_region, 't_frames':args.mtime * args.framerate}

    # tasks
    if not args.animate:
        msg_list = [(timeslice, 1, 1, meta)]
    else:
        frames = int(args.time * args.framerate)
        msg_list = [(timeslice, i, frames - 1, meta) for i in xrange(frames)]
    nproc = min(len(msg_list), args.nproc)
    # spawn slaves
    comm = MPI.COMM_WORLD.Spawn(
        sys.executable, args = [sys.argv[0]], maxprocs = nproc)
    # job tracking
    in_progress = [None] * nproc
    # distribute work to slaves who ask
    status = MPI.Status()
    while nproc:
        # previous job
        comm.recv(source = MPI.ANY_SOURCE, status = status)
        slave_id = status.Get_source()
        finished = in_progress[slave_id]

        if len(msg_list) == 0:
            # all jobs complete, kill off slaves
            msg_list.append(StopIteration)
            nproc -= 1

        # next job
        msg = msg_list[0]
        del(msg_list[0])
        comm.send(obj = msg, dest = slave_id)
        in_progress[slave_id] = msg

    # gather, print reports from slaves
    reports = comm.gather(None, root = MPI.ROOT)
    # stop mpi
    comm.Disconnect()

    basename = os.path.splitext(os.path.basename(meta['srf_file']))[0]
    if args.animate:
        # copy last frame
        frames_end = int(args.end * args.framerate)
        frames_start = int(args.delay * args.framerate)
        for i in xrange(frames_end):
            copy('%s/%s_perspective_%.4d.png' % (gmt_temp, basename, frames - 1), \
                    '%s/%s_perspective_%.4d.png' \
                        % (gmt_temp, basename, frames + frames_start + i))
        # shift sequence
        for i in xrange(frames - 1, -1, -1):
            move('%s/%s_perspective_%.4d.png' % (gmt_temp, basename, i), \
                    '%s/%s_perspective_%.4d.png' \
                        % (gmt_temp, basename, frames_start + i))
        # copy first frame
        for i in xrange(frames_start):
            copy('%s/%s_perspective_%.4d.png' % (gmt_temp, basename, frames_start), \
                    '%s/%s_perspective_%.4d.png' \
                        % (gmt_temp, basename, i))
        # output movie
        gmt.make_movie('%s/%s_perspective_%%04d.png' % (gmt_temp, basename), \
                basename, fps = args.framerate, codec = 'libx264')
    else:
        move('%s/%s_perspective.png' % (gmt_temp, basename), \
            '%s_perspective.png' % (basename))
    # cleanup
    rmtree(gmt_temp)

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
        elif task[0] is timeslice:
            timeslice(task[1], task[2], task[3])
            logbook.append(('timeslice', task[1], time() - t0))
        else:
            print('Slave recieved unknown task to complete: %s.' % (task))

    # reports to master
    comm.gather(sendobj = logbook, root = MASTER)
    # shutdown
    comm.Disconnect()

