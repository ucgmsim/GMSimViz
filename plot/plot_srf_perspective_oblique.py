#!/usr/bin/env python2

import math
import os
from shutil import copy, move, rmtree
import sys
from tempfile import mkdtemp
from time import time

from mpi4py import MPI
import numpy as np

import geo
import gmt
import srf
import xyts

MASTER = 0
TILT_MAX = 20
PAGE_WIDTH = 16
PAGE_HEIGHT = 9
SCALE_WIDTH = PAGE_WIDTH / 1.618
SCALE_SIZE = 0.25
SCALE_PAD = 0.1
OVERLAY_T = 40
# 100 for 1600x900
# 120 for 1920x1080, 16ix9i
# 240 for 3840x2160
DPI = 100
# borders on page
WINDOW_T = 0.8
WINDOW_B = 0.3
WINDOW_L = 0.5
WINDOW_R = 0.5

def load_xyts(meta):
    """
    Complete time-intensive tasks which aren't needed for beginning frames.
    """
    xfile = xyts.XYTSFile(meta['xyts_file'], meta_only = False)
    # higher res outline needed for mask to follow great circle line
    # not too high: becomes time intensive
    geo.path_from_corners(corners = xfile.corners(), min_edge_points = 15, \
            output = '%s/xyts/corners-hr.gmt' % (meta['wd']))
    gmt.grd_mask('%s/xyts/corners-hr.gmt' % (meta['wd']), \
            '%s/xyts/mask.nc' % (meta['wd']), \
            dx = meta['xyts_res'], dy = meta['xyts_res'], \
            region = meta['xyts_region'], wd = '%s/xyts' % (meta['wd']))
    # pgv used to generate cpt scale
    xfile.pgv(pgvout = '%s/xyts/pgv.bin' % (meta['wd']))
    cpt_max = gmt.xyv_cpt_range('%s/xyts/pgv.bin' % (meta['wd']))[2]
    gmt.makecpt('magma', '%s/xyts/gm.cpt' % (meta['wd']), 0, cpt_max, \
            invert = True, wd = '%s/xyts' % (meta['wd']), continuing = True)

def load_xyts_ts(meta, job):
    """
    Prepare xyts timeslice overlays.
    dependencies: xyts_cpt_max is available
    """
    # prevent gmt.conf clashes with unique working directory
    tmp = os.path.join(meta['wd'], 'xyts', '_%03d_' % (job['start']))
    os.makedirs(tmp)

    xfile = xyts.XYTSFile(meta['xyts_file'], meta_only = False)
    crop_grd = os.path.join(meta['wd'], 'xyts', 'mask.nc')
    # preload timeslice overlays
    for t in xrange(job['start'], xfile.nt, job['inc']):
        ts_prefix = os.path.join(meta['wd'], 'xyts', 'ts%04d' % (t))
        # load binary data
        xfile.tslice_get(t, outfile = '%s.bin' % (ts_prefix))
        # store as netCDF
        gmt.table2grd('%s.bin' % (ts_prefix), '%s.nc' % (ts_prefix), \
                grd_type = 'surface', region = meta['xyts_region'], \
                dx = meta['xyts_res'], climit = meta['xyts_cpt_max'] * 0.01, \
                wd = tmp, tension = '1.0')
        os.remove('%s.bin' % (ts_prefix))
        # crop values outside sim domain
        rc = gmt.grdmath(['%s.nc' % (ts_prefix), crop_grd, 'MUL', \
                '=', '%s.nc' % (ts_prefix)], wd = tmp)
        # don't show insignificant ground motion
        rc = gmt.grdclip('%s.nc' % (ts_prefix), '%s.nc' % (ts_prefix), \
                min_v = meta['xyts_cpt_max'] * 0.03, wd = tmp)
        # nothing to display
        if rc == gmt.STATUS_INVALID:
            os.remove('%s.nc' % (ts_prefix))

    # clean up
    rmtree(tmp)

def timeslice(job, meta):
    """
    Render image in animation.
    """

    if job['seq'] == None:
        i = 0
    else:
        i = job['seq']
    # working directory for current image
    swd = os.path.join(meta['wd'], '_%.4d_' % (i))
    os.makedirs(swd)
    gmt_ps = os.path.join(swd, '%s_perspective%s.ps' \
        % (os.path.splitext(os.path.basename(meta['srf_file']))[0], \
        '_%.4d' % (i) * (job['seq'] != None)))

    # 
    map_region = meta['map_region']
    map_region2 = (-50, 50, -30, 30)
    lon0 = sum(map_region[:2]) / 2.0
    lat0 = sum(map_region[2:4]) / 2.0
    projection = 'OA%s/%s/%s/%s' % (lon0, lat0, job['azimuth'] - 90, PAGE_WIDTH)
    map_region2 = gmt.fill_space_oblique(lon0, lat0, PAGE_WIDTH, \
            PAGE_HEIGHT / math.sin(math.radians(job['tilt'])), \
            map_region2, 'k', projection, DPI, swd)
    corners, llur = gmt.map_corners(projection = projection, \
            region = map_region2, region_units = 'k', return_region = 'llur', \
            wd = swd)

    p = gmt.GMTPlot(gmt_ps)
    # use custom page size
    gmt.gmt_defaults(wd = swd, \
            ps_media = 'Custom_%six%si' % (PAGE_WIDTH, PAGE_HEIGHT))
    z_scale = -0.1
    def proj(projected):
        if projected:
            p.spacial('OA', llur, lon0 = lon0, lat0 = lat0, \
                    z = 'z%s' % (z_scale), \
                    sizing = '%s/%s' % (job['azimuth'] - 90, PAGE_WIDTH), \
                    p = '180/%s/0' % (job['tilt']))
        else:
            p.spacial('X', (0, PAGE_WIDTH, 0, PAGE_HEIGHT), \
                    sizing = '%s/%s' % (PAGE_WIDTH, PAGE_HEIGHT))

    proj(True)
    p.basemap(topo = None, road = None)

    # simulation domain
    if os.path.isfile('%s/xyts/corners.gmt' % (meta['wd'])):
        p.path('%s/xyts/corners.gmt' % (meta['wd']), close = True, \
                width = '2p', split = '-', colour = '60/60/60')
    p.sites(gmt.sites_major)
    p.rose('C', 'M', '1.8i', pos = 'rel', dxp = PAGE_WIDTH / 2.0 - 1.8, \
            dyp = PAGE_HEIGHT / 2.0 - 2.2, \
            fill = 'white@80', clearance = '0.2i', pen = 'thick,red')

    # srf outline: underground, surface, hypocentre
    p.path(meta['gmt_bottom'], is_file = False, colour = 'black@30', width = '1p', \
            split = '-', close = True, z = True)
    p.path(meta['gmt_top'], is_file = False, colour = 'black', width = '2p', \
            z = True)
    p.points('%s %s %s\n' % (meta['hlon'], meta['hlat'], meta['hdepth']), \
            is_file = False, shape = 'a', size = 0.35, line = 'black', \
            line_thickness = '1p', z = True, clip = False)
    # load srf plane data
    if job['sim_time'] < 0:
        plot = 'slip'
        scale_p_final = 0.7
        srf_data = gmt.srf2map(meta['srf_file'], swd, prefix = 'plane', \
                value = 'slip', cpt_percentile = 95, wd = swd, \
                xy = True, pz = z_scale * math.cos(math.radians(job['tilt'])), \
                dpu = DPI)
    elif job['sim_time'] >= 0:
        plot = 'timeseries'
        scale_p_final = 1.0
        regions_sr = []
        # TODO: major refactoring, already exists within gmt.srf2map
        srt = int(round(job['sim_time'] / meta['srf_dt']))
        if srt >= meta['sr_len']:
            srt = meta['sr_len'] - 1
        for i in xrange(meta['n_plane']):
            # lon, lat, depth
            subfaults = np.fromfile('%s/subfaults_%d.bin' % (meta['wd'], i), \
                    dtype = '3f')
            # reproject
            xyv = np.empty((subfaults.shape[0], 3))
            xyv[:, :2] = gmt.mapproject_multi(subfaults[:, :2], wd = swd, \
                    z = '-Jz%s' % (z_scale), p = True)
            xyv[:, 1] += subfaults[:, 2] \
                    * z_scale * math.cos(math.radians(job['tilt']))
            xyv[:, 2] = np.fromfile('%s/sliptss_%d.bin' % (meta['wd'], i), \
                    dtype = 'f').reshape(len(subfaults), -1)[:, srt]
            # region
            x_min, y_min = np.min(xyv[:, :2], axis = 0)
            x_max, y_max = np.max(xyv[:, :2], axis = 0)
            regions_sr.append((x_min, x_max, y_min, y_max))
            # dump as binary
            xyv.astype(np.float32) \
                .tofile('%s/slip_%d.bin' % (swd, i))
            # search radius based on diagonal distance
            p2 = xyv[meta['planes'][i]['nstrike'], :2]
            search = math.sqrt(abs(xyv[0, 0] - p2[0]) ** 2 \
                    + abs(xyv[0, 1] - p2[1]) ** 2) * 1.1
            rc = gmt.table2grd('%s/slip_%d.bin' % (swd, i), \
                    '%s/slip_%d.grd' % (swd, i), \
                    file_input = True, grd_type = 'nearneighbor', \
                    region = regions_sr[i], dx = 1.0 / DPI, dy = 1.0 / DPI, \
                    wd = swd, geo = False, search = search, min_sectors = 2)
            if rc == gmt.STATUS_INVALID \
                    and os.path.exists('%s/slip_%d.grd' % (swd, i)):
                os.remove('%s/slip_%d.grd' % (swd, i))

    # slip distribution has been reprojected onto x, y of page area
    proj(False)
    if plot == 'slip':
        for s in xrange(len(srf_data[1])):
            if not os.path.exists('%s/plane_%d_slip_xy.grd' % (swd, s)):
                continue
            p.overlay('%s/plane_%d_slip_xy.grd' % (swd, s), \
                    '%s/slip.cpt' % (meta['srf_wd']), \
                    transparency = job['transparency'], \
                    crop_grd = '%s/plane_%d_mask_xy.grd' % (swd, s))
                    #custom_region = srf_data[1][s])
    elif plot == 'timeseries':
        for s in xrange(meta['n_plane']):
            if not os.path.exists('%s/slip_%d.grd' % (swd, s)):
                continue
            p.overlay('%s/slip_%d.grd' % (swd, s), \
                    '%s/slip.cpt' % (meta['srf_wd']), \
                    transparency = job['transparency'], \
                    #crop_grd = '%s/plane_%d_mask_xy.grd' % (swd, s), \
                    custom_region = regions_sr[s])
    # ground motion must be in geographic projection due to 3D Z scaling
    proj(True)
    try:
        xpos = int(round(job['sim_time'] / meta['xyts_dt']))
    except KeyError:
        # xyts file has not been given
        xpos = -1
    gm_file = os.path.join(meta['wd'], 'xyts', 'ts%04d.nc' % (xpos))
    if os.path.isfile(gm_file):
        p.overlay3d(gm_file, cpt = '%s/xyts/gm.cpt' % (meta['wd']), \
                transparency = job['transparency'], dpi = DPI, \
                z = '-Jz%s' % (1.5 / meta['xyts_cpt_max']), \
                mesh = True, mesh_pen = '0.1p')

    # calculate inner region for map ticks
    # manually adjusted y as mapproject -I not compatible with -p
    scale_p = job['scale_t'] * scale_p_final
    window_b = max(WINDOW_B, scale_p)
    # also include lr for better auto tick increment calculation
    llur_i = gmt.mapproject_multi([ \
            [WINDOW_L, window_b / math.sin(math.radians(job['tilt']))], \
            [PAGE_WIDTH - WINDOW_R, \
            (PAGE_HEIGHT - WINDOW_T) / math.sin(math.radians(job['tilt']))], \
            [PAGE_WIDTH - WINDOW_R, \
            window_b / math.sin(math.radians(job['tilt']))]], \
            inverse = True, wd = swd)

    ###
    ### map overlay border, labels, legends etc...
    ###
    proj(False)
    p.background(PAGE_WIDTH, PAGE_HEIGHT, spacial = False, \
            window = (WINDOW_L, WINDOW_R, WINDOW_T, window_b), \
            colour = 'white@50')
    # middle of scale
    cpt_y = scale_p - 0.5 * SCALE_SIZE - SCALE_PAD
    # space before scale starts
    scale_margin = (PAGE_WIDTH - SCALE_WIDTH) / 2.0
    if plot == 'slip':
        cpt_label = 'Slip (cm)'
        p.cpt_scale(PAGE_WIDTH / 2.0, scale_p, \
                '%s/slip.cpt' % (meta['srf_wd']), length = SCALE_WIDTH, \
                align = 'CT', dy = SCALE_PAD, thickness = SCALE_SIZE, \
                major = meta['slip_cpt_max'] / 5., \
                minor = meta['slip_cpt_max'] / 20., \
                cross_tick = meta['slip_cpt_max'] / 20.)
    elif plot == 'timeseries':
        cpt_label = ''
        y = scale_p
        if meta['xyts_file'] != None:
            x0 = WINDOW_L * 2.0
            x1 = scale_margin
            diff = x1 - x0
            x = x0 + diff * job['scale_x']
            length0 = PAGE_WIDTH / 2.0 - WINDOW_L * 3
            length1 = SCALE_WIDTH
            diff = length1 - length0
            length = length0 + diff * job['scale_x']
            p.cpt_scale(x, y, \
                '%s/xyts/gm.cpt' % (meta['wd']), \
                length = length, align = 'LT', dy = SCALE_PAD, \
                thickness = SCALE_SIZE, major = meta['xyts_cpt_max'] / 5., \
                minor = meta['xyts_cpt_max'] / 20., \
                cross_tick = meta['xyts_cpt_max'] / 20., label = 'Ground motion (cm/s)')
            x += length + WINDOW_L * 2 + scale_margin * job['scale_x']
        else:
            x = scale_margin
            length = SCALE_WIDTH
        p.cpt_scale(x, y, \
                '%s/slip.cpt' % (meta['srf_wd']), length = length, \
                align = 'LT', dy = SCALE_PAD, thickness = SCALE_SIZE, \
                major = meta['slip_cpt_max'] / 5., \
                minor = meta['slip_cpt_max'] / 20., \
                cross_tick = meta['slip_cpt_max'] / 20., \
                label = 'Slip (cm)')
    # cpt label
    if cpt_label != '':
        p.text(scale_margin, cpt_y, \
                cpt_label, align = 'RM', dx = - SCALE_PAD, size = 16)
    # title
    p.text(PAGE_WIDTH / 2.0, PAGE_HEIGHT, \
            os.path.basename(meta['srf_file']), align = 'RM', size = 26, \
            dy = WINDOW_T / -2.0, dx = - 0.2)
    # sim time
    if job['sim_time'] >= 0:
        p.text(PAGE_WIDTH - WINDOW_R, PAGE_HEIGHT - WINDOW_T, \
                '%.3fs' % (job['sim_time']), size = '24p', \
                align = "BR", font = 'Courier', dx = -0.2, dy = 0.1)

    if plot == 'slip':
        # box-and-whisker slip distribution
        scale_start = scale_margin
        scale_factor = 1.0 / meta['slip_cpt_max'] * SCALE_WIDTH
        # max point should not be off the page, leave space for label
        max_x = scale_start \
                + min(srf_data[2]['max'] * scale_factor, SCALE_WIDTH + 1.0)
        p.epoints('%s %s %s %s %s %s' \
                % (scale_start + srf_data[2]['50p'] * scale_factor, cpt_y, \
                scale_start + srf_data[2]['min'] * scale_factor, \
                scale_start + srf_data[2]['25p'] * scale_factor, \
                scale_start + srf_data[2]['75p'] * scale_factor, max_x), \
                is_file = False, xy = 'X', asymmetric = True, \
                width = SCALE_SIZE, colour = 'blue', line_width = '2p')
        # label max
        p.text(max_x, cpt_y, '%.1f' % (srf_data[2]['max']), size = '16p', \
                align = 'LM', dx = SCALE_PAD)

    # final projection of the inner area to draw map ticks
    p.spacial('OA%s/%s/%s/' % (lon0, lat0, job['azimuth'] - 90), \
            (str(llur_i[0][0]), str(llur_i[0][1]), \
            str(llur_i[1][0]), '%sr' % (llur_i[1][1])), \
            sizing = PAGE_WIDTH - WINDOW_L - WINDOW_R, \
            p = '180/%s/0' % (job['tilt']), \
            x_shift = WINDOW_L, y_shift = window_b)
    # MAP_FRAME_TYPE also changes cpt scales so cannot be globally set
    gmt.gmt_set(['MAP_FRAME_TYPE', 'inside', 'MAP_TICK_LENGTH_PRIMARY', '0.1i', 'MAP_ANNOT_OBLIQUE', '1'], wd = swd)
    p.ticks(major = '1d', minor = '0.2d')

    # finish, clean up
    p.finalise()
    p.png(dpi = DPI, clip = False, out_dir = meta['wd'])
    # temporary storage can get very large
    rmtree(swd)

###
### MASTER
###
if len(sys.argv) > 1:
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('srf_file', help = 'srf file to plot')
    parser.add_argument('-x', '--xyts', help = 'xyts file to plot')
    parser.add_argument('-a', '--animate', help = 'create animation', \
            action = 'store_true')
    parser.add_argument('-n', '--nproc', help = 'number of processes to run', \
            type = int, default = int(os.sysconf('SC_NPROCESSORS_ONLN')))
    parser.add_argument('-f', '--framerate', help = 'animation framerate', \
            type = int, default = 30)
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
    if args.framerate < 5:
        print('Framerate too low: %s' % (args.framerate))
        sys.exit(1)
    nproc = args.nproc
    # srf check
    srf_file = os.path.abspath(args.srf_file)
    try:
        assert(os.path.isfile(srf_file))
    except AssertionError:
        print('Could not find SRF: %s' % (srf_file))
        sys.exit(2)
    # xyts optional
    if args.xyts != None:
        xyts_file = os.path.abspath(args.xyts)
        if not os.path.isfile(xyts_file):
            print('Could not find XYTS: %s' % (xyts_file))
            sys.exit(2)
    else:
        xyts_file = None

    # list of work for slaves to do
    msg_list = []
    # dependencies for tasks yet to be added
    msg_deps = 0

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
    #print np.min(bounds, axis = 0)
    #exit()
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
            'map_region':map_region, 't_frames':int(args.mtime * args.framerate)}

    # TODO: sliprate preparation should be an early task
    slip_end = srf.srf2llv_py(srf_file, value = 'ttotal')
    rup_time = max([max(slip_end[p][:, 2]) for p in xrange(len(slip_end))])
    # internal dt
    srf_dt = srf.srf_dt(srf_file)
    # frames per slip rate increment
    fpdt = 1 / (srf_dt * args.framerate)
    # decimation of srf slip rate dt to show
    srf_ddt = max(1, math.floor(fpdt))
    # desimated dt
    ddt = srf_ddt * srf_dt
    ts_sr = int(math.ceil(rup_time / ddt))
    time_sr = ts_sr * ddt
    # frames containing slip rates
    frames_sr = int(time_sr * args.framerate)
    # TODO: possibly interpolate in future
    spec_sr = 'slipts-%s-%s' % (ddt, time_sr)
    slip_pos, slip_rate = srf.srf2llv_py(srf_file, value = spec_sr, depth = True)
    meta['planes'] = srf.read_header(srf_file, idx = True)
    for plane in xrange(len(slip_pos)):
        slip_pos[plane].astype(np.float32).tofile( \
                os.path.join(gmt_temp, 'subfaults_%d.bin' % (plane)))
        slip_rate[plane].astype(np.float32).tofile( \
                os.path.join(gmt_temp, 'sliptss_%d.bin' % (plane)))
    meta['srf_dt'] = ddt
    meta['sr_len'] = ts_sr
    meta['n_plane'] = len(slip_pos)
    # prepare cpt
    meta['srf_wd'] = os.path.join(gmt_temp, 'srf')
    os.makedirs(meta['srf_wd'])
    seg_slips = srf.srf2llv_py(srf_file, value = 'slip')
    all_vs = np.concatenate((seg_slips))[:, -1]
    percentile = np.percentile(all_vs, 95)
    # round percentile significant digits for colour pallete
    if percentile < 1000:
        # 1 sf
        cpt_max = round(percentile, \
                - int(math.floor(math.log10(abs(percentile)))))
    else:
        # 2 sf
        cpt_max = round(percentile, \
                1 - int(math.floor(math.log10(abs(percentile)))))
    meta['slip_cpt_max'] = cpt_max
    gmt.makecpt(gmt.CPTS['slip'], '%s/%s.cpt' % (meta['srf_wd'], 'slip'), 0, \
            cpt_max, max(1, cpt_max / 100))

    meta['xyts_file'] = xyts_file
    frames_gm = 0
    # xyts quick preparation
    if xyts_file != None:
        os.makedirs(os.path.join(gmt_temp, 'xyts'))
        xfile = xyts.XYTSFile(xyts_file, meta_only = True)
        xcnrs = xfile.corners(gmt_format = True)
        xregion = xfile.region(corners = xcnrs[0])
        # TODO: just use dx?
        xres = '%sk' % (xfile.hh * xfile.dxts * 3.0 / 5.0)
        with open('%s/xyts/corners.gmt' % (gmt_temp), 'w') as xpath:
            xpath.write(xcnrs[1])
        meta['xyts_region'] = xregion
        meta['xyts_res'] = xres
        meta['xyts_dt'] = xfile.dt
        if args.animate:
            msg_list.append([load_xyts, meta])
            msg_deps += 1
            frames_gm = int(xfile.dt * (xfile.nt - 0.6) * args.framerate)

    # tasks
    if not args.animate:
        msg_list = [(timeslice, {'azimuth':s_azimuth, 'tilt':map_tilt, \
                'scale_t':1, 'seq':None, 'transparency':OVERLAY_T, \
                'sim_time':-1}, meta)]
    else:
        # stage 1 position in animation
        frames = int(args.time * args.framerate)
        for i in xrange(frames):
            if s_azimuth <= 180:
                azimuth = 180 - (i / float(frames - 1)) * (180 - s_azimuth)
            else:
                azimuth = 180 + (i / float(frames - 1)) * (s_azimuth - 180)
            scale_t = min(float(i), meta['t_frames']) / meta['t_frames']
            tilt = 90 - (i / float(frames - 1)) * (90 - map_tilt)
            msg_list.append([timeslice, {'azimuth':azimuth, 'tilt':tilt, \
                    'scale_t':scale_t, 'seq':i, 'transparency':OVERLAY_T, \
                    'sim_time':-1}, meta])
        # slip fadeout - todo: copy last frame to begin with
        for i in xrange(meta['t_frames']):
            scale_t = 1 - i / (meta['t_frames'] - 1.0)
            over_t = 100 - (100 - OVERLAY_T) * scale_t
            msg_list.append([timeslice, {'azimuth':s_azimuth, 'tilt':map_tilt, \
                    'scale_t':scale_t, 'seq':frames + i, 'sim_time':-1, \
                    'transparency':over_t}, meta])
        frames += meta['t_frames']
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

        # dependency tracking
        if finished == None:
            pass

        elif finished[0] == load_xyts:
            meta['xyts_cpt_max'] = gmt.xyv_cpt_range('%s/xyts/pgv.bin' % (gmt_temp))[2]
            msg_deps -= 1

            # load xyts overlays
            for i in xrange(nproc):
                msg_list.append([load_xyts_ts, meta, {'start':i, 'inc':nproc}])
            msg_deps += nproc

        elif finished[0] == load_xyts_ts:
            ready = range(finished[2]['start'], xfile.nt, nproc)
            for i in xrange(frames_sr):
                # frames containing slip rate
                sim_time = float(i) / args.framerate
                xpos = int(round(sim_time / xfile.dt))
                if xpos in ready:
                    scale_t = min(float(i), meta['t_frames']) / meta['t_frames']
                    msg_list.append([timeslice, {'azimuth':s_azimuth, \
                            'tilt':map_tilt, 'scale_t':scale_t, 'scale_x':0.0, \
                            'seq':frames + i, 'transparency':OVERLAY_T, \
                            'sim_time':sim_time}, meta])
            for i in xrange(frames_sr, frames_gm):
                # frames containing only ground motion
                sim_time = float(i) / args.framerate
                xpos = int(round(sim_time / xfile.dt))
                scale_x = min(float(i - frames_sr), meta['t_frames']) \
                        / meta['t_frames']
                if xpos in ready:
                    msg_list.append([timeslice, {'azimuth':s_azimuth, \
                            'tilt':map_tilt, 'scale_t':1.0, 'scale_x':scale_x, \
                            'seq':frames + i, 'transparency':OVERLAY_T, \
                            'sim_time':sim_time}, meta])
            msg_deps -= 1

        if len(msg_list) == 0:
            if msg_deps == 0:
                # all jobs complete, kill off slaves
                msg_list.append(StopIteration)
                nproc -= 1
            else:
                # waiting for dependencies
                msg_list.append(None)

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
        # add dynamic frames to counter
        frames += max(frames_sr, frames_gm)
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

    # have to sleep if waiting for dependencies
    from time import sleep

    # ask for work until stop sentinel
    logbook = []
    for task in iter(lambda: comm.sendrecv(None, dest = MASTER), StopIteration):

        t0 = time()
        # no jobs available yet
        if task == None:
            sleep(1)
            logbook.append(('sleep', time() - t0))
        elif task[0] is load_xyts:
            load_xyts(task[1])
        elif task[0] is load_xyts_ts:
            load_xyts_ts(task[1], task[2])
        elif task[0] is timeslice:
            timeslice(task[1], task[2])
            logbook.append(('timeslice', task[1], time() - t0))
        else:
            print('Slave recieved unknown task to complete: %s.' % (task))

    # reports to master
    comm.gather(sendobj = logbook, root = MASTER)
    # shutdown
    comm.Disconnect()

