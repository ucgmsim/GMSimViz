#!/usr/bin/env python2

from glob import glob
import math
import os
from shutil import copy, move, rmtree
import sys
from tempfile import mkdtemp
from time import time

from mpi4py import MPI
import numpy as np
try:
    from h5py import File as h5open
except ImportError:
    print('Missing h5py module.')
    print('Will not be able to plot liquefaction or landslide.')

import qcore.geo as geo
import qcore.gmt as gmt
import qcore.srf as srf
import qcore.xyts as xyts

MASTER = 0
TILT_MIN = 40
TILT_MAX = 20
TILT_DIP = 1
PAGE_WIDTH = 16
PAGE_HEIGHT = 9
SCALE_WIDTH = PAGE_WIDTH / 1.618
SCALE_SIZE = 0.25
SCALE_PAD = 0.1
OVERLAY_T = 40
# 100 for 1600x900
# 120 for 1920x1080, 16ix9i
# 240 for 3840x2160
DPI = 240
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
    if not os.path.isdir(tmp):
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

def load_hdf5(h5file, basename):
    """
    Specific to liquefaction and landslide HDF5 files.
    """
    if not os.path.isdir(os.path.dirname(basename)):
        try:
            os.makedirs(os.path.dirname(basename))
        except IOError:
            pass

    # reformat data
    with h5open(h5file, 'r') as h:
        ylen, xlen = h['model'].shape
        data = np.empty((xlen, ylen, 3))
        data[:, :, 0] = np.repeat(h['x'], ylen).reshape(xlen, ylen)
        data[:, :, 1] = np.tile(h['y'][...][::-1], xlen).reshape(xlen, ylen)
        data[:, :, 2] = h['model'][...].T

    # clear low values because grid contains potentially rotated data (gaps)
    values = data[:, :, 2]
    low = np.nanpercentile(values, 1)
    np.nan_to_num(values)
    values[values <= low] = np.nan

    # calculate metadata
    x0, y0 = data[0, 0, :2]
    x1, y1 = data[1, 1, :2]
    dx = '%.2fk' % (geo.ll_dist(x0, y0, x1, y0) * 0.6)
    dy = '%.2fk' % (geo.ll_dist(x0, y0, x0, y1) * 0.6)
    region = (np.min(data[:, :, 0]), np.max(data[:, :, 0]), \
            np.min(data[:, :, 1]), np.max(data[:, :, 1]))

    # store data
    data.astype(np.float32).tofile('%s.bin' % (basename))
    gmt.table2grd('%s.bin' % (basename), '%s.nc' % (basename), \
            region = region, dx = dx, dy = dy, \
            climit = np.nanpercentile(data[:, :, 2], 10))
    # create corners, mask
    corners = [data[0, 0, :2], data[xlen - 1, 0, :2], \
            data[xlen - 1, ylen - 1, :2], data[0, ylen - 1, :2], data[0, 0, :2]]
    with open('%s.cnrs' % (basename), 'w') as c:
        c.write('> corners of data\n')
        c.write('\n'.join([' '.join(map(str, cnr)) for cnr in corners]))
    gmt.grd_mask('%s.cnrs' % (basename), '%s.mask' % (basename), \
            region = region, dx = dx, dy = dy)

    # cpt - liquefaction up to 0.6, landslide up to 0.4
    # assuming fixed file names
    if os.path.basename(basename) == 'liquefaction':
        cpt_max = 0.6
    elif os.path.basename(basename) == 'landslide':
        cpt_max = 0.4
    else:
        raise ValueError('Not implemented.')
    gmt.makecpt('hot', '%s.cpt' % (basename), 0, cpt_max, invert = True)

    # matching mask
    mask_path = '%s/xyts/corners-hr.gmt' % (meta['wd'])
    if os.path.exists(mask_path):
        gmt.grd_mask(mask_path, '%s_mask.nc' % (basename), \
                dx = dx, dy = dy, \
                region = region)

    return region

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
    gmt_ps = os.path.join(swd, '%s_perspective%s.ps' \
        % (os.path.splitext(os.path.basename(meta['srf_file']))[0], \
        '_%.4d' % (i) * (job['seq'] != None)))
    if os.path.exists('%s.png' % (os.path.splitext(gmt_ps)[0])):
        print('Sequence %d found.' % (i))
        return
    os.makedirs(swd)

    # allow fixed rotation (north azimuth)
    if meta['rot'] != 1000.0:
        job['azimuth'] = (meta['rot'] + 180) % 360

    # final view window
    lon1, lat1, dlon1, dlat1 = gmt.region_fit_oblique( \
            job['view'][0], job['azimuth'] - 90, wd = swd)
    dlon1 *= job['view'][1]
    dlat1 *= job['view'][1]
    # transitional view adjustment
    if len(job['view']) > 2:
        prog = job['view'][2]
        lon0, lat0, dlon0, dlat0 = gmt.region_fit_oblique( \
                job['view'][3], job['azimuth'] - 90, wd = swd)
        dlon0 *= job['view'][4]
        dlat0 *= job['view'][4]
        # adjust centre
        lon1 += (lon0 - lon1) * (1 - prog)
        lat1 += (lat0 - lat1) * (1 - prog)
        # adjust area
        dlon1 += (dlon0 - dlon1) * (1 - prog)
        dlat1 += (dlat0 - dlat1) * (1 - prog)

    km_region = (-dlon1, dlon1, -dlat1, dlat1)
    projection = 'OA%s/%s/%s/%s' % (lon1, lat1, job['azimuth'] - 90, PAGE_WIDTH)
    km_region = gmt.fill_space_oblique(lon1, lat1, PAGE_WIDTH, \
            PAGE_HEIGHT / math.sin(math.radians( \
            max(job['tilt'], meta['map_tilt']))), km_region, 'k', projection, \
            DPI / math.sin(math.radians(max(job['tilt'], meta['map_tilt']))), \
            swd)
    corners, llur = gmt.map_corners(projection = projection, \
            region = km_region, region_units = 'k', return_region = 'llur', \
            wd = swd)
    map_width, map_height = gmt.map_dimentions(projection = projection, \
            region = llur, wd = swd)
    map_height *= math.sin(math.radians(job['tilt']))
    # determine km per inch near centre of page for Z axis scaling
    mid_x = PAGE_WIDTH / 2.0
    mid_y1 = (PAGE_HEIGHT / 2.0 - 0.5) \
            / math.sin(math.radians(max(job['tilt'], meta['map_tilt'])))
    mid_y2 = (PAGE_HEIGHT / 2.0 + 0.5) \
            / math.sin(math.radians(max(job['tilt'], meta['map_tilt'])))
    mid_lls = gmt.mapproject_multi([[mid_x, mid_y1], [mid_x, mid_y2]], \
            wd = swd, projection = projection, region = llur, \
            inverse = True)
    km_inch = geo.ll_dist(mid_lls[0, 0], mid_lls[0, 1], \
            mid_lls[1, 0], mid_lls[1, 1])
    z_scale = -1.0 /  km_inch \
            * math.sin(math.radians(max(job['tilt'], meta['map_tilt'])))
    if job['tilt'] < meta['map_tilt']:
            # virtual tilt
            v_tilt = math.asin(job['tilt'] / meta['map_tilt'])
            a_tilt = math.sin(math.radians(meta['map_tilt']))
            z_scale = -1.0 / km_inch * (a_tilt + (1 - a_tilt) * math.cos(v_tilt))

    # begin plot
    p = gmt.GMTPlot(gmt_ps)
    # use custom page size
    gmt.gmt_defaults(wd = swd, \
            ps_media = 'Custom_%six%si' % (PAGE_WIDTH, PAGE_HEIGHT))
    # shortcut to switch between geographic and plot projections
    def proj(projected, shift = True):
        if projected:
            p.spacial('OA', llur, lon0 = lon1, lat0 = lat1, \
                    z = 'z%s' % (z_scale), \
                    sizing = '%s/%s' % (job['azimuth'] - 90, PAGE_WIDTH), \
                    p = '180/%s/0' % (job['tilt']), \
                    y_shift = (not shift) * PAGE_HEIGHT \
                        + shift * (PAGE_HEIGHT - map_height) / 2.0)
        else:
            p.spacial('X', \
                    (0, PAGE_WIDTH, (not shift) * -PAGE_HEIGHT, PAGE_HEIGHT), \
                    sizing = '%s/%s' % (PAGE_WIDTH, \
                        PAGE_HEIGHT + (not shift) * PAGE_HEIGHT), \
                    y_shift = (not shift) * -PAGE_HEIGHT \
                            + shift * (map_height - PAGE_HEIGHT) / 2.0)

    # draw sky and earth if map tilted beyond page filling point
    if job['tilt'] < meta['map_tilt']:
        p.spacial('X', (0, 1, -1, 1), \
                sizing = '%s/%s' % (PAGE_WIDTH, PAGE_HEIGHT))
        p.path('0 0\n1 0\n1 1\n0 1', is_file = False, close = True, \
                width = None, fill = 'p50+bskyblue+fwhite+r%s' % (DPI))
        p.path('0 0\n1 0\n1 -1\n0 -1', is_file = False, close = True, \
                width = None, fill = 'p30+bdarkbrown+fbrown+r%s' % (DPI))

    proj(True)
    p.basemap()

    # simulation domain
    if os.path.isfile('%s/xyts/corners.gmt' % (meta['wd'])):
        p.path('%s/xyts/corners.gmt' % (meta['wd']), close = True, \
                width = '2p', split = '-', colour = '60/60/60')
    p.sites(gmt.sites_major)

    # srf plane outline
    p.path(meta['gmt_bottom'], is_file = False, colour = 'black@30', width = '1p', \
            split = '-', close = True, z = True)
    # srf plane top edges
    p.path(meta['gmt_top'], is_file = False, colour = 'black', width = '2p', \
            z = True)

    # path plot such as road network
    if job['sim_time'] == -3:
        plot = 'paths'
        # there is no scale anyway
        scale_p_final = 1.0
        # draw path file
        if job['proportion'] >= 1:
            path_file = job['overlay']
        else:
            path_file = '%s/path_subselection.gmt' % (swd)
            gmt.proportionate_segs(job['overlay'], path_file, job['proportion'])
        # expected pen style to be defined in segment headers with -Wpen
        p.path(path_file)
    # plain surface plot for liquefaction and landslide data
    elif job['sim_time'] == -2:
        plot = 'surface'
        scale_p_final = 1.0
        crop_grd = '%s/overlay/%s_mask.nc' % (meta['wd'], job['overlay'])
        if not os.path.exists(crop_grd):
            crop_grd = None
        p.overlay('%s/overlay/%s.nc' % (meta['wd'], job['overlay']), \
                '%s/overlay/%s.cpt' % (meta['wd'], job['overlay']), \
                transparency = job['transparency'], \
                crop_grd = crop_grd, custom_region = job['region'])
    # load srf plane data
    elif job['sim_time'] == - 1:
        plot = 'slip'
        scale_p_final = 0.7
        if job['transparency'] < 100:
            srf_data = gmt.srf2map(meta['srf_file'], swd, prefix = 'plane', \
                    value = 'slip', cpt_percentile = 95, wd = swd, \
                    xy = True, dpu = DPI, \
                    pz = z_scale * math.cos(math.radians(job['tilt'])))
    elif job['sim_time'] >= 0:
        plot = 'timeseries'
        scale_p_final = 1.0
        regions_sr = []
        # TODO: major refactoring, already exists within gmt.srf2map
        srt = int(round(job['sim_time'] / meta['srf_dt']))
        if srt >= meta['sr_len']:
            srt = meta['sr_len'] - 1
        for i in xrange(meta['n_plane'] * (job['transparency'] < 100)):
            # lon, lat, depth
            subfaults = np.fromfile('%s/subfaults_%d.bin' % (meta['wd'], i), \
                    dtype = '3f')
            # reproject
            xyv = np.empty((subfaults.shape[0], 3))
            xyv[:, :2] = gmt.mapproject_multi(subfaults[:, :2], wd = swd, \
                    z = '-Jz%s' % (z_scale), p = True)
            xyv[:, 1] += subfaults[:, 2] \
                    * z_scale * math.cos(math.radians(job['tilt']))
            xyv[:, 2] = np.memmap('%s/sliptss_%d.bin' % (meta['wd'], i), \
                    dtype = 'f', shape = (len(subfaults), meta['sr_len']))[:, srt]
            # region
            x_min, y_min = np.min(xyv[:, :2], axis = 0)
            x_max, y_max = np.max(xyv[:, :2], axis = 0)
            regions_sr.append((x_min, x_max, y_min, y_max))
            # XY bounds
            bounds = []
            bounds_idx = [0, meta['planes'][i]['nstrike'] - 1, \
                    meta['planes'][i]['ndip'] * meta['planes'][i]['nstrike'] - 1, \
                    (meta['planes'][i]['ndip'] - 1) * meta['planes'][i]['nstrike']]
            for idx in bounds_idx:
                bounds.append(xyv[idx, :2])
            with open('%s/plane_%d_bounds.xy' % (swd, i), 'w') as bounds_f:
                for point in bounds:
                    bounds_f.write('%s %s\n' % tuple(point))
            # XY mask grid
            rc = gmt.grd_mask('%s/plane_%d_bounds.xy' % (swd, i), \
                    '%s/plane_%d_mask_xy.grd' % (swd, i), \
                    geo = False, dx = 1.0 / DPI, dy = 1.0 / DPI, \
                    region = regions_sr[i], wd = swd)
            if rc == gmt.STATUS_INVALID:
                # bounds are likely of area = 0, do not procede
                # caller should check if file below produced
                # attempted plotting could cause invalid postscript / crash
                continue
            # dump as binary
            xyv.astype(np.float32) \
                .tofile('%s/slip_%d.bin' % (swd, i))
            # search radius based on diagonal distance
            p2 = xyv[meta['planes'][i]['nstrike'] + 1, :2]
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
    proj(False, shift = False)
    if plot == 'slip' and job['transparency'] < 100:
        for s in xrange(len(srf_data[1])):
            if not os.path.exists('%s/plane_%d_slip_xy.grd' % (swd, s)):
                continue
            p.overlay('%s/plane_%d_slip_xy.grd' % (swd, s), \
                    '%s/slip.cpt' % (meta['srf_wd']), \
                    transparency = job['transparency'], \
                    crop_grd = '%s/plane_%d_mask_xy.grd' % (swd, s))
    elif plot == 'timeseries':
        for s in xrange(meta['n_plane']):
            if not os.path.exists('%s/slip_%d.grd' % (swd, s)):
                continue
            p.overlay('%s/slip_%d.grd' % (swd, s), \
                    '%s/slip.cpt' % (meta['srf_wd']), \
                    transparency = job['transparency'], \
                    crop_grd = '%s/plane_%d_mask_xy.grd' % (swd, s))
    # ground motion must be in geographic projection due to 3D Z scaling
    proj(True, shift = False)
    # srf hypocentre above slip, below ground motion
    p.points('%s %s %s\n' % (meta['hlon'], meta['hlat'], meta['hdepth']), \
            is_file = False, shape = 'a', size = 0.35, line = 'black', \
            line_thickness = '1p', z = True, clip = False)
    try:
        xpos = int(round(job['sim_time'] / meta['xyts_dt']))
    except KeyError:
        # xyts file has not been given
        xpos = -1
    gm_file = os.path.join(meta['wd'], 'xyts', 'ts%04d.nc' % (xpos))
    if os.path.isfile(gm_file):
        p.overlay3d(gm_file, cpt = '%s/xyts/gm.cpt' % (meta['wd']), \
                transparency = job['transparency'], dpi = DPI, \
                z = '-Jz%s' \
                % ((meta['gm_z_km'] * -z_scale) / meta['xyts_cpt_max']), \
                mesh = True, mesh_pen = '0.1p')
    p.rose('C', 'M', '1.8i', pos = 'rel', dxp = PAGE_WIDTH / 2.0 - 1.8, \
            dyp = map_height / 2.0 - 2.2 * (map_height / PAGE_HEIGHT), \
            fill = 'white@80', clearance = '0.2i', pen = 'thick,red')

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
    elif plot == 'surface':
        p.cpt_scale(PAGE_WIDTH / 2.0, scale_p, \
                '%s/overlay/%s.cpt' % (meta['wd'], job['overlay']), \
                length = SCALE_WIDTH, align = 'CT', dy = SCALE_PAD, \
                thickness = SCALE_SIZE, major = 0.1, minor = 0.02, \
                cross_tick = 0.02, label = job['cpt_label'])
    # cpt label
    try:
        assert(cpt_label != '')
        p.text(scale_margin, cpt_y, \
                cpt_label, align = 'RM', dx = - SCALE_PAD, size = 16)
    except (AssertionError, NameError):
        pass

    # title
    p.text(PAGE_WIDTH / 2.0, PAGE_HEIGHT, \
            os.path.basename(meta['srf_file']), align = 'RM', size = 26, \
            dy = WINDOW_T / -2.0, dx = - 0.2)
    # sim time
    if job['sim_time'] >= 0 and job['transparency'] < 100:
        p.text(PAGE_WIDTH - WINDOW_R, PAGE_HEIGHT - WINDOW_T, \
                '%.3fs' % (job['sim_time']), size = '24p', \
                align = "BR", font = 'Courier', dx = -0.2, dy = 0.1, \
                colour = 'black@%s' % (job['transparency']))

    if plot == 'slip' and job['transparency'] < 100:
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

    ###
    ### projection of the inner area to draw map ticks, legend
    ###
    p.spacial('OA%s/%s/%s/' % (lon1, lat1, job['azimuth'] - 90), \
            (str(llur_i[0][0]), str(llur_i[0][1]), \
            str(llur_i[1][0]), '%sr' % (llur_i[1][1])), \
            sizing = PAGE_WIDTH - WINDOW_L - WINDOW_R, \
            p = '180/%s/0' % (job['tilt']), \
            x_shift = WINDOW_L, y_shift = window_b)

    # MAP_FRAME_TYPE also changes cpt scales so cannot be globally set
    gmt.gmt_set(['MAP_FRAME_TYPE', 'inside', 'MAP_TICK_LENGTH_PRIMARY', '0.1i', 'MAP_ANNOT_OBLIQUE', '1'], wd = swd)
    if job['tilt'] < meta['map_tilt']:
        # don't show tick marks but still draw border line
        p.ticks(major = 0, minor = 0)
    else:
        p.ticks(major = '1d', minor = '0.2d')

    # draw legend on top
    if plot == 'paths':
        legend_file = '%s.legend' % (os.path.splitext(job['overlay'])[0])
        if os.path.exists(legend_file):
            gmt.gmt_set(['FONT_ANNOT_PRIMARY', '12p'])
            p.legend(legend_file, 'L', 'T', '4i', pos = 'rel', \
                    dx = 0.4, dy = + 0.2, transparency = job['transparency'], \
                    frame_fill = 'white@%s' % (65 + 0.35 * job['transparency']))

    # finish, clean up
    p.finalise()
    p.png(dpi = DPI * meta['downscale'], clip = False, out_dir = meta['wd'], \
            downscale = meta['downscale'])
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
    parser.add_argument('-p', '--ptime', help = 'animation pause time (s)', \
            type = float, default = 5.0)
    parser.add_argument('-d', '--delay', help = 'animation start delay (s)', \
            type = float, default = 1.5)
    parser.add_argument('-e', '--end', help = 'animation end delay (s)', \
            type = float, default = 3.0)
    parser.add_argument('-r', '--rot', help = 'fixed rotation (north deg)', \
            type = float, default = 1000.0)
    parser.add_argument('--downscale', type = int, default = 1, \
            help = 'ghostscript downscale factor (prevents jitter)')
    parser.add_argument('-q', '--liquefaction', help = 'liquefaction hdf5 filepath')
    parser.add_argument('-s', '--landslide', help = 'landslide hdf5 filepath')
    parser.add_argument('--temp', help = 'continue from previous temp dir')
    parser.add_argument('--paths', help = 'standard road network input')
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
    # liquefaction / landslide also optional
    if args.liquefaction != None:
        if not os.path.exists(args.liquefaction):
            print('Could not find liquefaction file: %s' % (args.liquefaction))
            sys.exit(2)
    if args.landslide != None:
        if not os.path.exists(args.landslide):
            print('Could not find landslide file: %s' % (args.landslide))
            sys.exit(2)
    # paths optional
    if args.paths != None:
        if not os.path.isdir(args.paths):
            print('Could not find path directory: %s' % (args.paths))
            sys.exit(2)
        path_sort = lambda n : int(os.path.basename(n).split('_')[0])
        path_files = sorted(glob('%s/*.gmt' % (args.paths)), key = path_sort)
        legend_files = sorted(glob('%s/*.legend' % (args.paths)), key = path_sort)
        if len(path_files) != len(legend_files) and len(legend_files) != 0:
            print('Inconsintent number of path files and legend files.')
            sys.exit(2)
        for i in xrange(len(path_files) * len(legend_files) != 0):
            if '%s.legend' % (path_files[i][:-4]) not in legend_files:
                print('Filenames are not matching between ".gmt" and ".legend".')
                sys.exit(2)

    # list of work for slaves to do
    msg_list = []
    # list of work prepared earlier but needs dependencies before run
    msg_list_post_xyts = []
    # post processing instructions
    op_list = []
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
    map_tilt = min(map_tilt, TILT_MIN)
    # plane domains
    bounds = srf.get_bounds(srf_file, depth = True)
    poi_srf = []
    for plane in bounds:
        for point in plane:
            poi_srf.append(point)
    hlon, hlat, hdepth = srf.get_hypo(srf_file, depth = True)
    top_left = bounds[0][0]
    top_right = bounds[-1][1]
    top_mid = geo.ll_mid(top_left[0], top_left[1], top_right[0], top_right[1])
    gmt_bottom = '\n>\n'.join(['\n'.join([' '.join(map(str, b)) \
            for b in p]) for p in bounds])
    gmt_top = '\n>\n'.join(['\n'.join([' '.join(map(str, b)) \
            for b in p[:2]]) for p in bounds])

    # working directory
    if args.temp != None:
        gmt_temp = os.path.abspath(args.temp)
        print('Resuming from: %s' % (gmt_temp))
    else:
        gmt_temp = mkdtemp(prefix = '_GMT_WD_PERSPECTIVE_', \
                dir = os.path.dirname(srf_file))

    meta = {'wd':gmt_temp, 'srf_file':srf_file, 's_azimuth':s_azimuth, \
            'map_tilt':map_tilt, 'hlon':hlon, 'hlat':hlat, 'hdepth':hdepth, \
            'gmt_bottom':gmt_bottom, 'gmt_top':gmt_top, 'animate':args.animate, \
            't_frames':int(args.mtime * args.framerate), 'srf_bounds':bounds, \
            'rot':args.rot, 'downscale':args.downscale}

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
    if not os.path.isdir(meta['srf_wd']):
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
    final_azimuth = 180
    # xyts quick preparation
    if xyts_file != None:
        if not os.path.isdir(os.path.join(gmt_temp, 'xyts')):
            os.makedirs(os.path.join(gmt_temp, 'xyts'))
        xfile = xyts.XYTSFile(xyts_file, meta_only = True)
        xcnrs = xfile.corners(gmt_format = True)
        xregion = xfile.region(corners = xcnrs[0])
        # TODO: just use dx?
        xres = '%sk' % (xfile.hh * xfile.dxts * 3.0 / 5.0)
        with open('%s/xyts/corners.gmt' % (gmt_temp), 'w') as xpath:
            xpath.write(xcnrs[1])
        poi_gm = xcnrs[0]
        meta['xyts_region'] = xregion
        meta['xyts_res'] = xres
        meta['xyts_dt'] = xfile.dt
        if args.animate:
            msg_list.append([load_xyts, meta])
            msg_deps += 1
            frames_gm = int(xfile.dt * (xfile.nt - 0.6) * args.framerate)
        # ground motion 3D Z extent based on sim domain size
        # final size also depends on map tilt angle
        xlen1 = geo.ll_dist(xcnrs[0][0][0], xcnrs[0][0][1], \
                xcnrs[0][1][0], xcnrs[0][1][1])
        xlen2 = geo.ll_dist(xcnrs[0][1][0], xcnrs[0][1][1], \
                xcnrs[0][2][0], xcnrs[0][2][1])
        meta['gm_z_km'] = max(xlen1, xlen2) \
                * math.sin(math.radians(meta['map_tilt']))
        # final view angle based on longer region size
        if xlen1 > xlen2:
            bearing_xl = geo.ll_bearing(xcnrs[0][0][0], xcnrs[0][0][1], \
                    xcnrs[0][1][0], xcnrs[0][1][1], midpoint = True) + 90
        else:
            bearing_xl = geo.ll_bearing(xcnrs[0][1][0], xcnrs[0][1][1], \
                    xcnrs[0][2][0], xcnrs[0][2][1], midpoint = True) + 90
        if abs(geo.angle_diff(s_azimuth, bearing_xl)) > 90:
            final_azimuth = (bearing_xl + 180) % 360
        else:
            final_azimuth = bearing_xl % 360
    else:
        poi_gm = poi_srf
    # transform over longer direction
    diff_azimuth = geo.angle_diff(s_azimuth, final_azimuth)
    if diff_azimuth < 0:
        diff_azimuth += 360
    else:
        diff_azimuth -= 360
    frames_azimuth = int(max(frames_sr, frames_gm) * 0.85)
    diff_azimuth /= frames_azimuth

    # prepare other data
    if args.liquefaction != None:
        region_liquefaction = load_hdf5(args.liquefaction, \
                '%s/overlay/liquefaction' % (gmt_temp))
        poi_liquefaction = [[region_liquefaction[0], region_liquefaction[2]], \
                [region_liquefaction[1], region_liquefaction[2]], \
                [region_liquefaction[1], region_liquefaction[3]], \
                [region_liquefaction[0], region_liquefaction[3]]]
    if args.landslide != None:
        region_landslide = load_hdf5(args.landslide, \
                '%s/overlay/landslide' % (gmt_temp))
        poi_landslide = [[region_landslide[0], region_landslide[2]], \
                [region_landslide[1], region_landslide[2]], \
                [region_landslide[1], region_landslide[3]], \
                [region_landslide[0], region_landslide[3]]]

    # tasks
    if not args.animate:
        msg_list.append((timeslice, {'azimuth':s_azimuth, 'tilt':map_tilt, \
                'scale_t':1, 'seq':None, 'transparency':OVERLAY_T, \
                'sim_time':-1, 'view':(poi_srf, 1.618)}, meta))
    else:
        # how long camera should pause when making a pause
        pause_frames = int(round(args.ptime * args.framerate))

        # stage 1 face slip
        frames_slip = int(args.time * args.framerate)
        for i in xrange(frames_slip):
            if s_azimuth <= 180:
                azimuth = 180 - (i / float(frames_slip - 1)) * (180 - s_azimuth)
            else:
                azimuth = 180 + (i / float(frames_slip - 1)) * (s_azimuth - 180)
            scale_t = min(float(i), meta['t_frames']) / meta['t_frames']
            tilt = 90 - (i / float(frames_slip - 1)) * (90 - map_tilt)
            msg_list.append([timeslice, {'azimuth':azimuth, 'tilt':tilt, \
                    'scale_t':scale_t, 'seq':i, 'transparency':OVERLAY_T, \
                    'sim_time':-1, 'view':(poi_srf, 1.618, \
                    i / (frames_slip - 1.0), poi_gm, 1.2)}, meta])
        frames2now = frames_slip

        # stage 2 dip below surface
        frames_dip = pause_frames
        for i in xrange(frames_dip):
            tilt2 = 90 - (90 - TILT_DIP) * (i + 1.0) / frames_dip
            tiltv = math.sin(math.radians(tilt2)) * map_tilt
            msg_list.append([timeslice, {'azimuth':s_azimuth, 'tilt':tiltv, \
                    'scale_t':1, 'seq':frames2now + i, 'sim_time':-1, \
                    'transparency':OVERLAY_T, 'view':(poi_srf, 1.618)}, meta])
        # pause, then reverse frames
        op_list.append(['DUP', frames2now + i, pause_frames])
        op_list.append(['REV', frames2now, frames_dip, pause_frames])
        # frames will be duplicated in reverse and paused
        frames2now += frames_dip * 3

        # stage 3 slip fadeout
        frames_fade = meta['t_frames']
        for i in xrange(frames_fade):
            scale_t = 1 - i / (frames_fade - 1.0)
            over_t = 100 - (100 - OVERLAY_T) * scale_t
            msg_list.append([timeslice, {'azimuth':s_azimuth, 'tilt':map_tilt, \
                    'scale_t':scale_t, 'seq':frames2now + i, 'sim_time':-1, \
                    'transparency':over_t, 'view':(poi_srf, 1.618)}, meta])
        frames2now += frames_fade

        # for tasks added later (with xyts file), reference to time = 0
        frames2sim = frames2now
        # frame index to zoom progression (zoom out within 12.5 seconds)
        if frames_gm / args.framerate >= 12.5:
            zfac = 12.5 * args.framerate * 0.2
        else:
            zfac = frames_gm * 0.618 * 0.2
        def i2p(i):
            return 0.5 + 0.5 * math.tanh(i / zfac - 3)
        # stage 4 slip animation if no xyts
        for i in xrange(frames_sr * (frames_gm == 0)):
            sim_time = float(i) / args.framerate
            if i >= frames_azimuth:
                azimuth = final_azimuth
            else:
                azimuth = s_azimuth + i * diff_azimuth
            scale_t = min(float(i), meta['t_frames']) / meta['t_frames']
            msg_list.append([timeslice, {'azimuth':azimuth, \
                    'tilt':map_tilt, 'scale_t':scale_t, \
                    'seq':frames2now + i, 'transparency':OVERLAY_T, \
                    'sim_time':sim_time, 'view':(poi_srf, 1.618)}, meta])
        # must work if there will be later tasks or not
        frames2now += max(frames_gm, frames_sr)

        # stage 5 camera reset
        frames_return = meta['t_frames'] * 4
        diff_azimuth_return = geo.angle_diff(final_azimuth, 180) / frames_return
        if frames_gm == 0:
            msgs = msg_list
        else:
            msgs = msg_list_post_xyts
        for i in xrange(frames_return):
            scale_t = 1 - min(1, i / (meta['t_frames'] - 1.0))
            over_t = 100 - (100 - OVERLAY_T) * scale_t
            azimuth = final_azimuth + i * diff_azimuth_return
            tilt = 90 - (1 - (i / float(frames_return - 1))) * (90 - map_tilt)
            msgs.append([timeslice, {'azimuth':azimuth, 'tilt':tilt, \
                    'scale_t':scale_t, 'seq':frames2now + i, \
                    'sim_time':max(frames_sr, frames_gm) / args.framerate, \
                    'transparency':over_t, 'scale_x':1.0, \
                    'view':(poi_gm, 1.2)}, meta])
        frames2now += frames_return

        # stage 6 fade into liquefaction and landslide
        frames2liquefaction = frames2now
        frames_liquefaction = meta['t_frames'] * (args.liquefaction != None)
        frames2landslide = frames2liquefaction \
                + (frames_liquefaction * 2 + pause_frames) \
                * (args.liquefaction != None)
        frames_landslide = meta['t_frames'] * (args.landslide != None)
        for i in xrange(meta['t_frames'] \
                * (args.liquefaction != None or args.landslide != None)):
            scale_t = i / (meta['t_frames'] - 1.0)
            over_t = 100 - (100 - OVERLAY_T) * scale_t
            if args.liquefaction != None:
                msgs.append([timeslice, {'azimuth':azimuth, 'tilt':tilt, \
                        'scale_t':scale_t, 'seq':frames2liquefaction + i, \
                        'sim_time':-2, 'region':region_liquefaction, \
                        'transparency':over_t, 'overlay':'liquefaction', \
                        'cpt_label':'Liquefaction Hazard Probability', \
                        'view':(poi_liquefaction, 1.2, scale_t, poi_gm, 1.2)}, \
                        meta])
            if args.landslide != None:
                msgs.append([timeslice, {'azimuth':azimuth, 'tilt':tilt, \
                        'scale_t':scale_t, 'seq':frames2landslide + i, \
                        'sim_time':-2, 'region':region_landslide, \
                        'transparency':over_t, 'overlay':'landslide', \
                        'cpt_label':'Landslide Hazard Probability', \
                        'view':(poi_landslide, 1.2, scale_t, poi_gm, 1.2)}, \
                        meta])
        # pause, reverse frames for animation
        if args.liquefaction != None:
            op_list.append(['DUP', frames2liquefaction + i, pause_frames])
            op_list.append(['REV', frames2liquefaction, frames_liquefaction, pause_frames])
        if args.landslide != None:
            op_list.append(['DUP', frames2landslide + i, pause_frames])
            op_list.append(['REV', frames2landslide, frames_landslide, pause_frames])
        frames2now = frames2landslide \
                + (frames_landslide * 2 + pause_frames) \
                * (args.landslide != None)

        # stage 7 paths such as road network
        if args.paths != None:
            # first status
            poi_paths0 = gmt.simplify_segs(path_files[0])
            poi_paths = []
            # kaikoura roads, bottom isn't interesting
            for poi in poi_paths0:
                if poi[1] > -44.5:
                    poi_paths.append([poi[0], poi[1]])
            for i in xrange(meta['t_frames']):
                scale_t = i / (meta['t_frames'] - 1.0)
                over_t = 100 - 100 * scale_t
                msg_list.append([timeslice, {'azimuth':azimuth, 'tilt':tilt, \
                        'scale_t':0, 'seq':frames2now + i, 'sim_time':-3, \
                        'transparency':over_t, 'overlay':path_files[0], \
                        'proportion':scale_t, \
                        'view':(poi_paths, 1.5, scale_t, poi_gm, 1.2)}, meta])
            op_list.append(['DUP', frames2now + i, pause_frames])
            frames2now = frames2now + meta['t_frames'] + pause_frames
            # mid statuses
            for i in xrange(1, len(path_files)):
                msg_list.append([timeslice, {'azimuth':azimuth, 'tilt':tilt, \
                        'scale_t':0, 'seq':frames2now, \
                        'sim_time':-3, 'transparency':0, \
                        'overlay':path_files[i], 'proportion':1, \
                        'view':(poi_paths, 1.5)}, meta])
                op_list.append(['DUP', frames2now, pause_frames - 1])
                frames2now += pause_frames
            # end fadeout
            for i in xrange(meta['t_frames']):
                scale_t = 1 - min(1, i / (meta['t_frames'] - 1.0))
                over_t = 100 - 100 * scale_t
                msg_list.append([timeslice, {'azimuth':azimuth, 'tilt':tilt, \
                        'scale_t':0, 'seq':frames2now + i, 'sim_time':-3, \
                        'transparency':over_t, 'overlay':path_files[-1], \
                        'proportion':scale_t, \
                        'view':(poi_paths, 1.5, scale_t, poi_gm, 1.2)}, meta])
            frames2now += meta['t_frames']

        # frames of pause at beginning / end of movie
        frames_end = int(args.end * args.framerate)
        frames_start = int(args.delay * args.framerate)
        # overall movie frame alterations
        op_list.append(['DUP', frames2now - 1, frames_end])
        op_list.append(['SHIFT', 0, frames2now + frames_end, frames_start])
        op_list.append(['DUP', frames_start, - frames_start])

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
                    if i >= frames_azimuth:
                        azimuth = final_azimuth
                    else:
                        azimuth = s_azimuth + i * diff_azimuth
                    scale_t = min(float(i), meta['t_frames']) / meta['t_frames']
                    msg_list.append([timeslice, {'azimuth':azimuth, \
                            'tilt':map_tilt, 'scale_t':scale_t, 'scale_x':0.0, \
                            'seq':frames2sim + i, 'transparency':OVERLAY_T, \
                            'sim_time':sim_time, \
                            'view':(poi_gm, 1.2, i2p(i), poi_srf, 1.618)}, \
                            meta])
            for i in xrange(frames_sr, frames_gm):
                # frames containing only ground motion
                sim_time = float(i) / args.framerate
                xpos = int(round(sim_time / xfile.dt))
                scale_x = min(float(i - frames_sr), meta['t_frames']) \
                        / meta['t_frames']
                if xpos in ready:
                    if i >= frames_azimuth:
                        azimuth = final_azimuth
                    else:
                        azimuth = s_azimuth + i * diff_azimuth
                    msg_list.append([timeslice, {'azimuth':azimuth, \
                            'tilt':map_tilt, 'scale_t':1.0, 'scale_x':scale_x, \
                            'seq':frames2sim + i, 'transparency':OVERLAY_T, \
                            'sim_time':sim_time, \
                            'view':(poi_gm, 1.2, i2p(i), poi_srf, 1.618)}, \
                            meta])
            msg_deps -= 1
            if msg_deps == 0:
                # return frames require xyts colour palette and last timeslice
                # enough other jobs, don't need to start asap
                msg_list.extend(msg_list_post_xyts)

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

    # output files prefix
    basename = os.path.splitext(os.path.basename(meta['srf_file']))[0]
    # frame operations shortcut
    def frame_op(op, seq_from, seq_to):
        op('%s/%s_perspective_%.4d.png' % (gmt_temp, basename, seq_from), \
                '%s/%s_perspective_%.4d.png' % (gmt_temp, basename, seq_to))
    # frame operation instructions
    for op in op_list:
        if op[0] == 'DUP':
            # 'DUP', dup_frame, dup_count (negative to dup backwards)
            for i in xrange(1, op[2] + 1) if (op[2] >= 0) else \
                    xrange(op[2], 0, 1):
                frame_op(copy, op[1], op[1] + i)
        elif op[0] == 'REV':
            # 'REV', seq_from, seq_len, rev_gap
            for i in xrange(op[2]):
                frame_op(copy, op[1] + i, op[1] + op[2] * 2 + op[3] - 1 - i)
        elif op[0] == 'SHIFT':
            # 'SHIFT', seq_from, seq_len, shift
            for i in xrange(op[1] + op[2] - 1, op[1] - 1, -1) \
                    if (op[3] > 0) else xrange(op[1], op[1] + op[2]):
                frame_op(move, i, i + op[3])

    if args.animate:
        # output movie
        gmt.make_movie('%s/%s_perspective_%%04d.png' % (gmt_temp, basename), \
                basename, fps = args.framerate, codec = 'libx264')
    else:
        # take result out of temp
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

