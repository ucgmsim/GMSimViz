#!/usr/bin/env python2

"""
Created: 4 January 2017
Purpose: Generate visualisations of obs/sim ratios, PGA, PGV, PSA.
Authors: Viktor Polak <viktor.polak@canterbury.ac.nz>

USAGE:
Execute with python: "$ ./plot_stations.py" or "$ python2 plot_stations.py"
First parameter is the file to plot.
Second parameter is optional override for output folder.

INPUT FORMAT:
File to plot must be in the following format:
Note numbers are lines of the file.
1. Plot Title (blank line for no title)
2. Legend Title for the colour palette scale
3. cpt source, station size
cpt source and station size can be followed by comma separated properties after a ':'
..::cpt source examples::..
hot:invert,t-40 will invert the hot palette and display with 40% transparency
hot:fg-black,bg-white set foreground and background colour
..::station size examples::..
0.2:shape-c will make the station size 0.2 with a circle shape
1k:g-nearneighbor will make grid spacing 1km and use the nearneighbor algorithm
4. min cpt, max cpt, cpt inc, cpt legend tick
all optionoal but cpt min and max must both be provided or none at all
optional parameters must be in order
5. number of data colums excluding longitude and latitude, optional label colour
6. comma separated column labels. placed inside top left corner of map.
Optional but must either be of length 0 or number of columns

7 - END. Data are longitude, latitude, col_1, col_2 ... col_N

ISSUES:
"""

from glob import glob
import os
from shutil import copy, rmtree
import sys
sys.path.append('.')
from tempfile import mkdtemp
from time import time

from mpi4py import MPI
import numpy as np

import qcore_path
import gmt
from shared import get_corners
import geo
from srf import srf2corners

MASTER = 0

# check for sim data
try:
    import params_base as sim_params
    SIM_DIR = True
except ImportError:
    SIM_DIR = False

# find files in standard folder structure
# working directory must be base of structure for this functionality
base_dir = os.path.abspath(os.getcwd())
event_name = os.path.basename(base_dir.rstrip(os.sep))
event_name = '_'.join(event_name.split('_')[:3])
try:
    sfs_modelparams = glob('%s/VM/Model/%s/*/model_params_*' % (base_dir, event_name))[0]
except IndexError:
    sfs_modelparams = None
try:
    sfs_srf = glob('%s/Src/Model/*/Srf/*.srf' % (base_dir))[0]
except IndexError:
    sfs_srf = None

# temporary working directories for gmt are within here
# prevent multiprocessing issues by isolating processes
gmt_temp = os.path.join(os.path.abspath('.'), 'GMT_WD_STATIONS')
ps_template = '%s/template.ps' % (gmt_temp)
txt_cnrs = '%s/srf_cnrs.txt' % (gmt_temp)
model_path = '%s/sim.modelpath_hr' % (gmt_temp)
mask = '%s/mask.grd' % (gmt_temp)

# process file header
def load_file(station_file):
    # load values
    val_pool = np.atleast_2d( \
            np.loadtxt(station_file, dtype = 'f', skiprows = 6)[:, 2:].T)
    ncol = val_pool.shape[0]

    with open(station_file) as statf:
        head = [next(statf).strip() for _ in xrange(6)]
        # 1st - title
        title = head[0]

        # 2nd line - legend title
        legend = head[1]

        # 3rd line - cpt description 1
        # src, point size, foreground colour, background colour
        cpt_info = head[2].split()
        cpt = cpt_info[0].split(':')[0]
        # default properties
        transparency = 0
        cpt_fg = None
        cpt_bg = None
        if os.path.exists(cpt):
            # assuming it is a built in cpt if not matching filename
            cpt = os.path.abspath(cpt)
        try:
            # src:invert will add the 'invert' property to invert cpt
            cpt_properties = cpt_info[0].split(':')[1].split(',')
            for p in cpt_properties:
                if p[:2] == 't-':
                    transparency = p[2:]
                elif p[:3] == 'fg-':
                    cpt_fg = p[3:]
                elif p[:3] == 'bg-':
                    cpt_bg = p[3:]
        except IndexError:
            cpt_properties = []
        if len(cpt_info) > 1:
            stat_size = cpt_info[1].split(':')[0]
            # also default nearneighbor search distance
            nn_search = stat_size
            # default properties
            shape = 't'
            grid = None
            grd_mask_dist = None
            try:
                stat_properties = cpt_info[1].split(':')[1].split(',')
                for p in stat_properties:
                    if p[:6] == 'shape-':
                        shape = p[6]
                    elif p[:2] == 'g-':
                        grid = p[2:]
                    elif p[:3] == 'nns-':
                        nn_search = p[3:]
                    elif p[:6] == 'gmask-':
                        grd_mask_dist = p[6:]
            except IndexError:
                stat_properties = []

        # 4th line - cpt description 2
        # cpt_min, cpt_max, cpt_inc, cpt_tick
        cpt_info2 = head[3].split()
        if len(cpt_info2) > 1:
            usr_min, usr_max = map(float, cpt_info2[:2])
            cpt_min = [usr_min] * ncol
            cpt_max = [usr_max] * ncol
        else:
            cpt_min = []
            cpt_max = np.percentile(val_pool, 99.5, axis = 1)
            cpt_inc = []
            cpt_tick = []
            for i in xrange(len(cpt_max)):
                if cpt_max[i] > 115:
                    # 2 significant figures
                    cpt_max[i] = round(cpt_max[i], 1 - int(np.floor(np.log10(cpt_max[i]))))
                else:
                    # 1 significant figures
                    cpt_max[i] = round(cpt_max[i], - int(np.floor(np.log10(cpt_max[i]))))
                if val_pool[i].min() < 0:
                    cpt_min.append(-cpt_max)
                else:
                    cpt_min.append(0)
                cpt_inc.append(cpt_max[i] / 10.)
                cpt_tick.append(cpt_inc[i] * 2.)
        if len(cpt_info2) > 2:
            cpt_inc = [float(cpt_info2[2])] * ncol
        if len(cpt_info2) > 3:
            cpt_tick = [float(cpt_info2[3])] * ncol

        # 5th line ncols and optional column label prefix
        col_info = head[4].split()
        ncol = int(col_info[0])
        if len(col_info) > 1:
            label_colour = col_info[1]
        else:
            label_colour = 'black'

        # 6th line - column labels
        col_labels = map(str.strip, head[5].split(','))
        if col_labels == ['']:
            col_labels = []
        if len(col_labels) != ncol and len(col_labels) != 0:
            print('%d column labels found when there are %d columns.' \
                    % (len(col_labels), ncol))
            exit(1)

    return {'title':title, 'legend':legend, 'stat_size':stat_size, \
            'nn_search':nn_search, 'shape':shape, 'grid':grid, \
            'grd_mask_dist':grd_mask_dist, 'cpt':cpt, 'cpt_fg':cpt_fg, \
            'cpt_bg':cpt_bg, 'cpt_min':cpt_min, 'cpt_max':cpt_max, \
            'cpt_inc':cpt_inc, 'cpt_tick':cpt_tick, \
            'transparency':transparency, 'ncol':ncol, \
            'label_colour':label_colour, 'col_labels':col_labels}

###
### boundaries
###
def determine_sizing(station_file, meta):
    if statplot.region == None and (SIM_DIR or sfs_modelparams != None):
        try:
            corners, cnr_str = get_corners(sim_params.MODELPARAMS, gmt_format = True)
        except (NameError, IOError):
            corners, cnr_str = get_corners(sfs_modelparams, gmt_format = True)
        # path following sim domain curved on mercator like projections
        fine_path = geo.path_from_corners(corners = corners, output = None)
        # fit simulation region
        x_min = min([xy[0] for xy in fine_path])
        x_max = max([xy[0] for xy in fine_path])
        y_min = min([xy[1] for xy in fine_path])
        y_max = max([xy[1] for xy in fine_path])
    elif statplot.region == None:
        # fit all values
        xy = np.loadtxt(station_file, skiprows = 6, usecols = (0, 1), dtype = 'f')
        x_min, y_min = np.min(xy, axis = 0) - 0.1
        x_max, y_max = np.max(xy, axis = 0) + 0.1
    else:
        x_min, x_max, y_min, y_max = statplot.region
    # combined region
    ll_region = (x_min, x_max, y_min, y_max)
    # avg lon/lat (midpoint of plotting region)
    ll_avg = sum(ll_region[:2]) / 2, sum(ll_region[2:]) / 2

    # create masking if using grid overlay
    if meta['grid'] != None:
        try:
            path_from_corners(corners = corners, min_edge_points = 100, \
                    output = model_path)
        except NameError:
            path_from_corners(corners = [ \
                    [x_min, y_min], [x_max, y_min], \
                    [x_max, y_max], [x_min, y_max]], \
                    min_edge_points = 100, \
                    output = model_path)
        grd_mask(model_path, mask, dx = meta['stat_size'], \
                dy = meta['stat_size'], region = ll_region)

    # work out an ideal tick increment (ticks per inch)
    # x axis is more constrainig
    if statplot.tick_major == None:
        try:
            width = float(statplot.width)
        except ValueError:
            # expecting a unit suffix even though formula only works for inches
            width = float(statplot.width[:-1])
        statplot.tick_major, statplot.tick_minor = \
                gmt.auto_tick(x_min, x_max, width)
    elif statplot.tick_minor == None:
        statplot.tick_minor = statplot.tick_major / 5.
    # cities/locations to plot
    if statplot.sites == None:
        region_sites = []
    if statplot.sites == 'auto':
        if x_max - x_min > 3:
            region_sites = gmt.sites_major
        else:
            region_sites = gmt.sites.keys()
    elif statplot.sites == 'major':
        region_sites = gmt.sites_major
    elif statplot.sites == 'all':
        region_sites = gmt.sites.keys()
    else:
        region_sites = statplot.sites

    return {'ll_region':ll_region, 'll_avg':ll_avg, \
            'major_tick':statplot.tick_major, \
            'minor_tick':statplot.tick_minor, 'sites':region_sites}

def template(plot):
    """
    Initialises incomplete template base file.
    Upon completion other threads can start on copies before basemap completes.
    """
    # incomplete template for common working GMT conf/history files
    t = gmt.GMTPlot(ps_template)
    # background can be larger as whitespace is later cropped
    t.background(11, 15)
    t.spacial('M', plot['ll_region'], sizing = statplot.width, \
            x_shift = 1, y_shift = 2.5)
    t.leave()

def template_2(cnr_str = None):
    """
    Finish teplate file which will be placed below rest of images.
    Top layer of template is slowest but can run simultaneous to others.
    """
    t0 = time()
    # work in separate folder
    # prevents GMT IO errors on its conf/history files
    wd = '%s/basemap' % (gmt_temp)
    if not os.path.isdir(wd):
        os.makedirs(wd)
    # name of complete template postscript
    ps = '%s%s%s' % (swd, os.sep, os.path.basename(ps_template))
    # copy GMT setup and basefile
    copy('%s/gmt.conf' % (gmt_temp), wd)
    copy('%s/gmt.history' % (gmt_temp), wd)
    copy(ps_template, ps)
    t = GMTPlot(ps, append = True)

    # topo, water, overlay cpt scale (slow)
    t.basemap()
    # simulation domain - optional
    if cnr_str != None:
        t.path(cnr_str, is_file = False, split = '-', \
                close = True, width = '0.4p', colour = 'black')
    t.leave()
    print('Basemap prepared (%.2f).' % (time() - t0))

def fault_prep(srf_file = None, cnrs_file = None):
    """
    Makes sure simple SRF bounds description is available.
    Prevents re-loading the SRF file.
    """
    # fault path - boundaries already available
    if cnrs_file != None:
        copyfile(cnrs_file, txt_cnrs)
        print('Corners file found.')
    # fault path - determine boundaries (slow)
    elif srf_file != None:
        t0 = time()
        srf2corners(srf_file, cnrs = txt_cnrs)
        print('Corners prepared (%.2f).' % (time() - t0))
    else:
        print('No fault will be plotted.')

def column_overlay(n, meta, plot):
    """
    Produces map for a column of data.
    """
    t0 = time()

    # prepare resources in separate folder
    # prevents GMT IO errors on its conf/history files
    swd = '%s/c%.3dwd' % (gmt_temp, n)
    if not os.path.isdir(swd):
        os.makedirs(swd)
    # name of slice postscript
    ps = '%s/c%.3d.ps' % (swd, n)

    # copy GMT setup and append top layer to blank file
    copy('%s/gmt.conf' % (gmt_temp), swd)
    copy('%s/gmt.history' % (gmt_temp), swd)
    p = gmt.GMTPlot(ps, append = True)

    # prepare cpt
    cpt_stations = '%s/stations.cpt' % (swd)
    # overlay colour scale
    gmt.makecpt(meta['cpt'], cpt_stations, meta['cpt_min'][n], meta['cpt_max'][n], \
            inc = meta['cpt_inc'][n], invert = 'invert' in meta['cpt_properties'], \
            fg = meta['cpt_fg'], bg = meta['cpt_bg'], \
            transparency = meta['transparency'])

    # common title
    if len(meta['title']):
        p.text(plot['ll_avg'][0], plot['ll_region'][3], meta['title'], \
                colour = 'black', align = 'CB', size = 28, dy = 0.2)

    # add ratios to map
    if meta['grid'] == None:
        p.points(station_file, shape = meta['shape'], \
                size = meta['stat_size'], fill = None, line = None, \
                cpt = cpt_stations, cols = '0,1,%d' % (n + 2), header = 6)
    else:
        grd_file = '%s/overlay.grd' % (swd)
        if meta['grd_mask_dist'] != None:
            col_mask = '%s/column_mask.grd' % (swd)
            mask = col_mask
        else:
            col_mask = None
        gmt.table2grd(station_file, grd_file, file_input = True, \
                grd_type = meta['grid'], region = plot['ll_region'], \
                dx = meta['stat_size'], climit = meta['cpt_inc'][n] * 0.5, \
                wd = swd, geo = True, sectors = 4, min_sectors = 1, \
                search = meta['nn_search'], cols = '0,1,%d' % (n + 2), \
                header = 6, automask = col_mask, \
                mask_dist = meta['grd_mask_dist'])
        p.overlay(grd_file, cpt_stations, dx = meta['stat_size'], \
                dy = meta['stat_size'], crop_grd = mask, land_crop = False, \
                transparency = meta['transparency'])
    # add locations to map
    p.sites(plot['sites'])

    # title for this data column
    if len(meta['col_labels']):
        p.text(plot['ll_region'][0], plot['ll_region'][3], \
                meta['col_labels'][n], colour = meta['label_colour'], \
                align = 'LB', size = '18p', dx = 0.2, dy = -0.35)

    # ticks on top otherwise parts of map border may be drawn over
    p.ticks(major = plot['major_tick'], minor = plot['minor_tick'], sides = 'ws')

    # colour scale
    p.cpt_scale(3, -0.5, cpt_stations, meta['cpt_tick'][n], \
        meta['cpt_inc'][n], label = meta['legend'], \
        arrow_f = meta['cpt_max'][n] > 0, arrow_b = meta['cpt_min'][n] < 0)

    print('Column %d complete in %.2fs' % (n, time() - t0))

def render(n):
    """
    Rendering postscript can be slow (by amount of details)
    """
    ps = '%s/c%.3dwd/c%.3d.ps' % (gmt_temp, n, n)
    p = gmt.GMTPlot(ps, append = True)

    # reduces dependency on main creation function
    # add srf to map
    if os.path.exists(txt_cnrs):
        p.fault(txt_cnrs, is_srf = False, \
                plane_width = 0.5, top_width = 1, hyp_width = 0.5)

    p.finalise()
    p.png(dpi = statplot.dpi, out_dir = statplot.out_dir, clip = True)

###
### MASTER
###
if len(sys.argv) > 1:
    # timing
    t_start = MPI.Wtime()

    sysproc = int(os.sysconf('SC_NPROCESSORS_ONLN'))
    # TODO: allow override
    nproc_max = sysproc

    ###
    ### PREPARE
    ###
    script_dir = os.path.abspath(os.path.dirname(__file__))
    if not os.path.exists('params_plot.py'):
        copyfile('%s/params_plot.template.py' % (script_dir), 'params_plot.py')
    import params_plot
    statplot = params_plot.STATION
    try:
        station_file = os.path.abspath(sys.argv[1])
        assert(os.path.exists(station_file))
    except IndexError:
        print('First parameter must be input file. Parameter not found.')
        exit(1)
    except AssertionError:
        print('Cannot find input file: %s' % (station_file))
        exit(1)
    # allow overriding output directories more easily when scripting
    if len(sys.argv) > 2:
        statplot.out_dir = os.path.abspath(sys.argv[2])
    # clear output dirs
    for out_dir in [statplot.out_dir, gmt_temp]:
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
    # load station file
    meta = load_file(station_file)
    meta['out_dir'] = statplot.out_dir
    # calculate other parameters
    plot = determine_sizing(station_file, meta)
    # template postscript
    template(plot)

    # tasks
    task_list = [(template_2), (column_overlay)]
    nproc = min(len(task_list), nproc_max)
    # stop signals for each slave
    msg_list = task_list + ([StopIteration] * nproc)
    # spawn slaves
    comm = MPI.COMM_WORLD.Spawn(
        sys.executable, args = [sys.argv[0]], maxprocs = nproc)
    # distribute work to slaves who ask
    status = MPI.Status()
    for position, msg in enumerate(msg_list):
        comm.recv(source = MPI.ANY_SOURCE, status = status)
        comm.send(obj = msg, dest = status.Get_source())

    # gather reports from slaves
    reports = comm.gather(None, root = MPI.ROOT)

    # statistics
    print('\nFinished in %.2fs.' % (MPI.Wtime() - t_start))
    # clear all working files
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

    ###
    ### PREPARE
    ###
    from gmt import *
    import params_plot
    statplot = params_plot.STATION

    # ask for work until stop sentinel
    logbook = []
    for task in iter(lambda: comm.sendrecv(None, dest = MASTER), StopIteration):
        logbook.append(task)

        print(task)

    # reports to master
    comm.gather(sendobj = logbook, root = MASTER)
    # shutdown
    comm.Disconnect()
