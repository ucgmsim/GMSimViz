#!/usr/bin/env python2
"""
Works on a standard local folder structure.
"""

from glob import glob
import os
import shutil
import sys

import numpy as np

import qcore_path
from shared import get_corners
import gmt
import timeseries as ts
from geo import ll_dist, path_from_corners

###
### verify parameters
###
try:
    # must expand eg: '.' becomes actual path
    base_dir = os.path.abspath(sys.argv[1])
    assert os.path.isdir(base_dir)
except IndexError:
    print('First parameter is the event folder path.')
    exit(1)
except AssertionError:
    print('Folder not found: %s.' % (base_dir))
    exit(1)
# rstrip incase folder ends with separator eg: '/folder/path/'
event_name = os.path.basename(base_dir.rstrip(os.sep))
event_name = '_'.join(event_name.split('_')[:3])
obs_dir = os.path.join(base_dir, 'GM', 'Obs')
obs_velbb = os.path.join(obs_dir, 'Data', 'velBB')
if not os.path.exists(obs_velbb):
    print('vellBB folder does not exist: %s.' % (obs_velbb))
    exit(1)
event_stats = os.path.join(base_dir, 'Stat', event_name, \
        '%s.ll' % (event_name))
if not os.path.exists(event_stats):
    print('Event stats not found: %s.' % (event_stats))
    exit(1)
###
### continue environment setup and importing
###
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, base_dir)
params = '%s/params_plot.py' % (base_dir)
if not os.path.exists(params):
    shutil.copyfile('%s/params_plot.template.py' % (script_dir), params)
import params_plot as plot
seisplot = plot.SEISMO

###
### check for sim data
###
try:
    import params_base as sim_params
    SIM_DIR = True
except ImportError:
    SIM_DIR = False

# find files in standard folder structure
try:
    sfs_modelparams = glob('%s/VM/Model/%s/*/model_params_*' % (base_dir, event_name))[0]
except IndexError:
    sfs_modelparams = None
try:
    sfs_srf = glob('%s/Src/Model/*/Srf/*.srf' % (base_dir))[0]
except IndexError:
    sfs_srf = None

out_dir = '%s/GM/Obs/Figures' % (base_dir)
if not os.path.exists(seisplot.wd):
    os.makedirs(seisplot.wd)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
# start GMT environment
b = gmt.GMTPlot('%s/%s.ps' % (seisplot.wd, seisplot.name))

# subset of event_stats to plot
plot_stats = os.path.join(seisplot.wd, 'plot_stats.ll')
files = glob('%s/*.000' % (obs_velbb))
stats = []
for f in files:
    stat = os.path.splitext(os.path.basename(f))[0]
    stats.append(stat)

###
### boundaries
###
if seisplot.region == None and (SIM_DIR or sfs_modelparams != None):
    try:
        corners, cnr_str = get_corners(sim_params.MODELPARAMS, gmt_format = True)
    except (NameError, IOError):
        corners, cnr_str = get_corners(sfs_modelparams, gmt_format = True)
    # path following sim domain curved on mercator like projections
    fine_path = path_from_corners(corners = corners, output = None)
    # fit simulation region
    x_min = min([xy[0] for xy in fine_path])
    x_max = max([xy[0] for xy in fine_path])
    y_min = min([xy[1] for xy in fine_path])
    y_max = max([xy[1] for xy in fine_path])
elif seisplot.region == None:
    # fit all values
    xy = np.loadtxt(event_stats, usecols = (0, 1), dtype = 'f')
    x_min, y_min = np.min(xy, axis = 0) - 0.1
    x_max, y_max = np.max(xy, axis = 0) + 0.1
else:
    x_min, x_max, y_min, y_max = statplot.region
ll_region = (x_min, x_max, y_min, y_max)
ll_avg = sum(ll_region[:2]) / 2, sum(ll_region[2:]) / 2
# distance and size is used to calculate ts_xlen, ts_ymax, min_dist
lon_dist = ll_dist(x_min, ll_avg[1], x_max, ll_avg[1])
lon_size, lat_size = gmt.mapproject(x_max, y_max, wd = seisplot.wd, \
        projection = 'M%s' % (seisplot.width), region = ll_region)
if seisplot.min_dist == None:
    seisplot.min_dist = lon_dist / 12.
if seisplot.ts_xlen == None:
    seisplot.ts_xlen = 0.25 * lon_size
if seisplot.ts_ymax == None:
    seisplot.ts_ymax = 0.25 * seisplot.ts_xlen

###
### STEP 1: REDUCE STATION DENSITY
###
# stations to be plotted (not too crowded)
spaced = []

# add specific wanted stations first
with open(event_stats, 'r') as esf:
    for line in esf:
        if line.split()[-1] in seisplot.wanted_stats:
            x, y, n = line.split()
            spaced.append([float(x), float(y), n])

# add other stations sufficiently spaced out
with open(event_stats, 'r') as esf:
    for line in esf:
        x, y, n = line.split()
        x = float(x)
        y = float(y)
        # check if even in plotting region
        if not x_min < x < x_max or not y_min < y < y_max:
            continue
        # if the data file for this station exists
        if n in stats:
            good = True
            for s in spaced:
                if ll_dist(x, y, s[0], s[1]) < seisplot.min_dist:
                    good = False
                    break
            if good:
                spaced.append([float(x), float(y), n])

# result of wanted / spaced stations
with open(plot_stats, 'w') as out:
    for s in spaced:
        out.write('%s %s %s\n' % (s[0], s[1], s[2]))

###
### STEP 2: PRODUCE SEISMOGRAMS
###
# generate seismo file
lons = []
lats = []
stats = []
obs_vts = []
max_v = 0
with open(plot_stats, 'r') as sf:
    for line in sf:
        lon, lat, stat = line.split()
        lons.append(float(lon))
        lats.append(float(lat))
        stats.append(stat)
        obs_vts.append(ts.read_ascii('%s/%s.090' % (obs_velbb, stat), t0 = True))
        max_v = max(np.max(np.abs(obs_vts[-1])), max_v)
yfac = float(seisplot.ts_ymax) / max_v
xfac = float(seisplot.ts_xlen) / len(obs_vts[0])
# must start fresh
if os.path.exists(seisplot.obs_src):
    os.remove(seisplot.obs_src)
for s, stat in enumerate(stats):
    x0, y0 = gmt.mapproject(lons[s], lats[s], wd = seisplot.wd, \
            projection = 'M%s' % (seisplot.width), region = ll_region)
    gmt.make_seismo(seisplot.obs_src, obs_vts[s], x0, y0, xfac, yfac)

###
### STEP 3: PLOT SEISMOGRAMS
###
# background can be larger as whitespace is later cropped
b.background(11, 11)
b.spacial('M', ll_region, sizing = seisplot.width, \
        x_shift = 1, y_shift = 1)
# title, fault model and velocity model subtitles
b.text(ll_avg[0], ll_region[3], seisplot.title, size = 20, dy = 0.2)
# topo, water, overlay cpt scale
b.basemap()
# stations
b.points(event_stats, shape = 't', size = 0.08, \
        fill = None, line = 'white', line_thickness = 0.8)
b.points(plot_stats, shape = 't', size = 0.08, \
        fill = 'blue', line = 'white', line_thickness = 0.8)
with open(plot_stats, 'r') as psf:
    for line in psf:
        lon, lat, name = line.split()
        b.text(lon, lat, name, dy = - 0.05, align = 'CT', size = '8p', colour = 'blue')
if SIM_DIR or sfs_srf != None:
    # fault file - creating direct from SRF is slower
    # OK if only done in template
    if SIM_DIR and os.path.exists(sim_params.srf_files[0]):
        b.fault(sim_params.srf_files[0], is_srf = True, \
                plane_width = 0.5, top_width = 1, hyp_width = 0.5, \
                plane_colour = 'red', top_colour = 'red', hyp_colour = 'red')
    elif sfs_srf != None:
        b.fault(sfs_srf, is_srf = True, \
                plane_width = 0.5, top_width = 1, hyp_width = 0.5, \
                plane_colour = 'red', top_colour = 'red', hyp_colour = 'red')
    else:
        print('SRF file not found, not adding fault planes to plot.')
b.coastlines()

if SIM_DIR or sfs_modelparams != None:
    # simulation domain
    b.path(cnr_str, is_file = False, split = '-', \
            close = True, width = '0.4p', colour = 'black')
# ticks on top otherwise parts of map border may be drawn over
b.ticks(major = '60m', minor = '30m', sides = 'ws')

# seismograms are in xy format
b.spacial('X', (0, lon_size, 0, lat_size), sizing = '%s/%s' % (lon_size, lat_size))
# add seismograms - should be 'inc' format for static images
if seisplot.obs_src != None:
    b.seismo(seisplot.obs_src, seisplot.max_ts, fmt = 'inc', \
            width = seisplot.seis_width, colour = seisplot.obs_colour)
#if seisplot.sim_src != None:
#    b.seismo(seisplot.sim_src, seisplot.max_ts, fmt = 'inc', \
#            width = seisplot.seis_width, colour = seisplot.sim_colour)

# render and do not continue
b.finalise()
b.png(dpi = seisplot.dpi, out_dir = out_dir, clip = True)

shutil.rmtree(seisplot.wd)
