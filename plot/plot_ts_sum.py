#!/usr/bin/env python2
"""
Plot timeslice summaries.
Generates LON,LAT,VEL IEEE754 32bit (single precision) binary file.
Reduces TimeSlice files for PGV and MMI.
Creates MMI input for USGS PAGER

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 9 June 2016

USAGE: 

ISSUES: currently assumes ABSMAX TS input.
"""

from datetime import datetime
from glob import glob
from math import ceil, sqrt
import os
from shutil import copy, rmtree
import sys

import h5py as h5
import numpy as np

import qcore_path
from geo import path_from_corners, wgs_nztm2000x
import gmt
from xyts import XYTSFile
# to import local parameters add them to path
sys.path.append(os.path.abspath(os.path.curdir))
from params_base import *
srf_cnr = srf_cnrs[0]
xyts_file = xyts_files[0]
import params_plot as plot
pgvplot = plot.PGV
mmiplot = plot.MMI

# working directory for intermediate files
gmt_temp = os.path.abspath('GMT_WD_TSSUM')
# output directory for PNG outputs
png_dir = os.path.abspath('PNG_tssum')

# clear output dirs
for out_dir in [gmt_temp, png_dir]:
    if os.path.isdir(out_dir):
        rmtree(out_dir)
    os.makedirs(out_dir)

###
### PART 1: generate inputs to plots
###
#########################################################
print('Reducing data for PGV and MMI...')
scenarios = len(xyts_files)
if scenarios > 1:
    print('%d scenarios detected.' % (scenarios))
    pgvs = []
    mmis = []
    for s in xrange(scenarios):
        print('PGV and MMI for scenario %d...' % (s + 1))
        xyts = XYTSFile(xyts_files[s])
        pgv_file = '%s/PGV_%.3d.bin' % (gmt_temp, s)
        mmi_file = '%s/MMI_%.3d.bin' % (gmt_temp, s)
        xyts.pgv(mmi = True, pgvout = pgv_file, mmiout = mmi_file)
        pgvs.append(np.fromfile(pgv_file, dtype = '3f')[:, 2])
        mmis.append(np.fromfile(mmi_file, dtype = '3f')[:, 2])
    print('Producing minimum and maximum for PGV and MMI...')
    # retrieve longitude and latitude for values
    template = np.fromfile(pgv_file, dtype = '3f')
    prefix = 'PGV'
    for values in pgvs, mmis:
        mn = np.min(np.array(values), axis = 0)
        template[:, 2] = mn
        template.astype(np.float32).tofile('%s/%s_MIN.bin' % (gmt_temp, prefix))
        mx = np.max(np.array(values), axis = 0)
        template[:, 2] = mx
        template.astype(np.float32).tofile('%s/%s_MAX.bin' % (gmt_temp, prefix))
        prefix = 'MMI'
else:
    xyts = XYTSFile(xyts_files[0])
    xyts.pgv(mmi = True, pgvout = '%s/PGV.bin' % (gmt_temp), \
            mmiout = '%s/MMI.bin' % (gmt_temp))
print('PGV and MMI ready.')

# second step: reformat for lifeline analysis
# re-save PGV with NZTM2000 coordinates as CSV
print('Creating NZTM2000 PGV grid for lifeline analysis...')
# grab a template
if scenarios == 1:
    xyv = np.fromfile('%s/PGV.bin' % (gmt_temp), dtype = '3f')
else:
    xyv = np.fromfile('%s/PGV_000.bin' % (gmt_temp), dtype = '3f')
xyv[:, :2] = wgs_nztm2000x(xyv[:, :2])
if scenarios == 1:
    suffixes = [None]
    np.savetxt('%s/nztm2000pgv.txt' % (png_dir), xyv, \
            fmt = '%f', delimiter = ',')
else:
    suffixes = ['%.3d' % (i) for i in xrange(scenarios)]
    suffixes.extend(['MIN', 'MAX'])
    for suffix in suffixes:
        xyv[:, 2] = np.fromfile('%s/PGV_%s.bin' \
                % (gmt_temp, suffix), dtype = '3f')[:, 2]
        np.savetxt('%s/nztm2000pgv_%s.txt' % (png_dir, suffix), xyv, \
                fmt = '%f', delimiter = ',')
print('NZTM2000 PGV grid complete.')

# third step: reformat for PAGER input
# PAGER requires a square grid, aranged left -> right, top -> bottom
# uses GMT to resample data - most resources are used for plotting also
print('Reformatting MMI for USGS PAGER...')
# area info
corners, cnr_str = xyts.corners(gmt_format = True)
corners = np.array(corners)
x_min, y_min = np.min(corners, axis = 0)
x_max, y_max = np.max(corners, axis = 0)
plot_region = (x_min, x_max, y_min, y_max)
# grid - spacing
if pgvplot.grd_dx == None or pgvplot.grd_dy == None:
    dx = '%sk' % (xyts.dx / 2.0)
    dy = dx
else:
    dx = pgvplot.grd_dx
    dy = pgvplot.grd_dy
# masking
mask = '%s/modelmask.grd' % (gmt_temp)
path_from_corners(corners = corners.tolist(), min_edge_points = 100, \
        output = '%s/sim.modelpath_hr' % (gmt_temp))
gmt.grd_mask('%s/sim.modelpath_hr' % (gmt_temp), mask, \
        dx = dx, dy = dy, region = plot_region)

### above are shared for plotting section
### below is unique for PAGER
event_id = run_name
mag = run_name.split('_')[1].split('-')[0].lstrip('m').replace('pt', '.')
dep = 8.66
event_type = 'SCENARIO' # or 'SCENARIO'
origin_time = '2017-04-25T13:02:33.631Z'

for s in xrange(len(suffixes)):
    if scenarios == 1:
        grid_out = '%s/grid.xml' % (png_dir)
        mmi_in = '%s/MMI.bin' % (gmt_temp)
    else:
        grid_out = '%s/grid_%s.xml' % (png_dir, suffixes[s])
        mmi_in = '%s/MMI_%s.bin' % (gmt_temp, suffixes[s])

    # use GMT to resample data on square grid
    # create cropped grid with minimum MMI = 1 (roman numeral minimum)
    gmt.table2grd(mmi_in, '%s/PAGER_MMI.grd' % (gmt_temp), \
            region = plot_region, dx = dx, dy = dy, climit = 0.1)
    # AND returns B if A == NaN, else A
    gmt.grdmath(['%s/PAGER_MMI.grd' % (gmt_temp), mask, 'MUL', \
            1, 'AND', 1, 'MAX', '=', '%s/PAGER_MMI.grd' % (gmt_temp)], \
            dx = dx, dy = dy, region = plot_region)

    # retrieve grid data
    grd = h5.File('%s/PAGER_MMI.grd' % (gmt_temp))
    # lat stored min -> max but pager requires top -> bottom
    lons = grd['lon'][...]
    lats = grd['lat'][...]
    mmis = grd['z'][...]
    grd_ny, grd_nx = mmis.shape
    grd.close()

    try:
        with open(srf_cnrs[s], 'r') as scnrs:
            hyp = scnrs.readline()
            while hyp[0] == '>':
                hyp = scnrs.readline()
        hlon, hlat = hyp.split()
    except IndexError:
        # MIN and MAX don't have specified hypocentres
        hlon, hlat = '', ''

    pager_grid = open(grid_out, 'w')
    pager_grid.write('<?xml version="1.0" \
encoding="UTF-8" standalone="yes"?>\n')
    pager_grid.write('<shakemap_grid \
xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" \
xmlns="http://earthquake.usgs.gov/eqcenter/shakemap" \
xsi:schemaLocation="http://earthquake.usgs.gov \
http://earthquake.usgs.gov/eqcenter/shakemap/xml/schemas/shakemap.xsd" \
event_id="%s" shakemap_id="%s" \
shakemap_version="1" code_version="1" \
process_timestamp="%s" \
shakemap_originator="nz" \
map_status="RELEASED" \
shakemap_event_type="%s">\n' % (event_id, event_id, \
            datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"), event_type))
    pager_grid.write('<event event_id="%s" magnitude="%s" \
depth="%s" lat="%s" lon="%s" event_timestamp="%s" \
event_network="nz" event_description="%s" />\n' % (event_id, mag, \
            dep, hlat, hlon, origin_time, run_name.split('_')[0]))
    pager_grid.write('<grid_specification lon_min="%s" lat_min="%s" \
lon_max="%s" lat_max="%s" nominal_lon_spacing="%s" \
nominal_lat_spacing="%s" nlon="%d" nlat="%d" />\n' % (x_min, y_min, \
            x_max, y_max, abs(x_max - x_min) / grd_nx, \
            abs(y_max - y_min) / grd_ny, grd_nx, grd_ny))
    pager_grid.write('<event_specific_uncertainty name="mi" \
value="-1" numsta="0" />\n')
    pager_grid.write('<grid_field index="1" name="LON" units="dd" />\n')
    pager_grid.write('<grid_field index="2" name="LAT" units="dd" />\n')
    pager_grid.write('<grid_field index="3" name="MMI" units="intensity" />\n')
    pager_grid.write('<grid_data>\n')
    for a in xrange(len(lats) - 1, - 1, - 1):
        for b in xrange(len(lons)):
            pager_grid.write('%s %s %s\n' % (lons[b], lats[a], mmis[a, b]))
    pager_grid.write('</grid_data>\n</shakemap_grid>\n')
    pager_grid.close()

print('PAGER MMI input ready.')

###
### PART 2: plotting
###
#########################################################
print ('Plotting PGV and MMI...')
# resources
ll_avg = (x_min + x_max) / 2.0, (y_min + y_max) / 2.0
if (x_max - x_min) > 3:
    plot_sites = gmt.sites_major
else:
    plot_sites = gmt.sites.keys()
if pgvplot.major_tick == None:
    try:
        width = float(pgvplot.width)
    except ValueError:
        # expecting a unit suffix even though formula only works for inches
        width = float(pgvplot.width[:-1])
    pgvplot.tick_major, pgvplot.minor_tick = \
            gmt.auto_tick(x_min, x_max, width)
elif pgvplot.minor_tick == None:
    pgvplot.minor_tick = pgvplot.major_tick / 5.
# colour palettes
cpt_pgv = '%s/pgv.cpt' % (gmt_temp)
cpt_mmi = '%s/cpt/mmi.cpt' % (os.path.abspath(os.path.dirname(__file__)))
if scenarios == 1:
    reference = '%s/PGV.bin' % (gmt_temp)
else:
    reference = '%s/PGV_MAX.bin' % (gmt_temp)
if pgvplot.cpt_inc == None or pgvplot.cpt_max == None:
    pgvplot.cpt_inc, pgvplot.cpt_max = \
            gmt.xyv_cpt_range(reference, \
            my_inc = pgvplot.cpt_inc, my_max = pgvplot.cpt_max)[1:3]
if pgvplot.convergence_limit == None:
    pgvplot.convergence_limit = pgvplot.cpt_inc * 0.2
plot.vel_model = plot.vel_model.replace('<HH>', str(xyts.hh))
gmt.makecpt(pgvplot.cpt, cpt_pgv, pgvplot.cpt_min, pgvplot.cpt_max, \
        inc = pgvplot.cpt_inc, invert = pgvplot.cpt_inv, \
        bg = pgvplot.cpt_bg, fg = pgvplot.cpt_fg)

# create common plot first
p = gmt.GMTPlot('%s/PGV.ps' % (gmt_temp))
p.background(11, 11)
p.spacial('M', plot_region, sizing = pgvplot.width, \
        x_shift = 1, y_shift = 2.5)
# common titles
p.text(x_min, y_max, plot.fault_model, size = 14, align = 'LB', dy = 0.3)
p.text(x_min, y_max, plot.vel_model, size = 14, align = 'LB', dy = 0.1)
# common features
p.basemap()
# stations - split into real and virtual
with open(stat_file, 'r') as sf:
    stations = sf.readlines()
stations_real = []
stations_virtual = []
for i in xrange(len(stations)):
    if len(stations[i].split()[-1]) == 7:
        stations_virtual.append(stations[i])
    else:
        stations_real.append(stations[i])
p.points(''.join(stations_real), is_file = False, \
        shape = 't', size = 0.08, fill = None, \
        line = 'white', line_thickness = 0.8)
p.points(''.join(stations_virtual), is_file = False, \
        shape = 'c', size = 0.02, fill = 'black', line = None)

# diverge into plot sub-types
p.leave()
copy('%s/PGV.ps' % (gmt_temp), '%s/MMI.ps' % (gmt_temp))

for s in xrange(len(suffixes)):
    if scenarios == 1:
        pgv_ps = '%s/PGV.ps' % (gmt_temp)
        mmi_ps = '%s/MMI.ps' % (gmt_temp)
        pgv_in = '%s/PGV.bin' % (gmt_temp)
        mmi_in = '%s/MMI.bin' % (gmt_temp)
    else:
        pgv_ps = '%s/PGV_%s.ps' % (gmt_temp, suffixes[s])
        mmi_ps = '%s/MMI_%s.ps' % (gmt_temp, suffixes[s])
        pgv_in = '%s/PGV_%s.bin' % (gmt_temp, suffixes[s])
        mmi_in = '%s/MMI_%s.bin' % (gmt_temp, suffixes[s])
        copy('%s/PGV.ps' % (gmt_temp), pgv_ps)
        copy('%s/MMI.ps' % (gmt_temp), mmi_ps)

    p = gmt.GMTPlot(pgv_ps, append = True)
    m = gmt.GMTPlot(mmi_ps, append = True)

    p.overlay(pgv_in, cpt_pgv, dx = dx, \
            dy = dy, climit = pgvplot.convergence_limit, \
            min_v = pgvplot.lowcut, crop_grd = mask, \
            contours = pgvplot.cpt_inc, land_crop = pgvplot.land_crop)
    m.overlay(mmi_in, cpt_mmi, dx = dx, \
            dy = dy, climit = mmiplot.convergence_limit, min_v = 0, \
            crop_grd = mask, contours = 1, land_crop = mmiplot.land_crop)
    p.sites(plot_sites)
    m.sites(plot_sites)
    if scenarios == 1:
        p.fault(srf_cnrs[0], is_srf = False, plane_width = 0.5, \
                top_width = 1, hyp_width = 0.5)
        m.fault(srf_cnrs[0], is_srf = False, plane_width = 0.5, \
                top_width = 1, hyp_width = 0.5)
    p.coastlines()
    m.coastlines()
    p.text(ll_avg[0], y_max, pgvplot.title, size = 20, dy = 0.6, align = 'CB')
    m.text(ll_avg[0], y_max, mmiplot.title, size = 20, dy = 0.6, align = 'CB')
    p.cpt_scale(3, -0.5, cpt_pgv, pgvplot.cpt_inc, pgvplot.cpt_inc, \
            label = pgvplot.cpt_legend, arrow_f = pgvplot.cpt_max > 0, \
            arrow_b = pgvplot.cpt_min < 0)
    m.cpt_scale(3, -0.5, cpt_mmi, 1, 1, \
            label = mmiplot.cpt_legend)
    p.ticks(major = pgvplot.major_tick, minor = pgvplot.minor_tick)
    m.ticks(major = pgvplot.major_tick, minor = pgvplot.minor_tick)
    p.path(cnr_str, is_file = False, split = '-', \
            close = True, width = '0.4p', colour = 'black')
    m.path(cnr_str, is_file = False, split = '-', \
            close = True, width = '0.4p', colour = 'black')

    # save plots
    p.finalise()
    m.finalise()
    p.png(out_dir = png_dir, dpi = pgvplot.dpi)
    m.png(out_dir = png_dir, dpi = mmiplot.dpi)

# clear temporary files
rmtree(gmt_temp)
