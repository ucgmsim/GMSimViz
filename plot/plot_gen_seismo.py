#!/usr/bin/env python2
"""
Creates seismograms for use in GMT (see shared_gmt.py::seismo).
Inputs are ASCII Vel files.
For ideal/fast calculations:
    x and y scales are based on degrees longitude and latitude
    output map should have parallel lon/lat lines such as 'M'ercator

TODO:
could have displacements in length units by projecting origin into length
these would then be added on top after changing to an 'X' projection
"""

from math import sin, cos
import os
import sys
from time import time

import numpy as np

from params import *
import params_plot as plot
from qcore.timeseries import *
from qcore.geo import *

seismoplot = plot.SEISMO

for output in [seismoplot.obs_src, seismoplot.sim_src]:
    if os.path.exists(output):
        os.remove(output)

###
### read values, calculate maximums
###
lons = []
lats = []
stats = []
obs_vts = []
sim_vts = []
max_v = 0
xy = open(seismoplot.plot_stats, 'r')
for line in xy:
    lon, lat, stat = line.split()
    lat = float(lat)
    lon = float(lon)
    data_obs = seismoplot.obs_ts % (stat)
    data_sim = seismoplot.sim_ts % (stat)

    # first try simulated as the station may be outside domain
    # TODO: optionally read lower resolution data from XYTS file
    try:
        svt = read_ascii(data_sim, t0 = True)
        sim_vts.append(svt)
    except:
        print('Could not open SIM data for station: %s' % (data_sim))
        continue

    # add corresponding observed seismogram
    ovt = read_ascii(data_obs, t0 = True)
    obs_vts.append(ovt)

    # store metadata
    lats.append(lat)
    lons.append(lon)
    stats.append(stat)

    # PGV
    max_v = max(np.max(np.abs(ovt)), np.max(np.abs(svt)), max_v)
    #max_v = max(np.max(np.abs(ovt)), max_v)
xy.close()

print('Max Velocity: %s' % (max_v))

# y scaling factor based on overall max (km)
yfac = float(yamp) / max_v
# x scaling based on simulated time
xfac = float(tlen) / len(sim_vts[0])

yfac = float(seismoplot.ts_ymax) / max_v
xfac = float(seismoplot.ts_xlen) / len(sim_vts[0])
###
### generate XY data based on read values and max / other params
### simple version where seismo line is extended
###
for s, stat in enumerate(stats):

    offset_obs, oazim_obs = 0, 45
    offset_sim, oazim_sim = 0, 135
    # offset of beginning of seismogram
    lat0_obs, lon0_obs = ll_shift(lats[s], lons[s], offset_obs, oazim_obs)
    lat0_sim, lon0_sim = ll_shift(lats[s], lons[s], offset_sim, oazim_sim)

    # rotate as wanted (do not visually flip seismo)
    if seismoplot.ts_xaz % 360 - 180 < 0:
        yazim = seismoplot.ts_xaz - 90
    else:
        yazim = seismoplot.ts_xaz + 90

    # show the extending line to source
    lls_obs = ['%f %f\n' % (lons[s], lats[s])]
    lls_sim = ['%f %f\n' % (lons[s], lats[s])]
    for i, value in enumerate(sim_vts[s]):
        dy = value * yfac
        dx = i * xfac
        lls_sim.append('%f %f\n' % (lon0_sim + dx, lat0_sim + dy))
    for i, value in enumerate(obs_vts[s]):
        dy = value * yfac
        dx = i * xfac
        lls_obs.append('%f %f\n' % (lon0_obs + dx, lat0_obs + dy))

    with open(seismoplot.sim_src, 'a') as sp:
        sp.write('> station at %f %f\n' % (lons[s], lats[s]))
        sp.write(''.join(lls_sim))
    with open(seismoplot.obs_src, 'a') as sp:
        sp.write('> station at %f %f\n' % (lons[s], lats[s]))
        sp.write(''.join(lls_obs))
