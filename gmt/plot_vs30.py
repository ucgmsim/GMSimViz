#!/usr/bin/env python2
"""
Prepares vs30 data to plot (specific) and plots with plot_stations.py (generic).
"""

from glob import glob
import os
from shutil import rmtree
from subprocess import Popen
import sys

vs_script = '/home/kmf76/Vs30-mapping/GMsim_Vs30_extract/extractVs30_raster.R'
vs_data = '/home/kmf76/VsMap/Rdata/20170404/hyb_NZGD00_allNZ.Rdata'
plot_script = os.path.join(os.path.abspath(os.path.dirname(__file__)), \
        'plot_stations.py')
if not os.path.exists(vs_script):
    print('Kevin\'s vs30 script has moved.')
    exit(1)
if not os.path.exists(vs_data):
    print('Kevin\'s vs30 data file has moved.')
    exit(1)

# first / only parameter is base directory
try:
    base_dir = os.path.abspath(sys.argv[1])
    assert(os.path.exists(base_dir))
except IndexError:
    print('Base directory parameter required.')
    exit(1)
except AssertionError:
    print('Base directory not found: %s' % (base_dir))
base_name = os.path.basename(base_dir)
base_name = '_'.join(base_name.split('_')[:3])
figures_dir = os.path.join(base_dir, 'VM', 'Figures')
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

# station file has a standard location
stat_file = os.path.join(base_dir, 'Stat', base_name, '%s.ll' % (base_name))
if not os.path.exists(stat_file):
    print('Station file missing, expected at %s' % (stat_file))
    exit(1)

# input to Kevin's script is tab-delimited and requires a header
vs30_tmp = '%s.vs30_tmp' % (os.path.splitext(stat_file)[0])
with open(vs30_tmp, 'w') as inf:
    inf.write('Site.Longitude\tSite.Latitude\tSite.Code\n')
    with open(stat_file, 'r') as sf:
        for line in sf:
            inf.write('%s\n' % ('\t'.join(line.split())))

# run Kevin's code
vs30_out = os.path.join(figures_dir, 'Vs30StationMap.ll')
Popen(['Rscript', vs_script, vs_data, vs30_tmp, vs30_out], \
        cwd = os.path.dirname(vs_script)).wait()
os.remove(vs30_tmp)

# reformat for plot_stations.py
with open(vs30_tmp, 'w') as inf:
    # plot_stations.py header
    inf.write('Velocity Model\nVs30 (m/s)\nhot:invert 0.2\n0 1000 50 250\n1\n\n')
    with open(vs30_out, 'r') as vsf:
        # Kevin header
        vsf.readline()
        for line in vsf:
            stat, lon, lat, vs30 = line.split()
            if vs30 == 'NA':
                vs30 = '500'
            inf.write('%s %s %s\n' % (lon, lat, vs30))
os.remove(vs30_out)
os.rename(vs30_tmp, vs30_out)

# call plotting script
Popen([plot_script, vs30_out], cwd = base_dir).wait()
# find output - questionable logic (could be rel/abs path)
sys.path.insert(0, base_dir)
from params_plot import STATION
out_dir = os.path.join(base_dir, os.path.basename(STATION.out_dir))
output = glob('%s/*.png' % (out_dir))[0]
os.rename(output, os.path.join(figures_dir, 'Vs30StationMap.png'))
# original output location not needed
rmtree(out_dir)
