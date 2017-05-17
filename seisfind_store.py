#!/usr/bin/env python2
"""
Retrieves X, Y, Z components from SEIS files.
Processes ready for BB.
Saves to Binary, ready for BB.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 13 September 2016

ISSUES: All options of EMOD3D are not considered. eg: scale.
    https://github.com/h5py/h5py/issues/428
        - strings must be fixed length for endianness interoperability
"""

from glob import glob
from math import sin, cos, radians
from multiprocessing import Pool
from os import remove, path, makedirs
from struct import unpack
import sys

import numpy as np
import h5py as h5

from shared_bin import *
from shared_ts import *
from shared_xyts import XYTSFile
from siteamp_models import *

from params import *
from params_base import *
from params_base_bb import hf_sim_basedir

# process some parameters to be more useful
srf_base = os.path.splitext(os.path.basename(srf_files[0]))[0]
lf_output = os.path.join(lf_sim_root_dir, srf_base, 'OutBin')
hf_output = os.path.join(hf_sim_basedir, srf_base, 'Acc')
hf_base = os.path.join(hf_output, hf_prefix)
sf_output = os.path.join(bb_dir, 'seisfinder')
bb_output = os.path.join(sf_output, 'bin')
# xyts file used to retrieve simulation parameters
# less likely to change format or missmatch
xyts_file = os.path.join(lf_output, '%s_xyts.e3d' % (run_name))
# velocity model resources
vmod_rho = os.path.join(vel_mod_dir, 'rho3dfile.d')
vmod_vs = os.path.join(vel_mod_dir, 'vs3dfile.s')
vmod_vp = os.path.join(vel_mod_dir, 'vp3dfile.p')


procs = 64
if len(sys.argv) > 1:
    procs = int(sys.argv[1])

# check dependencies
seis_file_list = glob(path.join(lf_output, '*_seis-?????.e3d'))
assert(len(seis_file_list) > 0)
for dependency in [xyts_file, hf_output, vmod_rho, vmod_vs, vmod_vp]:
    try:
        assert(os.path.exists(dependency))
    except AssertionError:
        print('FATAL: could not find %s' % (dependency))
        exit(1)
if not os.path.exists(bb_output):
    os.makedirs(bb_output)

# python data-type format string
INT_S = 'i'
FLT_S = 'f'
# automatically swap bytes if required
if get_seis_swap(seis_file_list[0]):
    swapping_char = get_byteswap_char()
    INT_S = swapping_char + INT_S
    FLT_S = swapping_char + FLT_S

# read common data once
nt, dt, hh, mrot = get_seis_common(seis_file_list[0], INT_S, FLT_S)
# additional metadata retrieval
xytf = XYTSFile(xyts_file, meta_only = True)
rot_matrix = xytf.rot_matrix
mlon = xytf.mlon
mlat = xytf.mlat
nx = xytf.nx_sim
ny = xytf.ny_sim
dx_ts = xytf.dxts
dy_ts = xytf.dyts

# TODO: should have values for each location
vsite = 250.0

# read every station in seis files
# keep virtual ones (name of 8 c-string bytes)
def process_seis_file(index):
    filename = seis_file_list[index]
    print('processing %s...' % (filename))
    fp = open(seis_file_list[index], 'rb')
    # speed optimise
    seek = fp.seek
    read = fp.read
    # shortcuts for processing input
    read_int = lambda : unpack(INT_S, read(SIZE_INT))[0]
    read_flt = lambda : unpack(FLT_S, read(SIZE_FLT))[0]

    # number of stations within this seis file
    num_stat = read_int()

    # map of array doesn't use RAM unless accessed
    comp_data = np.memmap(seis_file_list[index], \
            dtype = FLT_S, mode = 'r', \
            offset = SIZE_INT + num_stat * SIZE_SEISHEAD, \
            shape = (nt, num_stat, N_COMPS))
    rho = np.memmap(vmod_rho, dtype = '<f4', \
            shape = (ny, int(nz), nx))
    vs = np.memmap(vmod_vs, dtype = '<f4', \
            shape = (ny, int(nz), nx))
    vp = np.memmap(vmod_vp, dtype = '<f4', \
            shape = (ny, int(nz), nx))

    v_stat_info = []
    for stat_i in xrange(num_stat):
        seek(SIZE_INT + (stat_i + 1) * SIZE_SEISHEAD - STAT_CHAR)
        stat = str(read(STAT_CHAR)).rstrip('\0')
        if len(stat) == VSTAT_LEN:
            # station is a virtual entry
            seek( - STAT_CHAR - 5 * SIZE_FLT - 4 * SIZE_INT, 1)
            x = read_int()
            y = read_int()
            # location in format (XXXXYYYY)
            g_name = '%s%s' % (str(x).zfill(4), str(y).zfill(4))

            # read HF pairs
            try:
                hf = np.array([ \
                        np.fromfile('%s_%s.%s' % \
                                (hf_base, stat, MY_COMPS[0]), '>f4'), \
                        np.fromfile('%s_%s.%s' % \
                                (hf_base, stat, MY_COMPS[1]), '>f4'), \
                        np.fromfile('%s_%s.%s' % \
                                (hf_base, stat, MY_COMPS[2]), '>f4')])
            except IOError:
                print('Warning: Cannot find HF files for %s.' % g_name)
                print('Skipping point %s.' % (g_name))
                continue

            pga = [0, 0, 0]
            ampf = [0, 0, 0]
            vref = [vs[y][0][x], vs[y][0][x], vp[y][0][x]]
            vpga = vref
            # process HF TODO: len(hf) -> hf.shape[0] ?
            for tsi in xrange(len(hf)):
                # get pga
                pga[tsi] = np.max(np.abs(hf[tsi])) / 981.0
                # amplification factors
                ampf[tsi] = cb_amp(dt, get_ft_len(nt), \
                        vref[tsi], vsite, vpga[tsi], pga[tsi], \
                        version = "2008")
                # apply amplification factors
                hf[tsi] = ampdeamp(hf[tsi], ampf[tsi], amp = True)
                # filter
                hf[tsi] = bwfilter(hf[tsi], dt, 1.0, 'highpass', \
                        match_powersb = match_powersb)

            # pass on station info
            seek(2 * SIZE_INT + 3 * SIZE_FLT, 1)
            v_stat_info.append({'NAME':stat, 'X':x, 'Y':y, \
                    'LAT':read_flt(), 'LON':read_flt(), 'VSITE':vsite, \
                    'PGA_0':pga[0], 'PGA_1':pga[1], 'PGA_2':pga[2], \
                    'RHO':rho[y][0][x], 'VS':vref[0], 'VP':vref[-1]})

            # read LF pairs, rotate and reorientate
            lf = np.dot(np.dstack(( \
                    comp_data[..., stat_i, 0], \
                    comp_data[..., stat_i, 1], \
                    comp_data[..., stat_i, 2])), rot_matrix)[0].T

            # difference in start timestep of LF, HF
            ddt = int(1.0 / dt)
            bb = np.empty(shape = (N_MY_COMPS, nt + ddt))
            # process LF, HF
            for tsi in xrange(len(lf)):
                # convert to acceleration
                lf[tsi] = vel2acc(lf[tsi], dt)
                # apply amplification factors
                lf[tsi] = ampdeamp(lf[tsi], ampf[tsi], amp = True)
                # filter
                lf[tsi] = bwfilter(lf[tsi], dt, 1.0, 'lowpass', \
                        match_powersb = match_powersb)

                # broadband
                bb[tsi] = np.cumsum((np.hstack(([0] * ddt, hf[tsi])) + \
                        np.hstack((lf[tsi], [0] * ddt))) * dt)
            # save broadband to file
            bb.T.astype(np.float32).tofile('%s/bb_%s' % \
                    (bb_output, g_name), format = '<f4')

    fp.close()
    return v_stat_info

###
### debug single threaded
###
#seis_results = []
#for i in xrange(len(seis_file_list)):
#    seis_results.append(process_seis_file(i))

###
### multiprocessing
###
p = Pool(procs)
seis_results = p.map(process_seis_file, range(len(seis_file_list)))

print('Broadband Files Complete')
print('Storing metadata for export in HDF5')

h5p = h5.File('%s/metadata.hdf5' % (sf_output), 'w')
# add common metadata
h5p.attrs['NSTAT'] = 0
h5p.attrs['MLON'] = mlon
h5p.attrs['MLAT'] = mlat
h5p.attrs['MROT'] = mrot
h5p.attrs['NX'] = nx
h5p.attrs['NY'] = ny
h5p.attrs['DX_TS'] = dx_ts
h5p.attrs['DY_TS'] = dy_ts
h5p.attrs['NT'] = nt
h5p.attrs['DT'] = dt
h5p.attrs['HH'] = hh
h5p.attrs['NCOMPS'] = N_MY_COMPS
for n, key in enumerate(MY_COMPS):
    h5p.attrs['COMP_%d' % (n)] = np.string_(MY_COMPS[key])

# add station data
for file_results in seis_results:
    h5p.attrs['NSTAT'] += len(file_results)
    for stat_i in file_results:
        g_name = '%s%s' % (str(stat_i['X']).zfill(4), str(stat_i['Y']).zfill(4))
        f_name = 'bb_%s' % (g_name)
        h5group = h5p.create_group(g_name)
        # have to save strings as fixed width because of h5py bug
        h5group.attrs['NAME'] = np.string_(stat_i['NAME'])
        for key in ['X', 'Y', 'LAT', 'LON', 'VSITE', \
                'PGA_0', 'PGA_1', 'PGA_2', 'RHO', 'VS', 'VP']:
            h5group.attrs[key] = stat_i[key]

h5p.close()
