#!/usr/bin/env python2
"""
Function A (prepare):
Prepares and Packages data on Fitzroy

Function B (retrieve):
Retrieves data from Fitzroy ready for plotting on Hypocentre

Function C (auto):
Runs function A remotely, Runs function B locally, cleans up.
"""

import os
from shutil import copyfile
from subprocess import Popen
import sys
import tarfile as tf
from tempfile import mkstemp

# srf corners used repeatedly - precalculate
from qcore.srf import srf2corners

def prepare(sd, tp):
    """
    Prepare plotting dependencies in `sd` sim_dir for transport.
    sd: sim_dir base dir for simulation to prepare
    tp: tar filepath, where to store output
    """
    # tar file will have 1 folder containing resources
    # ie. do not create a tarbomb
    prefix = os.path.basename(os.path.abspath(sd))
    # import parameters from sim_dir (resource locations)
    sys.path.insert(0, sd)
    from params_base import run_name, srf_files, FD_STATLIST

    # create tar file
    tar = tf.open(tp, 'w', format = tf.PAX_FORMAT)
    # params_base resources with relative paths
    p_base = mkstemp()[1]
    copyfile('%s/params_base.py' % (sd), p_base)
    pbf = open(p_base, 'a')
    pbf.write('# FOLLOWING APPENDED BY TRANSFER SCRIPT\n')

    # plotting resources
    tar.add(FD_STATLIST, arcname = '%s/%s' % (prefix, os.path.basename(FD_STATLIST)))
    pbf.write('stat_file = \'%s\'\n' % (os.path.basename(FD_STATLIST)))

    # each SRF file will also produce different TSFiles
    srf2cnr = lambda path : '%s/cnrs.txt' \
            % (os.path.basename(os.path.splitext(path)[0]))
    srf2xyts = lambda path : '%s/xyts.e3d' \
            % (os.path.basename(os.path.splitext(path)[0]))
    for srf in srf_files:
        srf_name = os.path.basename(srf)
        # only corners needed to plot SRF
        corners = mkstemp()[1]
        srf2corners(srf, corners)
        tar.add(corners, arcname = '%s/%s/cnrs.txt' % (prefix, srf_name[:-4]))
        os.remove(corners)
        tar.add('%s/LF/%s/OutBin/%s_xyts.e3d' % (sd, srf_name[:-4], run_name), \
                arcname = '%s/%s/xyts.e3d' % (prefix, srf_name[:-4]))
    pbf.write('srf_cnrs = %s\n' % (str(map(srf2cnr, srf_files))))
    pbf.write('xyts_files = %s\n' % (str(map(srf2xyts, srf_files))))

    # can also set this as absolute path upon retrieval
    pbf.write('sim_dir = \'.\'\n')

    # store parameter files
    pbf.close()
    tar.add(p_base, arcname = '%s/params_base.py' % (prefix))
    os.remove(p_base)

    tar.close()
    print('Tar prepared at: %s' % (os.path.abspath(tp)))

def retrieve(src, bd):
    """
    Retrieve packaged plotting dependencies from `src` to `wd`.
    src: remote package file
        eg: fitzroy:/nesi/projects/nesi00213/scratch/temp.tar
        or: user@fitzroy:/nesi/projects/nesi00213/scratch/temp.tar
        same username on fitzroy as local is ideal
        ssh keys setup (no password entry) is ideal
    bd: base directory to extract mini sim_dir in
    """
    tar_download = '%s/%s' % (bd, os.path.basename(src))

    # retrieve tar file from fitzroy
    cmd = ['scp', src, tar_download]
    Popen(cmd).wait()
    if not os.path.exists(tar_download):
        print('Cannot find downloaded tar file.')
        exit(1)
    if not tf.is_tarfile(tar_download):
        print('Cannot read downloaded tar file.')
        exit(1)

    # extract tar file
    with tf.open(tar_download) as tar:
        tar.extractall(path = bd)

###
### Process Arguements
###
if len(sys.argv) < 4:
    print('USAGE: %s <prepare|retrieve|auto> <src> <dest>' % (sys.argv[0]))
    exit(1)
if sys.argv[1] == 'prepare':
    from srf import srf2corners

    if os.path.isdir(sys.argv[2]):
        prepare(os.path.abspath(sys.argv[2]), sys.argv[3])
    else:
        print('source directory doesn\'t exist!')
        exit(1)

elif sys.argv[1] == 'retrieve':
    # plotting base directory (like RunFolder which contains mini sim_dir)
    bd = os.path.abspath(sys.argv[3])
    if not os.path.isdir(bd):
        os.makedirs(bd)
    retrieve(sys.argv[2], bd)

elif sys.argv[1] == 'auto':
    # in case someone likes having different usernames everywhere
    if len(sys.argv) > 4:
        ruser = sys.argv[4]
    else:
        ruser = os.environ.get('USER')

    ###
    ### REMOTE PREPARE
    ###

    # location of scripts
    lbase = os.path.abspath(os.path.dirname(__file__))
    # remote temp dir
    # TODO: remove this hardcoded value and add to external config
    rtemp = '/nesi/projects/nesi00213/scratch/%s/data_transfers' % (ruser)
    # make sure user has temp directory
    Popen(['ssh', '%s@fitzroy.nesi.org.nz' % (ruser), \
            'mkdir -p %s' % (rtemp)]).wait()
    # everything is put in same location, no need for qcore_path.py
    Popen(['ssh', '%s@fitzroy.nesi.org.nz' % (ruser), \
            'touch %s/qcore_path.py' % (rtemp)]).wait()
    # copy scripts to remote temp
    Popen(['scp', __file__, '%s/../srf.py' % (lbase), \
            '%s/../geo.py' % (lbase),\
            '%s@fitzroy.nesi.org.nz:%s/' % (ruser, rtemp)]).wait()
    # run prepare in temp
    # TODO: change those hardcoded values to fitzroy
    Popen(['ssh', '%s@fitzroy.nesi.org.nz' % (ruser), \
            '/opt/niwa/Python/AIX/2.7.5/bin/python %s/%s prepare %s %s' \
            % (rtemp, os.path.basename(__file__), \
            sys.argv[2], '%s/data.tar' % (rtemp)) \
            ]).wait()

    ###
    ### LOCAL RETRIEVE
    ###

    bd = os.path.abspath(sys.argv[3])
    if not os.path.isdir(bd):
        os.makedirs(bd)
    retrieve('%s@fitzroy.nesi.org.nz:%s/data.tar' % (ruser, rtemp), bd)

    ###
    ### REMOVE TEMP
    ###

    # temporary directory on Fitzroy
    Popen(['ssh', '%s@fitzroy.nesi.org.nz' % (ruser), \
            'rm -r %s' % (rtemp)]).wait()
    # local tar file
    if os.path.exists('%s/data.tar' % (bd)):
        os.remove('%s/data.tar' % (bd))

    ###
    ### START PLOTS
    ###

    # local working directory
    wd = '%s/%s' % (bd, os.path.basename(sys.argv[2].rstrip('/')))
    # if params_plot.py doesn't exist, copy default
    if not os.path.exists('%s/params_plot.py' % (wd)):
        copyfile('%s/params_plot.template.py' % (lbase), \
                '%s/params_plot.py' % (wd))
    # start PGV/MMI plot
    Popen(['%s/plot_ts_sum.py' % (lbase)], cwd = wd).wait()
    Popen(['%s/plot_ts.py' % (lbase), '32'], cwd = wd).wait()
