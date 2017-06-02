#!/usr/bin/env python2
"""
MPI pattern based on
github.com/jbornschein/mpi4py-examples/blob/master/09-task-pull.py
"""

from hashlib import sha1
import os
from shutil import rmtree
from time import time

from scipy.misc import imread
from mpi4py import MPI

# first test is loading library
import gmt

# MPI stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
status = MPI.Status()
hostname = MPI.Get_processor_name()
MASTER = 0
class tags:
    READY = 15
    DONE = 85
    START = 204
    EXIT = 240

# GMT versions to test
GMT_PATHS = { \
    '5.1.2':'/opt/gmt-5.1.2/bin/gmt', \
    '5.3.2':'/opt/gmt-5.3.2/bin/gmt', \
    '5.3.3':'/opt/gmt-5.3.3/bin/gmt', \
    '5.4.1':'/opt/gmt-5.4.1/bin/gmt'}

test_dir = os.path.abspath('GMT_TESTING')
if rank == MASTER:
    if os.path.exists(test_dir):
        rmtree(test_dir)
    os.makedirs(test_dir)

###
### TESTING FUNCTIONS
###

def test_coastlines(pf):
    p = gmt.GMTPlot(pf)
    p.spacial('M', (170, 175, -44, -38), sizing = 3)
    p.coastlines(width = '0.4p')
    p.finalise()
    p.png(dpi = 200, clip = True)

def test_land(pf):
    p = gmt.GMTPlot(pf)
    p.spacial('M', (170.1, 179.91, -37, -34), sizing = 7)
    p.land(fill = 'darkred')
    p.finalise()
    p.png(dpi = 100, clip = True)

def test_ticks(pf):
    p = gmt.GMTPlot(pf)
    p.spacial('M', (160.992, 174.9122, -44, -34.01), sizing = 2)
    p.ticks(major = '1d', minor = '20m', sides = 'ew')
    p.finalise()
    p.png(dpi = 100, clip = True)

def test_ticks2(pf):
    p = gmt.GMTPlot(pf)
    p.spacial('T', (160.992, 174.9122, -44, -34.01), sizing = 5, \
            x_shift = 0.5, y_shift = 0.5, lon0 = 174.9122)
    p.ticks(major = '1d', minor = '20m', sides = 'ew')
    p.finalise()
    p.png(dpi = 100, clip = True)

def test_cpt(pf):
    cptf = '%s.cpt' % (os.path.splitext(pf)[0])
    gmt.makecpt('hot', cptf, 0, 120, inc = 0.1, invert = True, \
            wd = os.path.dirname(pf))
    p = gmt.GMTPlot(pf)
    p.spacial('X', (0, 15, 0, 4), sizing = '15/4')
    p.cpt_scale(6.1, 2.05, cptf, 20, 5, label = 'test_scale', \
            length = 3.05, thickness = '0.3i')
    p.finalise()
    p.png(dpi = 320, clip = True)

def test_cpt2(pf):
    cptf = '%s.cpt' % (os.path.splitext(pf)[0])
    gmt.makecpt('polar', cptf, -1.5, 1.5, inc = 0.25, invert = False, \
            wd = os.path.dirname(pf), bg = '0/0/80', fg = '80/0/0')
    p = gmt.GMTPlot(pf)
    p.spacial('X', (0, 4, 0, 2), sizing = '4/2', x_shift = 1, y_shift = 1)
    p.cpt_scale(0, 0, cptf, 0.5, 0.25, cross_tick = 0.5, align = 'LB', \
            length = 3, thickness = '0.3i', arrow_f = True, arrow_b = True)
    p.finalise()
    p.png(dpi = 222, clip = True)

def test_fill(pf):
    p = gmt.GMTPlot(pf)
    gmt.gmt_defaults(wd = os.path.dirname(pf), ps_media = 'A5')
    p.background(1.5, 1)
    p.background(3, 2, x_margin = 1.5, colour = 'blue')
    p.background(1.5, 1, y_margin = 1, colour = 'red')
    p.background(3, 2, x_margin = 4.5, colour = 'firebrick')
    p.finalise()
    p.png(dpi = 100, clip = False)

###
### LIST OF FUNCTIONS, EXPECTED HASH RESULT, MINIMUM VERSION
###
TESTS = ( \
    (test_coastlines, '61efdfbe4cd9bfd5a90ad86c67013c8d9494abc6', 5.0), \
    (test_land, '8c31ef6345aaad068cc4eb3e1d4b3f78cd6fc9a3', 5.0), \
    (test_ticks, 'd625bfb10464664d38f60537686999945eafafe7', 5.0), \
    (test_ticks2, '903830d2ba5d19a415a6d4e15fa62175bd7b7242', 5.0), \
    (test_cpt, '960032a10ef5fef8f013aab8892082fc4a606c9c', 5.2), \
    (test_cpt2, 'd83a3f811e8e514a6e88dacf38fe24bc4dea3280', 5.2), \
    (test_fill, 'eaa069d015eefb20f9826061a6de09d17a420a91', 5.0)
)

###
### ALL TESTS ARE RUN FROM HERE
###
def run_test(test, gmt_version):
    iwd = os.path.join(test_dir, test[0].__name__, gmt_version)
    if not os.path.exists(iwd):
        os.makedirs(iwd)
    gmt.update_gmt_path(GMT_PATHS[gmt_version])

    pf = '%s/%s-%s.ps' % (iwd, test[0].__name__, gmt_version)
    t0 = time()
    test[0](pf)
    t = time() - t0
    try:
        pp = '%s.png' % (os.path.splitext(pf)[0])
        ph = sha1(imread(pp)).hexdigest()
        os.symlink(pp, os.path.join(test_dir, os.path.basename(pp)))
        os.symlink(pp, os.path.join(test_dir, \
                test[0].__name__, os.path.basename(pp)))
    except IOError:
        print('%s [%s] FAIL (NO OUTPUT) %.2fs' \
                % (gmt_version, test[0].__name__, t))
        return False

    if ph == test[1]:
        print('%s [%s] PASS %.2fs' \
                % (gmt_version, test[0].__name__, t))
        return True
    print('%s [%s] FAIL (%s) %.2fs' \
            % (gmt_version, test[0].__name__, ph, t))
    return False

###
### DISTRIBUTE TESTS VIA MPI PROCESSES
###
if rank == MASTER:
    # DISTRIBUTE WORK
    gmt_versions = list(GMT_PATHS)
    jobs = len(GMT_PATHS) * len(TESTS)
    jobs_run = 0
    job = 0
    passed = 0

    workers = size - 1
    workers_closed = 0

    while workers_closed < workers:
        data = comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, \
                status = status)
        source = status.Get_source()
        tag = status.Get_tag()

        if tag == tags.READY:
            found_job = False
            while job < jobs:
                workload = (TESTS[job % len(TESTS)], \
                        gmt_versions[job // len(TESTS)])
                # check that this test is compatible with this version
                major_v = float('.'.join(workload[1].split('.')[:2]))
                if major_v >= workload[0][2]:
                    found_job = True
                    comm.send(workload, dest = source, tag = tags.START)
                    jobs_run += 1
                    job += 1
                    break
                # incompatible gmt version / test combination. try next one
                job += 1
            if not found_job:
                # no more valid combinations
                comm.send(None, dest = source, tag = tags.EXIT)
        elif tag == tags.DONE:
            passed += data
        elif tag == tags.EXIT:
            workers_closed += 1

    print('=============================')
    print('TESTS: %d' % (jobs_run))
    print('PASSED: %d' % (passed))
    print('FAILED: %d' % (jobs_run - passed))

    if jobs == passed:
        rmtree(test_dir)

else:
    # ASK FOR WORK
    while True:
        comm.send(None, dest = MASTER, tag = tags.READY)
        task = comm.recv(source = MASTER, tag = MPI.ANY_TAG, status = status)
        tag = status.Get_tag()

        if tag == tags.START:
            result = run_test(task[0], task[1])
            comm.send(result, dest = MASTER, tag = tags.DONE)
        elif tag == tags.EXIT:
            break
    comm.send(None, dest = MASTER, tag = tags.EXIT)
