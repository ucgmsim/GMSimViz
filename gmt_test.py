#!/usr/bin/env python2

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

###
### LIST OF FUNCTIONS AND EXPECTED HASH RESULT
###
TESTS = ( \
    (test_coastlines, '61efdfbe4cd9bfd5a90ad86c67013c8d9494abc6'), \
    (test_land, '8c31ef6345aaad068cc4eb3e1d4b3f78cd6fc9a3'), \
    (test_ticks, 'd625bfb10464664d38f60537686999945eafafe7')
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
    ph = sha1(imread('%s.png' % (os.path.splitext(pf)[0]))).hexdigest()

    if ph == test[1]:
        print('%s [%s] PASS %.2fs' % (gmt_version, test[0].__name__, t))
        return True
    print('%s [%s] FAIL (%s) %.2fs' % (gmt_v, test[0].__name__, ph, t))
    return False

###
### DISTRIBUTE TESTS VIA MPI PROCESSES
###
if rank == MASTER:
    # DISTRIBUTE WORK
    gmt_versions = list(GMT_PATHS)
    jobs = len(GMT_PATHS) * len(TESTS)
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
            if job < jobs:
                workload = (TESTS[job % len(TESTS)], \
                        gmt_versions[job % len(gmt_versions)])
                comm.send(workload, dest = source, tag = tags.START)
                job += 1
            else:
                comm.send(None, dest = source, tag = tags.EXIT)
        elif tag == tags.DONE:
            passed += data
        elif tag == tags.EXIT:
            workers_closed += 1

    print('=============================')
    print('TESTS: %d' % (jobs))
    print('PASSED: %d' % (passed))
    print('FAILED: %d' % (jobs - passed))

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
