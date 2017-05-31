#!/usr/bin/env python2

from hashlib import sha1
import os
from shutil import rmtree

from scipy.misc import imread

###
### 1. Import GMT Library
###
print('Loading GMT Library...')
try:
    import gmt
except:
    print('Loading GMT Library Failed.')
    raise
print('GMT Library Loaded.')

GMT_PATHS = { \
    '5.1.2':'/opt/gmt-5.1.2/bin/gmt', \
    '5.3.2':'/opt/gmt-5.3.2/bin/gmt', \
    '5.3.3':'/opt/gmt-5.3.3/bin/gmt', \
    '5.4.1':'/opt/gmt-5.4.1/bin/gmt'}

test_dir = os.path.abspath('GMT_TESTING')
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

###
### LIST OF FUNCTIONS AND EXPECTED HASH RESULT
###
tests = { \
    test_coastlines : '61efdfbe4cd9bfd5a90ad86c67013c8d9494abc6', \
    test_land : '8c31ef6345aaad068cc4eb3e1d4b3f78cd6fc9a3'
}

###
### ALL TESTS ARE RUN FROM HERE
###
def run_tests():
    for test_func, test_hash in tests.items():
        for gmt_v in GMT_PATHS:
            iwd = os.path.join(test_dir, test_func.__name__, gmt_v)
            if not os.path.exists(iwd):
                os.makedirs(iwd)
            gmt.update_gmt_path(GMT_PATHS[gmt_v])

            pf = '%s/coastlines-%s.ps' % (iwd, gmt_v)
            test_func(pf)
            ph = sha1(imread('%s.png' % (os.path.splitext(pf)[0]))).hexdigest()

            if ph == test_hash:
                print('%s [%s] PASS' % (gmt_v, test_func.__name__))
            else:
                print('%s [%s] FAIL (%s)' % (gmt_v, test_func.__name__, ph))

if __name__ == '__main__':
    run_tests()
