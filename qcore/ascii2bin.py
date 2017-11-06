#!/usr/bin/env python2
"""
Convert 3 ASCII component files into binary versions.
Works on a single set of files by specifying prefix, or
Works on all files in a directory.

Expects .090, .000 and .ver extensions.

If second parameter exists it's the output folder.
"""

from glob import glob
import os
import sys

import numpy as np

import timeseries as ts

# x, y, z
suffixes = ['.090', '.000', '.ver']

try:
    path = os.path.abspath(sys.argv[1])
except IndexError:
    print('First parameter is the prefix or path.')
    exit(1)

def convert_station(station_prefix):
    """
    convert 3 component files (ascii) to binary LE SP float
    """
    a = []
    for s in suffixes:
        a.append(ts.read_ascii('%s%s' % (station_prefix, s)))
    if len(sys.argv) > 2:
        station_prefix = '%s/%s' \
                % (sys.argv[2], os.path.basename(station_prefix))
    np.array(a).T.astype('<f4').tofile('%s.bin' % (station_prefix))


# input is for a single station
if os.path.exists('%s.090' % (path)):
    convert_station(path)

# input is a whole folder
elif os.path.isdir(path):
    for c in glob('%s/*%s' % (path, suffixes[0])):
        try:
            convert_station(c[:-len(suffixes[0])])
        except (ValueError, TypeError, IOError):
            print('Failed converting station %s.' % (c[:-len(suffixes[0])]))

else:
    print('Artificial inteligence does not yet understand: %s' % (path))
