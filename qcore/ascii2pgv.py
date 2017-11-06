#!/usr/bin/env python2
"""
First parameter is station file
Second parameter is Vel folder
"""
import os

import numpy as np

from timeseries import read_ascii

comps = ['090', '000', 'ver']

if __name__ == "__main__":
    import sys

    try:
        stat_file = sys.argv[1]
        vel_dir = sys.argv[2]
        assert(os.path.exists(stat_file))
        assert(os.path.exists(vel_dir))
    except IndexError:
        print("1st param is ll stat file and 2nd is vel dir")
        exit(1)
    except AssertionError:
        print("stat file or vel dir not found")
        exit(1)

    # load station locations
    stat_pos = {}
    with open(stat_file, 'r') as sf:
        for line in sf:
            lon, lat, stat = line.split()
            stat_pos[stat] = (float(lon), float(lat))

    # store PGVs for all stations
    pf = open('pgvs.txt', 'w')
    for stat, lonlat in stat_pos.items():
        try:
            # x, y, z
            data = (read_ascii('%s/%s.%s' % (vel_dir, stat, c)) \
                    for c in comps)
            # PGV = max(sqrt(x^2 + y^2 + z^2))
            pf.write('%f %f %f\n' % (lonlat[0], lonlat[1], \
                    max(np.sqrt(sum((d ** 2 for d in data))))))
        except IOError:
            print('Skipping station with missing files: %s/%s.xxx' % (vel_dir, stat))
            continue
        except KeyboardInterrupt:
            print('Exiting on KeyboardInterrupt.')
            pf.close()
            exit(1)
        except:
            pf.close()
            raise
