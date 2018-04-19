#!/usr/bin/env python2

import datetime
import os
import sys

import sqlite3

from qcore.config import qconfig

class WallClockDB:

    def __init__(self, db_path = qconfig['wallclock']):
        self.db_path = db_path

    def max_avg_min(self):
        """
        Retrieve maximum, average and minimum estimator factor.
        """
        db = sqlite3.connect(self.db_path)
        try:
            cursor = db.cursor()
            cursor.execute('''SELECT
                              max(estimator), avg(estimator), min(estimator)
                              FROM wall_clock''')
            return cursor.fetchone()
        except Exception:
            raise
        finally:
            db.close()

    def insert(self, values):
        """
        Insert new estimator factor.
        values: (path, nx, ny, nz, sim_duration, wall_clock_sec,
                estimator_factor, num_procs)
        """
        db = sqlite3.connect(self.db_path)
        try:
            cursor = db.cursor()
            cursor.execute('''INSERT INTO wall_clock
                              VALUES (?,?,?,?,?,?,?,?)''', values)
            db.commit()
        except Exception:
            db.rollback()
            raise
        finally:
            db.close()

    def estimate_wall_clock_time(self, nx, ny, nz, sim_duration, \
            num_procs=512, string = True):
        """
        Display, return max, avg, and min wall clock estimate for given params.
        """
        computation = nx * ny * nz * sim_duration * 512.0 / num_procs
        wall_clock_sec = [x * computation for x in self.max_avg_min]

        if string:
            return [str(datetime.timedelta(seconds = x)) \
                    for x in wall_clock_sec]
        return wall_clock_sec

    def add(self, path, nx, ny, nz, sim_duration, num_procs, start, end):
        """
        Add a wall clock time estimation factor given times as string.
        """
        wall_clock = datetime.datetime.strptime(end, '%H:%M:%S_%d/%m/%Y') \
                - datetime.datetime.strptime(start, '%H:%M:%S_%d/%m/%Y')
        wall_clock_sec = wall_clock.seconds + wall_clock.days * 86400

        values = (path, nx, ny, nz, sim_duration, wall_clock_sec, \
                float(wall_clock_sec * num_procs) \
                    / (nx * ny * nz * sim_duration * 512.0), \
                num_procs)
        self.insert(values)

def print_query(sec_str, nx, ny, nz, sim_duration, num_procs):
    print('nx=%d ny=%d nz=%d sim_duration=%d num_procs=%d' \
            % (nx, ny, nz, sim_duration, num_procs))
    print('Minimum: %s' % sec_str[2])
    print('Average: %s' % sec_str[1])
    print('Maximum: %s' % sec_str[0])

def usage():
    print('''
%s -q
    Query wall-clock estimation using params.py in path
%s -q nx ny nz sim_duration num_procs
    Query wall-clock estimation with specified parameters
%s -a start_time end_time
    Add wall-clock factor to db using params.py in path and given times

start_time, end_time must be in %%H:%%M:%%S_%%d/%%m/%%Y format.
Try `date +%%H:%%M:%%S_%%d/%%m/%%Y` for an example.''' \
        % (sys.argv[0], sys.argv[0], sys.argv[0]))

if __name__ == '__main__':
    db = WallClockDB()

    # cannot complete any further functions
    if len(sys.argv) == 1 or sys.argv[1] == '-h':
        usage()
        sys.exit(1)

    # query database
    if sys.argv[1] == '-q':

        # all parameters provided
        if len(sys.argv) == 7:
            # nx ny nz sim_duration num_procs
            factors = (int(sys.argv[2]), \
                    int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), \
                    int(sys.argv[6]))

        # parameters must be read from params.py
        elif len(sys.argv) == 2:
            sys.path.insert(0, os.curdir)
            try:
                import params as p
            except ImportError:
                print('params.py not found')
                sys.exit(1)
            # nx ny nz sim_duration num_procs
            factors = (int(p.nx), int(p.ny), \
                    int(p.nz), int(float(p.sim_duration)), int(p.n_proc))

        # invalid number of parameters
        else:
            usage()
            sys.exit(1)

        # run and display
        result = db.estimate_wall_clock_time(*factors)
        print_query(result, *factors)

    # add entry to database
    elif sys.argv[1] == '-a':

        # required arguments
        if len(sys.argv) != 4:
            usage()
            sys.exit(1)

        # required parameters
        sys.path.insert(0, os.curdir)
        try:
            import params as p
        except ImportError:
            print('params.py not found')
            sys.exit(1)

        # path, nx, ny, nz, sim_duration, num_procs, start, end
        db.add(os.path.abspath(os.curdir), int(p.nx), int(p.ny), int(p.nz), \
                int(float(p.sim_duration)), int(p.n_proc), \
                sys.argv[2], sys.argv[3])
