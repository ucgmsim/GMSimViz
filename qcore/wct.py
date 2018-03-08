#!/usr/bin/env python2

import datetime
import os
import sys

import sqlite3 as lite

from qcore.config import qconfig

class WallClockDB:
    def __init__(self):
        con = self.connect()
        self.est = None
        with con:
            cur = con.cursor()
            cur.execute("SELECT max(estimator),avg(estimator),min(estimator) FROM wall_clock")
            self.est = list(cur.fetchone())

    def connect(self):
        # TODO: remove this hardcoded value here. See the GM slurm based repo for the same file
        con = lite.connect(qconfig['wallclock'])
        return con

    def estimate_wall_clock_time(self, nx, ny, nz, sim_duration,num_procs=512):
        prod = 1.0* nx*ny*nz*sim_duration*512.0/num_procs

        wall_clock_est = [x * prod for x in self.est]
        #print wall_clock_est

        est = [str(datetime.timedelta(seconds=x)) for x in wall_clock_est]
        print "nx=%d ny=%d nz=%d sim_duration=%d num_procs=%d" %(nx,ny,nz,sim_duration,num_procs)
        print "Maximum: %s" %est[0]
        print "Average: %s" %est[1]
        print "Minimum: %s" %est[2]

        return est #max, avg, min

    def add(self, path, nx,ny,nz,sim_duration,num_procs,start,end):
        start = datetime.datetime.strptime(start,"%H:%M:%S_%d/%m/%Y")
        end = datetime.datetime.strptime(end,"%H:%M:%S_%d/%m/%Y")
        wall_clock = end-start        
        wall_clock_sec = wall_clock.seconds
        print "Wall clock time used:", wall_clock
        con = self.connect()
        values = ( (path,nx,ny,nz,sim_duration,wall_clock_sec,float(wall_clock_sec*num_procs)/(nx*ny*nz*sim_duration*512.0),num_procs))
        with con:
            cur = con.cursor()
            cur.execute("INSERT INTO wall_clock VALUES (?,?,?,?,?,?,?,?)", values)

def usage():
    print "Usage: %s -q : Query for wall-clock estimation. Uses params.py in the current path" %sys.argv[0]
    print "Usage: %s -q nx ny nz sim_duration num_procs: Query for wall-clock estimation with specified parameters" %sys.argv[0]
    print "Usage: %s -a start end : Add the measured wall-clock time to db. params.py must be in the current path"%sys.argv[0]
    print "        start, end : Must be in %H:%M:%S_%d/%m/%Y format. Try `date +%H:%M:%S_%d/%m/%Y`"
    sys.exit()

if __name__ == '__main__':
    db = WallClockDB()
    print db.est

    if len(sys.argv)==1 or sys.argv[1]=='-h':
        usage()

    if sys.argv[1]=='-q':
        if len(sys.argv)==7:
            nx = int(sys.argv[2])
            ny = int(sys.argv[3])
            nz = int(sys.argv[4])
            sim_duration = int(sys.argv[5])
            num_procs = int(sys.argv[6])
            db.estimate_wall_clock_time(nx,ny,nz,sim_duration,num_procs)
        elif len(sys.argv)==2:
            curdir = os.curdir
            sys.path.append(os.curdir)
            try:
                import params
            except ImportError:
                print "Error: Make sure the current directory has params.py"
                sys.exit()

            nx = int(params.nx)
            ny = int(params.ny)
            nz = int(params.nz)
            sim_duration = int(float(params.sim_duration))
            num_procs = int(params.n_proc)
            db.estimate_wall_clock_time(nx,ny,nz,sim_duration,num_procs)

        else:
            print "Check arguments:", sys.argv
            usage()

    elif sys.argv[1]=='-a':
        curdir = os.curdir
        sys.path.append(os.curdir)
        try:
            import params
        except ImportError:
            print "Error: Make sure the current directory has params.py"
            sys.exit()

        if len(sys.argv)==4:
            path = os.path.abspath(curdir)
            start = sys.argv[2]
            end = sys.argv[3]
            nx = int(params.nx)
            ny = int(params.ny)
            nz = int(params.nz)
            sim_duration = int(float(params.sim_duration))
            num_procs = int(params.n_proc)
            db.add(path, nx,ny,nz,sim_duration,num_procs,start,end)
        else:
            print "Check arguments:", sys.argv
            usage()

