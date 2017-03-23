#!/bin/bash
mydir=`dirname $0`
mpirun -np 8 python2 $mydir/parallel_download.py $1 $2 $3


