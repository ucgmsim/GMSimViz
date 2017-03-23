#!/bin/bash
mydir=`dirname $0`
mpirun -np 4 python2 $mydir/parallel_upload.py $1 $2 $3


