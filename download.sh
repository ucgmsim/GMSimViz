#!/bin/bash
mydir=`dirname $0
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 local_dir fitzroy_dir fitzroy_username"
    echo "   Eg. $0 /home/seb56 /a/b/x baes"
    echo "       This program creates /home/seb/x and copies everything under x recursively"
    exit 1
fi

mpirun -np 8 python2 $mydir/parallel_download.py $1 $2 $3


