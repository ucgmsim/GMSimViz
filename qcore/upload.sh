#!/bin/bash
mydir=`dirname $0`
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 local_dir fitzroy_dir fitzroy_username"
    echo "   Eg. $0 ./x /home/baes baes"
    echo "       where abspath of x is /a/b/x"
    echo "       This program creates /home/baes/x and copies everything under x recursively"
    exit 1
fi
mpirun -np 4 python2 $mydir/parallel_upload.py $1 $2 $3


