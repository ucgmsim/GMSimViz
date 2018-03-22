"""
Miscellaneous non-specific functions.
Module which contains shared functions/values.

@date 8 April 2016
"""

from __future__ import print_function

import subprocess
import sys

# returns a list of stations
# sample line in source file:
#   171.74765   -43.90236 ADCS
def get_stations(source_file, locations = False):
    stations = []
    station_lats = []
    station_lons = []
    with open(source_file, 'r') as sp:
        for line in sp.readlines():
            if line[0] not in  ['#', '%']:
                info = line.split()
                stations.append(info[2])
                if locations:
                    station_lons.append(info[0])
                    station_lats.append(info[1])
    if not locations:
        return stations
    return (stations, station_lats, station_lons)

def get_corners(model_params, gmt_format = False):
    """
    Retrieve corners of simulation domain from model params file.
    model_params: file path to model params
    gmt_format: if True, also returns corners in GMT string format
    """
    # with -45 degree rotation:
    #   c2
    # c1  c3
    #   c4
    corners = []
    with open(model_params, 'r') as vmpf:
        lines = vmpf.readlines()
        # make sure they are read in the correct order at efficiency cost
        for corner in ['c1=', 'c2=', 'c3=', 'c4=']:
            for line in lines:
                if corner in line:
                    corners.append(map(float, line.split()[1:3]))
                    break
    if not gmt_format:
        return corners
    # corners in GMT format
    cnr_str = '\n'.join([' '.join(map(str, cnr)) for cnr in corners])
    return corners, cnr_str

def exe(cmd, debug = True, shell = False, \
        stdout = True, stderr = True, stdin = None):
    """
    cmd: command as list starting with executable, followed by arguments.
         Strings will be split by whitespace even if this splits a parameter.
         This will cause issues when shell == False. List input is ideal.
    debug: print equivalent shell command, display output
    shell: execute command in shell environment (not recommended)
    stdout: True: return output | file: open file object
    stderr: True: return error | file: open file object
    stdin: None for no input or command input string
    """

    # always split for consistency
    if type(cmd) == str:
        cmd = cmd.split(' ')

    # display what command would look like if executed on a shell
    if debug:
        virtual_cmd = ' '.join(cmd)
        if type(stdout) == file:
            virtual_cmd = '%s 1>%s' % (virtual_cmd, stdout.name)
        if type(stderr) == file:
            virtual_cmd = '%s 2>%s' % (virtual_cmd, stderr.name)
        print(virtual_cmd, file = sys.stderr)

    # special cases for stderr and stdout
    if stdout == True:
        stdout = subprocess.PIPE
    if stderr == True:
        stderr = subprocess.PIPE

    p = subprocess.Popen(cmd, shell = shell, stdout = stdout, stderr = stderr)
    out, err = p.communicate(stdin)
    rc = p.wait()

    if debug:
        if out:
            print(out, file = sys.stderr)
        if err:
            print(err, file = sys.stderr)

    return out, err
