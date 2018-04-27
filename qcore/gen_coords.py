#!/usr/bin/env python2

import os
from subprocess import check_call
import sys

from qcore.config import qconfig

GEN_COORD_BIN = os.path.join(qconfig['tools_dir'], 'gen_model_cords')

def gen_coords(outdir = '.', debug = False, geoproj = '1', do_coords = '1', \
        centre_origin = '1'):
    """
    Generate coordinate files for an emod3d domain (set of grid points).
    outdir: directory to store coordinate files in
    debug: print additional info
    geoproj: 
    do_coords: 
    """

    # load params for velocity model
    sys.path.insert(0, '.')
    sys.path.insert(0, outdir)
    # parameters file historically not consistent
    try:
        import params_vel as vm
    except ImportError:
        try:
            import params_base as vm
        except ImportError:
            import params as vm
    # important to notify of location to prevent path order issues
    print('sim params: %s' % (os.path.abspath(vm.__file__)))

    # these are saved as strings, should be floats because that's what they are
    # YLEN == vm.extent_y etc. in params_vel but not elsewhere
    XLEN = float(vm.nx) * float(vm.hh)
    YLEN = float(vm.ny) * float(vm.hh)
    ZLEN = float(vm.nz) * float(vm.hh)

    # list of outputs that this function can create
    GRIDFILE = os.path.join(outdir, 'gridfile%s' % (vm.sufx))
    GRIDOUT = os.path.join(outdir, 'gridout%s' % (vm.sufx))
    MODEL_COORDS = os.path.join(outdir, 'model_coords%s' % (vm.sufx))
    MODEL_PARAMS = os.path.join(outdir, 'model_params%s' % (vm.sufx))
    MODEL_BOUNDS = os.path.join(outdir, 'model_bounds%s' % (vm.sufx))

    # generate gridfile
    try:
        with open(GRIDFILE, 'w') as gridf:
            gridf.write("xlen=%f\n" % (XLEN))
            gridf.write("%10.4f %10.4f %13.6e\n" % (0.0, XLEN, float(vm.hh)))
            gridf.write("ylen=%f\n" % (YLEN))
            gridf.write("%10.4f %10.4f %13.6e\n" % (0.0, YLEN, float(vm.hh)))
            gridf.write("zlen=%f\n" % (ZLEN))
            gridf.write("%10.4f %10.4f %13.6e\n" % (0.0, ZLEN, float(vm.hh)))
    except IOError:
        raise IOError('Cannot write GRIDFILE: %s' % (GRIDFILE))

    # generate model_params
    cmd = ("{GEN_COORD_BIN} "
            "geoproj={geoproj} gridfile='{GRIDFILE}' gridout='{GRIDOUT}' "
            "center_origin={centre_origin} do_coords={do_coords} "
            "nzout=1 name='{MODEL_COORDS}' gzip=0 latfirst=0 "
            "modellon={vm.MODEL_LON} modellat={vm.MODEL_LAT} "
            "modelrot={vm.MODEL_ROT} 1> '{MODEL_PARAMS}'") \
            .format(**dict(locals(), **globals()))
    if debug:
        print(cmd)
    else:
        cmd += ' 2>/dev/null'
    check_call(cmd, shell = True)

    # also generate coordinate related outputs
    if do_coords != '1':
        return

    # retrieve MODEL_BOUNDS
    x_bounds = [0, float(vm.nx) - 1]
    y_bounds = [0, float(vm.ny) - 1]
    try:
        with open(MODEL_COORDS, 'r') as coordf:
            with open(MODEL_BOUNDS, 'w') as boundf:
                for line in coordf:
                    x, y = map(float, line.split()[2:4])
                    if x in x_bounds or y in y_bounds:
                        boundf.write(line)
    except IOError:
        raise IOError('Cannot write MODEL_BOUNDS: %s' % (MODEL_BOUNDS))


# allow running from shell
if __name__ == "__main__":
    if len(sys.argv) > 1:
        gen_coords(outdir = sys.argv[1])
    else:
        gen_coords()

