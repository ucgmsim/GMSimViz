#!/usr/bin/env python2
import os
from subprocess import check_call
import sys
import os.path

# should already be in PYTHONPATH (qcore library)
from shared import *

GEN_CORD_BIN = "/nesi/projects/nesi00213/tools/gen_model_cords"
# hard coded defaults
CENTER_ORIGIN = '1'
GEOPROJ = '1'
DOCOORDS = '1'

# leave main even though it makes no sense unless you change all code calling this
def main(outdir = '.', debug = False):

    # load params for velocity model
    sys.path.insert(0, outdir)
    sys.path.append('.')
    try:
        import params_vel as vm
    except ImportError:
        try:
            import params_base as vm
        except ImportError:
            import params as vm

    # these are saved as strings, should be floats because that's what they are
    # YLEN == vm.extent_y etc. in params_vel but not elsewhere
    XLEN = float(vm.nx) * float(vm.hh)
    YLEN = float(vm.ny) * float(vm.hh)
    ZLEN = float(vm.nz) * float(vm.hh)

    GRIDFILE = os.path.join(outdir, 'gridfile%s' % (vm.sufx))
    GRIDOUT = os.path.join(outdir, 'gridout%s' % (vm.sufx))
    MODEL_COORDS = os.path.join(outdir, 'model_coords%s' % (vm.sufx))
    MODELPARAMS = os.path.join(outdir, 'model_params%s' % (vm.sufx))
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
        print('Cannot write to location: %s' % (GRIDFILE))
        sys.exit(1)

    # generate model_params
    cmd = ("{GEN_CORD_BIN} "
            "geoproj={GEOPROJ} gridfile='{GRIDFILE}' gridout='{GRIDOUT}' "
            "center_origin={CENTER_ORIGIN} do_coords={DOCOORDS} "
            "nzout=1 name='{MODEL_COORDS}' gzip=0 latfirst=0 "
            "modellon={vm.MODEL_LON} modellat={vm.MODEL_LAT} "
            "modelrot={vm.MODEL_ROT} 1> '{MODELPARAMS}'") \
            .format(**dict(locals(), **globals()))
    if not debug:
        cmd += ' 2>/dev/null'
    print(cmd)
    check_call(cmd, shell = True)

    # also generate coordinate related outputs
    if DOCOORDS != '1':
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
        print('Could not write bounds file.')
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(outdir = sys.argv[1])
    else:
        main()
