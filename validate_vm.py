#!/usr/bin/env python
"""
Checks VM folder for corectness. Compatible with Python 2.6+ and 3.0+.

Run script with VM folder location as first parameter. Returns 0 if successful.
or:
Import and use validate_vm directly.
"""

from __future__ import print_function

import os
import sys

try:
    import numpy as np
    numpy = True
except ImportError:
    numpy = False

def validate_vm(vm_dir, verbose = False, errors = True):
    """
    Go through rules of VM directories. Return False if invalid.
    vm_dir: folder path containing VM files
    verbose: print progress
    errors: print errors and warnings to stderr
    """
    SIZE_FLOAT = 4

    def vmfile(filename):
        return os.path.join(vm_dir, filename)

    def vprint(msg):
        if verbose:
            print(msg)
    def eprint(msg):
        if errors:
            print(msg, file = sys.stderr)

    vprint('Validating VM \'%s\'...' % (vm_dir))

    # 1: has to exist
    if not os.path.isdir(vm_dir):
        eprint('VM dir %s is not a directory.' % (vm_dir))
        return False

    # 2: fixed file names exist
    vm = {'s':'%s' % (vmfile('vs3dfile.s')), \
            'p':'%s' % (vmfile('vp3dfile.p')), \
            'd':'%s' % (vmfile('rho3dfile.d'))}
    for fixed_name in vm.values():
        if not os.path.exists(fixed_name):
            eprint('VM file not found: %s' % (fixed_file))
            return False
    if not os.path.exists(vmfile('params_vel.py')):
        eprint('VM configuration params_vel.py missing.')
        return False

    # 3: metadata files exist (made by gen_cords.py)
    sys.path.insert(0, vm_dir)
    import params_vel as vm_conf
    meta = {'gridfile':'%s' % (vmfile('gridfile%s' % (vm_conf.sufx))), \
            'gridout':'%s' % (vmfile('gridout%s' % (vm_conf.sufx))), \
            'bounds':'%s' % (vmfile('model_bounds%s' % (vm_conf.sufx))), \
            'coords':'%s' % (vmfile('model_coords%s' % (vm_conf.sufx))), \
            'params':'%s' % (vmfile('model_params%s' % (vm_conf.sufx)))}
    meta_created = True
    for meta_file in meta.values():
        if not os.path.exists(meta_file):
            eprint('WARNING: VM metadata not found: %s' % (meta_file))
            meta_created = False

    # 4: params_vel.py consistency
    try:
        nx, ny, nz = map(int, [vm_conf.nx, vm_conf.ny, vm_conf.nz])
        xlen, ylen, zmin, zmax, hh = map(float, \
                [vm_conf.extent_x, vm_conf.extent_y, \
                vm_conf.extent_zmin, vm_conf.extent_zmax, vm_conf.hh])
    except AttributeError:
        eprint('VM params_vel.py missing values.')
        return False
    except ValueError:
        eprint('VM params_vel.py contains invalid numbers.')
        return False
    zlen = zmax - zmin
    try:
        assert(nx == int(round(xlen / hh)))
        assert(ny == int(round(ylen / hh)))
        assert(nz == int(round(zlen / hh)))
    except AssertionError:
        eprint('VM params_vel.py missmatch between extents and nx, ny, nz.')
        return False

    # 5: binary file sizes
    vm_size = nx * ny * nz * SIZE_FLOAT
    for bin_file in vm.values():
        size = os.path.getsize(bin_file)
        if size != vm_size:
            eprint('VM filesize for %s expected: %d found: %d' \
                    % (bin_file, vm_size, size))
            return False

    # 6: binary contents
    if numpy:
        # check first zx slice (y = 0)
        smin = np.min(np.fromfile(vm['s'], dtype = '<f%d' % (SIZE_FLOAT), \
                count = nz * nx))
        pmin = np.min(np.fromfile(vm['p'], dtype = '<f%d' % (SIZE_FLOAT), \
                count = nz * nx))
        dmin = np.min(np.fromfile(vm['d'], dtype = '<f%d' % (SIZE_FLOAT), \
                count = nz * nx))
        # works even if min is np.nan
        if not min(smin, pmin, dmin) > 0:
            eprint('VM vs, vp or rho <= 0|nan found.')
            return False
    else:
        eprint('WARNING: VM check value sanity not completed, numpy missing.')

    # 7: contents of meta files
    if meta_created:
        # TODO: check individual file contents
        # not as important, can be re-created based on params_vel.py
        pass

    vprint('Validating VM success.')
    return True

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if validate_vm(sys.argv[1]):
            sys.exit(0)
    sys.exit(1)
