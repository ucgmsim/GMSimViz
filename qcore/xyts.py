"""
XYTS file related functions.
XYTS files contain time slices on the X-Y plane (z = 1, top level).
XYTS file processing

@author Viktor Polak
@date 19 January 2017
"""

from math import radians, cos, sin

import numpy as np

from geo import gp2ll_multi

###
### PROCESSING OF XYTS FILE
###
class XYTSFile:
    """
    Assumptions:
    dip = 0: simulation domain is flat
    t0 = 0: complete timeseries from t = 0
    """

    def __init__(self, xyts_path, meta_only = False):
        """
        Prepare file for reading.
        C structs available in WccFormat/src/structure.h
        Note: endianness is only important here (__init__)
                hardcoded big src '>', little would be '<'
                could check sanity in header for auto-endianness
                outputs are always native
        xyts_path: path to the xyts file
        meta_only: only using this to retrieve metadata
                doesn't enable timeslice capability
        """

        # read header
        xytf = open(xyts_path, 'rb')
        self.x0, self.y0, self.z0, self.t0, \
                self.nx, self.ny, self.nz, self.nt = \
                np.fromfile(xytf, dtype = '>i4', count = 8)
        self.dx, self.dy, self.hh, self.dt, \
                self.mrot, self.mlat, self.mlon = \
                np.fromfile(xytf, dtype = '>f4', count = 7)
        xytf.close()
        if self.nz != 1:
            print('Not an XY timeslice file.')
            print('Will not continue as data has different dimentions.')
            exit(1)

        # determine original sim parameters
        self.dxts = int(round(self.dx / self.hh))
        self.dyts = int(round(self.dy / self.hh))
        self.nx_sim = self.nx * self.dxts
        self.ny_sim = self.ny * self.dyts

        # orientation of components
        self.dip = 0
        self.comps = {'X':radians(90 + self.mrot), \
                'Y':radians(self.mrot), \
                'Z':radians(90 - self.dip)}
        # rotation of components so Y is true north
        self.cosR = cos(self.comps['X'])
        self.sinR = sin(self.comps['X'])
        # simulation plane always flat, dip = 0
        self.cosP = 0 # cos(self.comps['Z'])
        self.sinP = 1 # sin(self.comps['Z'])
        # xy dual component rotation matrix
        # must also flip vertical axis
        theta = radians(self.mrot)
        self.rot_matrix = np.array([[cos(theta), -sin(theta), 0], \
                [-sin(theta), -cos(theta), 0], \
                [0, 0, -1]])

        # save speed when only loaded to read metadata section
        if meta_only:
            return

        # memory map for data section
        self.data = np.memmap(xyts_path, dtype = '>f4', mode = 'r', \
                offset = 60, \
                shape = (self.nt, len(self.comps), self.ny, self.nx))

        # create longitude, latitude map for data
        xy_points = np.mgrid[0:self.nx_sim:self.dxts, \
                0:self.ny_sim:self.dyts] \
                .reshape(2, -1, order = 'F').T.tolist()
        ll_points = gp2ll_multi(xy_points, self.mlat, self.mlon, self.mrot, \
                self.nx_sim, self.ny_sim, self.hh)
        self.ll_map = np.array(ll_points).reshape(self.ny, self.nx, 2)

    def corners(self, gmt_format = False):
        """
        Retrieve corners of simulation domain.
        gmt_format: if True, also returns corners in GMT string format
        """
        # compared with model_params format:
        # c1 =   x0   y0
        # c2 = xmax   y0
        # c3 = xmax ymax
        # c4 =   x0 ymax
        # cannot just use self.ll_map as xmax, ymax for simulation domain
        # may have been decimated. sim nx 1400 (xmax 1399) with dxts 5 = 1395
        xy_cnrs = [[0, 0], [self.nx_sim - 1, 0], \
                [self.nx_sim - 1, self.ny_sim - 1], \
                [0, self.ny_sim - 1]]
        ll_cnrs = gp2ll_multi(xy_cnrs, self.mlat, self.mlon, self.mrot, \
                self.nx_sim, self.ny_sim, self.hh)

        if not gmt_format:
            return ll_cnrs
        gmt_cnrs = '\n'.join([\
                ' '.join(map(str, cnr)) for cnr in ll_cnrs])
        return ll_cnrs, gmt_cnrs

    def region(self, corners = None):
        """
        Returns simulation region as a tuple (x_min, x_max, y_min, y_max).
        corners: use pre-calculated corners if given
        """
        if corners == None:
            corners = self.corners()
        x_min, y_min = np.min(corners, axis = 0)
        x_max, y_max = np.max(corners, axis = 0)

        return (x_min, x_max, y_min, y_max)

    def tslice_get(self, step, comp = -1, outfile = None):
        """
        Retrieve timeslice data.
        Based on logic in WccFormat/src/ts2xyz.c
        step: timestep to retrieve data for
        """
        # loading x, y, z all the time is not significant
        # compared with logic unless you have repeated code
        y = self.data[step, 0, :, :] * self.cosR \
                - self.data[step, 1, :, :] * self.sinR
        x = self.data[step, 0, :, :] * self.sinR \
                + self.data[step, 1, :, :] * self.cosR
        z = self.data[step, 2, :, :] * -1

        if comp < 0:
            wanted = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        elif comp == 0:
            wanted = x
        elif comp == 1:
            wanted = y
        elif comp == 2:
            wanted = z

        # format as longitude, latitude, value columns
        wanted = np.dstack((self.ll_map, wanted)).reshape((-1, 3))
        if outfile == None:
            return wanted
        else:
            wanted.astype(np.float32).tofile(outfile)

    def pgv(self, mmi = False, pgvout = None, mmiout = None):
        """
        Retrieve PGV map.
        mmi: also calculate MMI
        pgvout: file to store pgv or None to return it
        mmiout: file to store mmi or None to return it
        NOTE: cannot save one of pgvout/mmiout to a file and return another
        """
        # PGV as timeslices reduced to maximum value at each point
        pgv = np.zeros(self.nx * self.ny)
        for ts in xrange(self.t0, self.nt):
            pgv = np.maximum( \
                    np.sqrt(np.sum(np.power(np.dot( \
                    self.data[ts, :, :, :].reshape(3, -1).T, self.rot_matrix \
                    ), 2), axis = 1)), \
                    pgv)

        # modified marcalli intensity formula
        if mmi:
            mmiv = np.where(np.log10(pgv) < 0.53, \
                    3.78 + 1.47 * np.log10(pgv), \
                    2.89 + 3.16 * np.log10(pgv))

        # transform to give longitude, latitude, pgv value
        pgv = np.vstack((self.ll_map.reshape(-1, 2).T, pgv)).T
        if mmi:
            mmiv = np.vstack((self.ll_map.reshape(-1, 2).T, mmiv)).T

        # store / output
        if pgvout != None:
            pgv.astype(np.float32).tofile(pgvout)
        if mmi and mmiout != None:
            mmiv.astype(np.float32).tofile(mmiout)

        if pgvout == None:
            if not mmi:
                return pgv
            else:
                return pgv, mmiv
