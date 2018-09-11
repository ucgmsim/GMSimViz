"""
Read xyts.e3d files.
C structs available in WccFormat/src/structure.h
XYTS file related functions.
XYTS files contain time slices on the X-Y plane (z = 1, top level).
XYTS file processing

@author Viktor Polak
@date 19 January 2017
"""

from math import radians, cos, sin

import numpy as np

from gmsimviz import geo

###
### PROCESSING OF XYTS FILE
###
class XYTSFile:
    """
    Assumptions:
    dip = 0: simulation domain is flat
    t0 = 0: complete timeseries from t = 0
    """

    def __init__(self, xyts_path, meta_only=False):
        """
        Load metadata and optionally prepare gridpoint datum locations.
        xyts_path: path to the xyts.e3d file
        meta_only: don't prepare gridpoint datum locations (slower)
                can't use timeslice (lon, lat, value) capability
        """

        xytf = open(xyts_path, "rb")

        # determine endianness, an x-y timeslice has 1 z value
        nz = np.fromfile(xytf, dtype=">i4", count=7)[-1]
        if nz == 0x00000001:
            endian = ">"
        elif nz == 0x01000000:
            endian = "<"
        else:
            xytf.close()
            raise ValueError("File is not an XY timeslice file: %s" % (xyts_path))
        xytf.seek(0)

        # read header
        self.x0, self.y0, self.z0, self.t0, self.nx, self.ny, self.nz, self.nt = np.fromfile(
            xytf, dtype="%si4" % (endian), count=8
        )
        self.dx, self.dy, self.hh, self.dt, self.mrot, self.mlat, self.mlon = np.fromfile(
            xytf, dtype="%sf4" % (endian), count=7
        )
        xytf.close()

        # determine original sim parameters
        self.dxts = int(round(self.dx / self.hh))
        self.dyts = int(round(self.dy / self.hh))
        self.nx_sim = self.nx * self.dxts
        self.ny_sim = self.ny * self.dyts

        # orientation of components
        self.dip = 0
        self.comps = {
            "X": radians(90 + self.mrot),
            "Y": radians(self.mrot),
            "Z": radians(90 - self.dip),
        }
        # rotation of components so Y is true north
        self.cosR = cos(self.comps["X"])
        self.sinR = sin(self.comps["X"])
        # simulation plane always flat, dip = 0
        self.cosP = 0  # cos(self.comps['Z'])
        self.sinP = 1  # sin(self.comps['Z'])
        # xy dual component rotation matrix
        # must also flip vertical axis
        theta = radians(self.mrot)
        self.rot_matrix = np.array(
            [[cos(theta), -sin(theta), 0], [-sin(theta), -cos(theta), 0], [0, 0, -1]]
        )

        # save speed when only loaded to read metadata section
        if meta_only:
            return

        # memory map for data section
        self.data = np.memmap(
            xyts_path,
            dtype="%sf4" % (endian),
            mode="r",
            offset=60,
            shape=(self.nt, len(self.comps), self.ny, self.nx),
        )

        # create longitude, latitude map for data
        grid_points = (
            np.mgrid[0 : self.nx_sim : self.dxts, 0 : self.ny_sim : self.dyts]
            .reshape(2, -1, order="F")
            .T
        )
        amat = geo.gen_mat(self.mrot, self.mlon, self.mlat)[0]
        self.ll_map = geo.xy2ll(
            geo.gp2xy(grid_points, self.nx_sim, self.ny_sim, self.hh), amat
        ).reshape(self.ny, self.nx, 2)

    def corners(self, gmt_format=False):
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
        gp_cnrs = np.array(
            [
                [0, 0],
                [self.nx_sim - 1, 0],
                [self.nx_sim - 1, self.ny_sim - 1],
                [0, self.ny_sim - 1],
            ]
        )
        amat = geo.gen_mat(self.mrot, self.mlon, self.mlat)[0]
        ll_cnrs = geo.xy2ll(
            geo.gp2xy(gp_cnrs, self.nx_sim, self.ny_sim, self.hh), amat
        ).tolist()
        if not gmt_format:
            return ll_cnrs

        gmt_cnrs = "\n".join([" ".join(map(str, cnr)) for cnr in ll_cnrs])
        return ll_cnrs, gmt_cnrs

    def region(self, corners=None):
        """
        Returns simulation region as a tuple (x_min, x_max, y_min, y_max).
        corners: use pre-calculated corners if given
        """
        if corners is None:
            corners = self.corners()
        x_min, y_min = np.min(corners, axis=0)
        x_max, y_max = np.max(corners, axis=0)

        return (x_min, x_max, y_min, y_max)

    def tslice_get(self, step, comp=-1, outfile=None):
        """
        Retrieve timeslice data.
        Based on logic in WccFormat/src/ts2xyz.c
        step: timestep to retrieve data for
        comp: timestep component -1:sqrt(x^2 + y^2 + z^2), 0:x, 1:y, 2:z
        """
        # loading x, y, z all the time not significant
        y = self.data[step, 0, :, :] * self.cosR - self.data[step, 1, :, :] * self.sinR
        x = self.data[step, 0, :, :] * self.sinR + self.data[step, 1, :, :] * self.cosR
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

    def pgv(self, mmi=False, pgvout=None, mmiout=None):
        """
        Retrieve PGV map.
        mmi: also calculate MMI
        pgvout: file to store pgv or None to return it
        mmiout: file to store mmi or None to return it
        """
        # PGV as timeslices reduced to maximum value at each point
        pgv = np.zeros(self.nx * self.ny)
        for ts in xrange(self.t0, self.nt):
            pgv = np.maximum(
                np.sqrt(
                    np.sum(
                        np.power(
                            np.dot(
                                self.data[ts, :, :, :].reshape(3, -1).T, self.rot_matrix
                            ),
                            2,
                        ),
                        axis=1,
                    )
                ),
                pgv,
            )

        # modified marcalli intensity formula
        if mmi:
            mmiv = np.where(
                np.log10(pgv) < 0.53,
                3.78 + 1.47 * np.log10(pgv),
                2.89 + 3.16 * np.log10(pgv),
            )

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
            elif mmiout == None:
                return pgv, mmiv
