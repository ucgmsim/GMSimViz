"""

GENERAL SRF NOTES:
SHYP: hypocentre position along strike from centre
DHYP: hypocentre position along dip from top
Common SRF functions.

SRF format:
https://scec.usc.edu/scecpedia/Standard_Rupture_Format
"""

from math import ceil, cos, floor, radians, sqrt
import os
from subprocess import Popen, PIPE

import numpy as np

# assumption that all srf files contain 6 values per line
VPL = 6.0


def read_header(sf, idx=False):
    """
    Parse header information.
    sf: open srf file at position 0 or filename
    """
    # detect if file already open
    try:
        sf.name
        close_me = False
    except AttributeError:
        sf = open(sf, "r")
        close_me = True

    version = float(sf.readline())
    line2 = sf.readline().split()
    # this function requires the optional PLANE header
    assert line2[0] == "PLANE"
    nseg = int(line2[1])
    planes = []
    for _ in range(nseg):
        # works for version 1.0 and 2.0
        elon, elat, nstk, ndip, ln, wid = sf.readline().split()
        stk, dip, dtop, shyp, dhyp = map(float, sf.readline().split())
        # return as dictionary
        if idx:
            planes.append(
                {
                    "centre": [float(elon), float(elat)],
                    "nstrike": int(nstk),
                    "ndip": int(ndip),
                    "length": float(ln),
                    "width": float(wid),
                    "strike": stk,
                    "dip": dip,
                    "shyp": shyp,
                    "dhyp": dhyp,
                    "dtop": dtop,
                }
            )
        else:
            # TODO: deprecated, causes confusion
            planes.append(
                (
                    float(elon),
                    float(elat),
                    int(nstk),
                    int(ndip),
                    float(ln),
                    float(wid),
                    stk,
                    dip,
                    dtop,
                    shyp,
                    dhyp,
                )
            )
    if close_me:
        sf.close()
    return planes


def skip_points(sf, np):
    """
    Skips wanted number of points entries in SRF.
    sf: open srf file at the start of a point
    np: number of points to read past
    """
    for _ in range(np):
        # header 1 not important
        sf.readline()
        # header 2 contains number of values which determines lines
        values = sum(map(int, sf.readline().split()[2::2]))
        for _ in range(int(ceil(values / VPL))):
            sf.readline()


def get_lonlat(sf, value=None, depth=False):
    """
    Returns only the longitude, latitude of a point.
    sf: open file at start of point
    value: also retrieve value
    depth: return value in 3d space (return lon, lat, depth, value)
    end: sf at start of next point
    """
    # header 1 contains:
    # LON, LAT, DEP, STK, DIP, AREA, TINIT, DT, VS (v2.0), DEN (v2.0)
    h1 = sf.readline().split()
    # header 2 contains:
    # RAKE, SLIP1, NT1, SLIP2, NT2, SLIP3, NT3
    h2 = sf.readline().split()

    # always returning lon, lat
    lon, lat, depth_v = map(float, h1[:3])
    if value == "slip":
        value = sqrt(float(h2[1]) ** 2 + float(h2[3]) ** 2 + float(h2[5]) ** 2)
    elif value == "tinit":
        value = float(h1[6])
    elif value == "trise":
        value = max(map(int, h2[2::2])) * float(h1[7])
    elif value == "ttotal":
        value = float(h1[6]) + max(map(int, h2[2::2])) * float(h1[7])
    elif value == "depth":
        value = float(h1[2])
    elif value == "rake":
        # rake in the file is actually the u1 axis
        # it just happens that this is usually adjusted to equal rake
        # but we shouldn't rely on this so check if there is a u2 component
        if float(h2[3]) == 0:
            # usually an int but allow float
            value = float(h2[0])
        else:
            # unadjusted u1 - full formula, expected to never be required
            u1 = radians(float(h2[0]))
            strike_slip = cos(u1) * float(h2[1]) - sin(u1) * float(h2[3])
            dip_slip = sin(u1) * float(h2[1]) + cos(u1) * float(h2[3])
            try:
                value = degrees(atan(strike_slip / dip_slip))
                if strike_slip < 0:
                    value += 180
            except ZeroDivisionError:
                if dip_slip == 0:
                    value = float(h2[0])
                else:
                    value = 90
            value = value % 360
    elif value == "dt":
        value = float(h1[7])

    # skip rest of point data
    # or return the time series
    values = sum(map(int, h2[2::2]))
    if type(value).__name__ != "str" or (
        value[:8] != "sliprate" and value[:6] != "slipts"
    ):
        for _ in range(int(ceil(values / VPL))):
            sf.readline()
    else:
        cumulative = False
        if value[:6] == "slipts":
            cumulative = True
        # store rest of point data
        srate = []
        for _ in range(int(ceil(values / VPL))):
            srate.extend(map(float, sf.readline().split()))
        # sliprate-dt-tend
        dt, t = value.split("-")[1:3]
        tinit = float(h1[6])
        srfdt = float(h1[7])
        # time series over wanted time
        # python thinks int(9.6/0.1) is 95, rounding is a must
        value = np.empty(int(round(float(t) / float(dt))))
        value.fill(np.nan)
        # fill with values during rupture period at this subfault
        i = 0
        for r in range(int(h2[2])):
            # time index as decimated
            i = int(floor((tinit + r * srfdt) / float(dt)))
            if cumulative:
                value[i] = sum(srate[: r + 1]) * srfdt
            elif np.isnan(value[i]):
                # first value at this point
                value[i] = srate[r]
                # repeating factor for averaging
                x = 1.0
            else:
                x += 1.0
                # average of all values up to now
                value[i] += (srate[r] - value[i]) / x
        if cumulative:
            # copy total sum to end
            value[i:] = value[i]

    if type(value).__name__ == "NoneType":
        if depth:
            return lon, lat, depth_v
        return lon, lat
    if depth:
        return lon, lat, depth_v, value
    return lon, lat, value


def get_bounds(srf, seg=-1, depth=False):
    """
    Return corners of segments.
    srf: srf source
    nseg: which segment (-1 for all)
    depth: also include depth if True
    """
    bounds = []
    with open(srf, "r") as sf:
        # metadata
        planes = read_header(sf)
        points = int(sf.readline().split()[1])

        # third value to retrieve after longitude, latitude
        if depth:
            value = "depth"
        else:
            value = None
        # each plane has a separate set of corners
        for n, plane in enumerate(planes):
            plane_bounds = []
            nstk, ndip = plane[2:4]
            # set of points starts at corner
            plane_bounds.append(get_lonlat(sf, value=value))
            # travel along strike, read last value
            skip_points(sf, nstk - 2)
            plane_bounds.append(get_lonlat(sf, value=value))
            # go to start of strike at bottom of dip
            skip_points(sf, (ndip - 2) * nstk)
            plane_bounds.append(get_lonlat(sf, value=value))
            # travel along strike at bottom of dip
            skip_points(sf, nstk - 2)
            plane_bounds.insert(2, get_lonlat(sf, value=value))
            # store plane bounds or return if only 1 wanted
            if n == seg:
                return plane_bounds
            bounds.append(plane_bounds)
    return bounds


def get_hypo(srf, lonlat=True, depth=False):
    """
    Return hypocentre.
    srf: srf source file path
    lonlat: in tearms of longitude and latitude (True), raw km offsets (False)
    depth: return longitude, latitude and depth. False to only return 2 values
    """
    # complexity increased by:
    # segments may be shared, shyp relative to centre of shared segments
    # have to be able to travel horizontally between shared segments
    with open(srf, "r") as sf:
        planes = read_header(sf)
        # keep only main segment set (containing hypocentre)
        if len(planes) > 1:
            for i in range(1, len(planes)):
                # negative dhyp if same segment set
                if planes[i][10] >= 0:
                    del planes[i:]
                    break
        # check for point source
        elif sum(planes[0][2:4]) == 2:
            # returning offsets doesn't make sense
            # should have already checked for point source before calling this
            assert lonlat
            points = int(sf.readline().split()[1])
            assert points == 1
            hlon, hlat, depth_km = get_lonlat(sf, value="depth")
            if not depth:
                return hlon, hlat
            return hlon, hlat, depth_km

        # dip will be constant along shared segments
        ndip = planes[0][3]
        wid = planes[0][5]
        shyp, dhyp = planes[0][9:11]
        hyp_dip = int(round(dhyp / (wid / (float(ndip) - 1))))
        try:
            assert 0 <= hyp_dip < ndip
        except AssertionError:
            # describe scenario for debugging
            print("Hypocentre determined to be outside fault width.")
            print(
                "Width is %skm with %skm subfault spacing."
                % (wid, (wid / (float(ndip) - 1)))
            )
            print(
                "Hypocentre dip subfault number %d, total subfaults: %d."
                % (hyp_dip, ndip)
            )
            raise

        # XXX: there are gaps between segments, ignored
        # total segment set length (along strike)
        ln_tot = sum(planes[p][4] for p in range(len(planes)))
        # distance from group start edge to hypocentre
        ln_shyp = ln_tot / 2.0 + planes[0][9]
        ln_shyp_rel = ln_shyp
        # calculate shyp within correct sub segment
        ln_segs = 0
        for p in range(len(planes)):
            ln_segs += planes[p][4]
            if ln_segs >= ln_shyp:
                break
            ln_shyp_rel -= planes[p][4]
        # give flat projection location
        if not lonlat:
            if depth:
                return p, ln_shyp_rel, dhyp
            return p, ln_shyp_rel
        # determine strike in correct sub segment
        nstk = planes[p][2]
        ln = planes[p][4]
        hyp_stk = int(round(ln_shyp_rel / (ln / (float(nstk) - 1))))
        try:
            assert 0 <= hyp_stk < nstk
        except AssertionError:
            print("Hypocentre determined to be outside fault length.")
            print(
                "Length of sub-segment %d of %d is %skm with %skm subfault spacing."
                % (p, len(planes), ln, ln / float(nstk) - 1)
            )
            print(
                "Hypocentre strike subfault number %d, total subfaults: %d."
                % (hyp_stk, nstk)
            )
            raise

        # retrieve value
        points = int(sf.readline().split()[1])
        for skip in range(p):
            nstk, ndip = planes[skip][2:4]
            skip_points(sf, nstk * ndip)
        nstk, ndip = planes[p][2:4]
        ln, wid = planes[p][4:6]
        skip_points(sf, hyp_dip * nstk + hyp_stk)
        hlon, hlat, depth_km = get_lonlat(sf, value="depth")

        if not depth:
            return hlon, hlat
        return hlon, hlat, depth_km


def srf2llv_py(srf, value="slip", seg=-1, lonlat=True, depth=False, flip_rake=False):
    """
    Return list of lon, lat, value for subfaults in each plane.
    # Reading all at once is faster than reading each separate because:
    #    reading through text is slowest part, seeking not possible.
    # speed for a large file (7 seg, 216k subfaults, slip):
    # All in python version: 3 seconds, result separated by planes
    # All in srf2xyz code: 6.5 seconds, result not separated by planes
    # Each plane with srf2xyz code: 6.5 seconds * 7 = 40 seconds
    # Should replace srf2llv.
    srf: srf file path
    nseg: which segment (-1 for all)
    lonlat: return lon lat (True) or x y (False)
    depth: return depth as well (lonlat must be true) (lon,lat,depth,value)
    flip_rake: angles given as -180 -> 180 instead of 0 -> 360 degrees
    """

    # if we want a whole series of values for each subfault
    multi = value[:8] == "sliprate" or value[:6] == "slipts"

    with open(srf, "r") as sf:
        # metadata
        planes = read_header(sf, idx=True)
        points = int(sf.readline().split()[1])

        # containers for planes
        values = []
        if multi:
            series = []

        # each plane has a separate set of subfaults
        for n, plane in enumerate(planes):
            # skip unwanted plane
            if seg >= 0 and seg != n:
                skip_points(sf, plane["nstrike"] * plane["ndip"])
                continue

            # numpy containers for plane
            plane_values = np.zeros(
                (plane["nstrike"] * plane["ndip"], 3 - multi + depth)
            )
            if multi:
                plane_series = [None] * (plane["nstrike"] * plane["ndip"])

            if not lonlat:
                # calculate x, y offsets
                # unlike srf2xyz offsets are relative to segment
                dx = plane["length"] / plane["nstrike"]
                dy = plane["width"] / plane["ndip"]
                # fill with x, y coord grid
                plane_values[:, :2] = (
                    np.mgrid[
                        0.5 * dx : plane["length"] : dx, 0.5 * dy : plane["width"] : dy
                    ]
                    .reshape(2, -1, order="F")
                    .T
                )
                # last item - values from SRF
                if multi:
                    for i in range(plane["nstrike"] * plane["ndip"]):
                        plane_series[i] = get_lonlat(sf, value=value)[-1]
                else:
                    for i in range(plane["nstrike"] * plane["ndip"]):
                        plane_values[i, 2] = get_lonlat(sf, value=value)[-1]

            else:
                if not multi:
                    for i in range(plane["nstrike"] * plane["ndip"]):
                        plane_values[i] = get_lonlat(sf, value=value, depth=depth)
                else:
                    for i in range(plane["nstrike"] * plane["ndip"]):
                        if depth:
                            plane_values[i][0], plane_values[i][1], plane_values[i][
                                2
                            ], plane_series[i] = get_lonlat(
                                sf, value=value, depth=depth
                            )
                        else:
                            plane_values[i][0], plane_values[i][1], plane_series[
                                i
                            ] = get_lonlat(sf, value=value)

            # adjust angles to -180 -> 180 degrees
            if flip_rake and value == "rake":
                np.where(
                    plane_values[:, 2 + depth] > 180,
                    plane_values[:, 2 + depth] - 360,
                    plane_values[:, 2 + depth],
                )

            # add plane specific dimention numpy arrays to return list
            values.append(plane_values)
            if multi:
                series.append(np.array(plane_series))

            # do not read rest of srf file (slow)
            if n == seg:
                break

    if not multi:
        return values
    return values, series


def srf_dxy(srf):
    """
    Retrieve SRF dx and dy.
    Assumes all planes have same dx, dy.
    srf: SRF file path to read from
    """
    with open(srf, "r") as sf:
        # version
        sf.readline()
        # planes definition
        sf.readline()
        # first plane
        elon, elat, nstk, ndip, length, width = sf.readline().split()
    return (
        float("%.2f" % (float(length) / int(nstk))),
        float("%.2f" % (float(width) / int(ndip))),
    )


def srf_dt(srf):
    """
    Retrieve SRF dt value.
    timestep in velocity function (sec)
    """
    with open(srf, "r") as sf:
        # skip metadata
        read_header(sf)
        # skip number of points
        sf.readline()
        # dt from first point
        return float(sf.readline().split()[7])
