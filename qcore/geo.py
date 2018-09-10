"""
Various tools which may be needed in various processes.
"""

from math import sin, asin, cos, atan, atan2, degrees, radians, sqrt, pi

import numpy as np

R_EARTH = 6378.139


class InputError(Exception):
    pass


def gen_mat(mrot, mlon, mlat):
    """
    Precursor for xy2ll and ll2xy functions.
    mrot: model rotation
    mlon: model centre longitude
    mlat: model centre latitude
    """
    arg = radians(mrot)
    cosA = cos(arg)
    sinA = sin(arg)

    arg = radians(90.0 - mlat)
    cosT = cos(arg)
    sinT = sin(arg)

    arg = radians(mlon)
    cosP = cos(arg)
    sinP = sin(arg)

    amat = np.array(
        [
            [
                cosA * cosT * cosP + sinA * sinP,
                sinA * cosT * cosP - cosA * sinP,
                sinT * cosP,
            ],
            [
                cosA * cosT * sinP - sinA * cosP,
                sinA * cosT * sinP + cosA * cosP,
                sinT * sinP,
            ],
            [-cosA * sinT, -sinA * sinT, cosT],
        ],
        dtype="f",
    )
    ainv = amat.T * 1.0 / np.linalg.det(amat)

    return amat.flatten(), ainv.flatten()


def xy2ll(xy_km, amat):
    """
    Converts km offsets to longitude and latitude.
    xy_km: 2D np array of [X, Y] offsets from origin (km)
    amat: from gen_mat function
    """
    x = xy_km[:, 0] / R_EARTH
    sinB = np.sin(x)
    y = xy_km[:, 1] / R_EARTH
    sinG = np.sin(y)
    z = np.sqrt(1.0 + sinB * sinB * sinG * sinG)
    xp = sinG * np.cos(x) * z
    yp = sinB * np.cos(y) * z
    zp = np.sqrt(1.0 - xp * xp - yp * yp)

    xg = xp * amat[0] + yp * amat[1] + zp * amat[2]
    yg = xp * amat[3] + yp * amat[4] + zp * amat[5]
    zg = xp * amat[6] + yp * amat[7] + zp * amat[8]

    lat = np.where(
        zg == 0.0, 0.0, 90.0 - np.degrees(np.arctan(np.sqrt(xg * xg + yg * yg) / zg))
    )
    lat[np.where(zg < 0.0)] -= 180.0

    lon = np.where(xg == 0.0, 0.0, np.degrees(np.arctan(yg / xg)))
    lon[np.where(xg < 0.0)] -= 180.0
    lon[np.where(lon < -180.0)] += 360.0

    return np.column_stack((lon, lat))


def ll2xy(ll, ainv):
    """
    Converts longitude and latitude to km offsets.
    ll: 2D np array of [lon, lat]
    ainv: from gen_mat function
    """
    lon = np.radians(ll[:, 0])
    lat = np.radians(90.0 - ll[:, 1])

    xg = np.sin(lat) * np.cos(lon)
    yg = np.sin(lat) * np.sin(lon)
    zg = np.cos(lat)

    xp = xg * ainv[0] + yg * ainv[1] + zg * ainv[2]
    yp = xg * ainv[3] + yg * ainv[4] + zg * ainv[5]
    zp = xg * ainv[6] + yg * ainv[7] + zg * ainv[8]

    # X km offsets from centre origin, Y km offsets from centre origin
    return np.column_stack(
        (
            R_EARTH * np.arcsin(yp / np.sqrt(1.0 - xp * xp)),
            R_EARTH * np.arcsin(xp / np.sqrt(1.0 - yp * yp)),
        )
    )


def xy2gp(xy, nx, ny, hh):
    """
    Converts km offsets to grid points.
    xy: 2D np array of [X, Y] offsets from origin (km)
    nx: number of X grid positions
    ny: number of Y grid positions
    hh: grid spacing
    """
    gp = np.copy(xy)

    # distance from corner
    gp[:, 0] += (nx - 1) * hh * 0.5
    gp[:, 1] += (ny - 1) * hh * 0.5

    # gridpoint from top corner
    gp /= hh

    return np.round(gp).astype(np.int32, copy=False)


def gp2xy(gp, nx, ny, hh):
    """
    Converts grid points to km offsets.
    xy: 2D np array of [X, Y] gridpoints
    nx: number of X grid positions
    ny: number of Y grid positions
    hh: grid spacing
    """
    xy = gp.astype(np.float32) * hh

    # shift for centre origin
    xy[:, 0] -= (nx - 1) * hh * 0.5
    xy[:, 1] -= (ny - 1) * hh * 0.5

    return xy


def ll_shift(lat, lon, distance, bearing):
    """
    Shift lat/long by distance at bearing.
    """
    # formula is for radian values
    lat, lon, bearing = map(radians, [lat, lon, bearing])

    shift = distance / R_EARTH
    lat2 = asin(sin(lat) * cos(shift) + cos(lat) * sin(shift) * cos(bearing))
    lon2 = lon + atan2(
        sin(bearing) * sin(shift) * cos(lat), cos(shift) - sin(lat) * sin(lat2)
    )

    return degrees(lat2), degrees(lon2)


def ll_mid(lon1, lat1, lon2, lat2):
    """
    Return midpoint between a pair of lat, long points.
    """
    # functions based on radians
    lon1, lat1, lat2, dlon = map(radians, [lon1, lat1, lat2, (lon2 - lon1)])

    Bx = cos(lat2) * cos(dlon)
    By = cos(lat2) * sin(dlon)

    lat3 = atan2(sin(lat1) + sin(lat2), sqrt((cos(lat1) + Bx) ** 2 + By ** 2))
    lon3 = lon1 + atan2(By, cos(lat1) + Bx)

    return degrees(lon3), degrees(lat3)


def ll_dist(lon1, lat1, lon2, lat2):
    """
    Return distance between a pair of lat, long points.
    """
    # functions based on radians
    lat1, lat2, dlon, dlat = map(radians, [lat1, lat2, (lon2 - lon1), (lat2 - lat1)])

    a = sin(dlat / 2.0) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2.0) ** 2
    return R_EARTH * 2.0 * atan2(sqrt(a), sqrt(1 - a))


def ll_bearing(lon1, lat1, lon2, lat2, midpoint=False):
    """
    Initial bearing when traveling from 1 -> 2.
    Direction facing from point 1 when looking at point 2.
    """
    if midpoint:
        lon1, lat1 = ll_mid(lon1, lat1, lon2, lat2)
    lat1, lat2, lon_diff = map(radians, [lat1, lat2, (lon2 - lon1)])
    return (
        degrees(
            atan2(
                cos(lat2) * sin(lon_diff),
                cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon_diff),
            )
        )
        % 360
    )


def angle_diff(b1, b2):
    """
    Return smallest difference (clockwise, -180 -> 180) from b1 to b2.
    """
    r = (b2 - b1) % 360
    if r > 180:
        return r - 360
    return r


def avg_wbearing(angles):
    """
    Return average angle given angles and weightings.
    NB: angles are clockwise from North, not anti-clockwise from East.
    angles: 2d list of (angle, weight)
    """
    x = 0
    y = 0
    for a in angles:
        x += a[1] * sin(radians(a[0]))
        y += a[1] * cos(radians(a[0]))
    q_diff = 0
    if y < 0:
        q_diff = pi
    elif x < 0:
        q_diff = 2 * pi
    return degrees(atan(x / y) + q_diff)


def path_from_corners(
    corners=None, output="sim.modelpath_hr", min_edge_points=100, close=True
):
    """
    corners: python list (4 by 2) containing (lon, lat) in order
        otherwise take from velocity model
    output: where to store path of (lon, lat) values
    min_edge_points: at least this many points wanted along edges
    """
    # input data using velocity model
    if corners is None:
        # don't fail importing if not needed
        from params import vel_mod_params

        # load model_params
        with open(vel_mod_params, "r") as mp:
            lines = mp.readlines()
        # find corners using tags
        for tag in ["c1=", "c2=", "c3=", "c4="]:
            for line in lines:
                if tag in line:
                    corners.append(map(float, line.split()[1:3]))

    # close the box by going back to origin
    if close:
        corners.append(corners[0])
    # until each side has at least wanted number of points
    while len(corners) < 4 * min_edge_points:
        # work backwards, insertions don't change later indexes
        for i in xrange(len(corners) - 1, 0, -1):
            val = ll_mid(
                corners[i][0], corners[i][1], corners[i - 1][0], corners[i - 1][1]
            )
            corners.insert(i, val)

    # write points the make the path
    if output != None:
        with open(output, "w") as mp:
            for point in corners:
                mp.write("%s %s\n" % (point[0], point[1]))
    else:
        return corners


def wgs_nztm2000x(points):
    """
    Coordinates Transform: between WGS84 lon, lat and NZ Tranverse Mercator 2000
    author: X. Bellagamba, Matlab -> Python by Viktor Polak

    points: series of decimal longitude, latitude or nztm2k x, y
    """
    if type(points).__name__ == "list":
        points = np.array(points)
    points = np.atleast_2d(points)

    if np.max(np.abs(points)) > 180:
        wgs_out = True
    else:
        wgs_out = False
        points = np.radians(points)

    # NZTM2000 definitions (PROJ4 naming)
    # origin
    lon_0 = radians(173.0)
    lat_0 = 0.0
    # false northing, easting (metres)
    y_0 = 10000000
    x_0 = 1600000
    # central meridian scaling factor
    k_0 = 0.9996
    # flattening of the ellipsoid (%)
    f = 1 / 298.257222101
    # semimajor, semiminor ellipsoid radius
    a = 6378137.0
    b = a * (1 - f)
    # eccentricity of the ellipsoid squared
    es = (f + f) - (f * f)

    # conversion specific
    A_0 = 1 - es / 4.0 - 3 * es ** 2 / 64.0 - 5 * es ** 3 / 256.0
    A_2 = 3 / 8.0 * (es + es ** 2 / 4.0 + 15 * es ** 3 / 128.0)
    A_4 = 15 / 256.0 * (es ** 2 + 3 * es ** 3 / 4.0)
    A_6 = 35 * es ** 3 / 3072.0
    m_0 = a * (
        A_0 * lat_0 - A_2 * sin(2 * lat_0) + A_4 * sin(4 * lat_0) - A_6 * sin(6 * lat_0)
    )
    n = (a - b) / (a + b)
    G = a * (1 - n) * (1 - n ** 2) * radians(1 + 9 * n ** 2 / 4.0 + 225 * n ** 4 / 64.0)

    ###
    ### process all points
    ###
    if wgs_out:
        N_prime = points[:, 1] - y_0
        E_prime = points[:, 0] - x_0
        m_prime = m_0 + N_prime / k_0
        sigma = np.radians(m_prime) * 1 / G
        phi_prime = (
            sigma
            + (3 * n / 2.0 - 27 * n ** 3 / 32.0) * np.sin(2 * sigma)
            + (21 * n ** 2 / 16.0 - 55 * n ** 4 / 32.0) * np.sin(4 * sigma)
            + (151 * n ** 3 / 96.0) * np.sin(6 * sigma)
            + (1097 * n ** 4 / 512.0) * np.sin(8 * sigma)
        )
        y_factors = phi_prime
    else:
        y_factors = points[:, 1]
        m = a * (
            A_0 * y_factors
            - A_2 * np.sin(2 * y_factors)
            + A_4 * np.sin(4 * y_factors)
            - A_6 * np.sin(6 * y_factors)
        )
        omega = points[:, 0] - lon_0

    # ellipsoid radius in the prime vertical
    nu = a / np.sqrt(1 - es * np.sin(y_factors) ** 2)
    # projection parameters
    rho = (a * (1 - es)) / ((1 - es * np.sin(y_factors) ** 2) ** 1.5)
    psi = nu / rho
    t = np.tan(y_factors)
    if wgs_out:
        rs = rho * nu * k_0 ** 2
        x = E_prime / (k_0 * nu)

        # terms for the north coordinates
        T1_N = t / (k_0 * rho) * E_prime * x / 2.0
        T2_N = (
            t
            / (k_0 * rho)
            * E_prime
            * x ** 3
            / 24.0
            * (-4 * psi ** 2 + 9 * psi * (1 - t ** 2) + 12 * t ** 2)
        )
        T3_N = (
            t
            / (k_0 * rho)
            * E_prime
            * x ** 5
            / 720.0
            * (
                8 * psi ** 4 * (11 - 24 * t ** 2)
                - 12 * psi ** 3 * (21 - 7 * t ** 2)
                + 15 * psi ** 2 * (15 - 98 * t ** 2 + 15 * t ** 4)
                + 180 * psi * (5 * t ** 2 - 3 * t ** 4)
                + 360 * t ** 4
            )
        )
        T4_N = (
            t
            / (k_0 * rho)
            * E_prime
            * x ** 7
            / 40320.0
            * (1385 + 3633 * t ** 2 + 4095 * t ** 4 + 1575 * t ** 6)
        )
        # north coordinates
        y = y_factors - T1_N + T2_N - T3_N + T4_N

        # terms for the east coordinates
        T1_E = x * 1 / np.cos(y_factors)
        T2_E = x ** 3 * 1 / np.cos(y_factors) / 6 * (psi + 2 * t ** 2)
        T3_E = (
            x ** 5
            * 1
            / np.cos(y_factors)
            / 120
            * (
                -4 * psi ** 3 * (1 - 6 * t ** 2)
                + psi ** 2 * (9 - 68 * t ** 2)
                + 72 * psi * t ** 2
                + 24 * t ** 4
            )
        )
        T4_E = (
            x ** 7
            * 1
            / np.cos(y_factors)
            / 5040
            * (61 + 662 * t ** 2 + 1320 * t ** 4 + 720 * t ** 6)
        )
        # east coordinates
        x = lon_0 + T1_E - T2_E + T3_E - T4_E

        return np.dstack((np.degrees(x), np.degrees(y)))[0]

    # terms for the north coordinates
    T1_N = (omega ** 2 / 2.0) * nu * np.sin(y_factors) * np.cos(y_factors)
    T2_N = (
        (omega ** 4 / 24.0)
        * nu
        * np.sin(y_factors)
        * np.cos(y_factors) ** 3
        * (4 * psi ** 2 + psi - t ** 2)
    )
    T3_N = (
        (omega ** 6 / 720.0)
        * nu
        * np.sin(y_factors)
        * np.cos(y_factors) ** 5
        * (
            8 * psi ** 4 * (11 - 24 * t ** 2)
            - 28 * psi ** 3 * (1 - 6 * t ** 2)
            + psi ** 2 * (1 - 32 * t ** 2)
            - psi * (2 * t ** 2)
            + t ** 4
        )
    )
    T4_N = (
        (omega ** 8 / 40320.0)
        * nu
        * np.sin(y_factors)
        * np.cos(y_factors) ** 7
        * (1385 - 3111 * t ** 2 + 543 * t ** 4 - t ** 6)
    )
    # north coordinates
    lat = y_0 + k_0 * (m - m_0 + T1_N + T2_N + T3_N + T4_N)

    # terms for the east coordinates
    T1_E = (omega ** 2 / 6.0) * np.cos(y_factors) ** 2 * (psi - t ** 2)
    T2_E = (
        (omega ** 4 / 120.0)
        * np.cos(y_factors) ** 4
        * (
            4 * psi ** 3 * (1 - 6 * t ** 2)
            + psi ** 2 * (1 + 8 * t ** 2)
            - psi * 2 * t ** 2
            + t ** 4
        )
    )
    T3_E = (
        (omega ** 6 / 5040.0)
        * np.cos(y_factors) ** 6
        * (61 - 479 * t ** 2 + 179 * t ** 4 - t ** 6)
    )
    # east coordinates
    lon = x_0 + k_0 * nu * omega * np.cos(y_factors) * (1 + T1_E + T2_E + T3_E)

    return np.dstack((lon, lat))[0]
