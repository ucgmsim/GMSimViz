# XYTS.e3d binary format
This file is produced by EMOD3D and contains a timeseries of ground motions on the XY plane. Unlike the LF seis files, this contains data at all grid points and may have a decimated resolution specified when running EMOD3D through the `e3d.par` file with the parameters `dxts` and `dyts`.

Numbers are 4 bytes in length and may be little or big endian.

File size can be derived knowing the format and the shape of time-series (all necessary values are at the beginning of the file).

The gridpoints are based on a model which is an area with equidistant gridpoints in the X, Y, and Z directions. It is centred on a position (longitude, latitude) and may be rotated.

- simulation metadata
    - **INTEGERS**
        - number of first x gridpoint
        - number of first y gridpoint
        - number of first z gridpoint
        - number of first timestep
        - number of x gridpoints
        - number of y gridpoints
        - number of z gridpoins (always 1 by definition of X-Y file)
        - number of timesteps
    - **FLOATS**
        - x spacing between given gridpoints (km)
        - y spacing between given gridpoints (km)
        - original (pre-decimated) grid spacing between gridpoints used in simulation (km)
        - timestep in timeseries (s)
        - model rotation of gridpoints (degrees)
        - model centre latitude (degrees)
        - model centre longitude (degrees)
- timeseries
    - float array of velocities in the dimentions of timesteps, components (x, y, z), y grid positions, x grid positions.

