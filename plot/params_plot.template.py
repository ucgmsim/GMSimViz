"""
Parameters to alter plotting.
Template should be copied to sim_dir as params_plot.py
"""
import os

event_title = 'Automatically Genenerated Event'
fault_model = 'Automatically generated source model'
# <HH> is converted automatically
vel_model = 'SIVM v1.65 h=<HH>km'
# pick a predefined region or model_params sim domain used ('')
# 'CANTERBURY', 'WIDERCANT', 'MIDNZ', 'SOUTHISLAND'
region = None

# topography at different resolutions
topo_file_low = '/nesi/projects/nesi00213/PlottingData/Topo/nztopo.grd'
topo_file_high = '/nesi/projects/nesi00213/PlottingData/Topo/srtm_all_filt_nz.grd'

# timeslice plotting
class TS:
    # which component to plot
    # -1 : geometric mean, 0 : X, 1 : Y, 2 : Z
    component = -1
    # width of map area
    width = '6i'
    major_tick = '1d'
    minor_tick = '30m'
    # timeslices should be relatively low DPI for video
    dpi = 96
    # larger overlay grd limit improves speed, too large = bad results
    # default is about 0.2 * cpt increment
    convergence_limit = None
    # colour scale for overlay
    cpt = 'hot'
    # override the background and foreground colours
    cpt_bg = None
    cpt_fg = None
    # invert colour range
    cpt_inv = True
    cpt_min = 0
    # max and increment override
    cpt_max = None
    cpt_inc = None
    # values below lowcut or above highcut removed
    # if lowcut > highcut then the range inbetween is removed
    # None - disable, 'auto' - automatic (about 2% of cpt max)
    # TODO: 'auto' only implemented for lowcut
    lowcut = 'auto'
    highcut = None
    cpt_legend = 'ground motion (cm/s)'
    # crop overlay to land
    land_crop = False
    # overlay grid resolution in longitude, latitude units
    # append 'k' for killometres, 'e' for metres
    # accurate at middle of grid, equal-longitude/latitude spacing
    # None is ideal smooth resolution for data (may be slower than wanted)
    # examples: '1k' 1 kilometre, '400e' 400 metres, '0.002' 0.002 lat or lon
    grd_dx = None
    grd_dy = None
    # seismogram overlay data, thickness, colour and format
    # only plotted if the file exists
    seis_data = 'gmt-seismo.xy'
    seis_line = '0.8p'
    seis_colour = 'red'
    # 'time' if data move, 'inc' if data extends
    seis_fmt = 'time'

# PGV plotting
class PGV(TS):
    dpi = 300
    title = event_title + ' PGV'
    # cpt for PGV should be kept static for the most part
    cpt = 'hot'
    cpt_min = 0
    # don't exclude values close to 0
    lowcut = None
    cpt_legend = 'peak ground velocity (cm/s)'
    # crop overlay to land
    land_crop = False

# MMI plotting
class MMI:
    dpi = 300
    title = event_title + ' MMI'
    cpt_legend = 'modified mercalli intensity'
    # MMI deals with +- 1 values, needs smaller convergence limit
    convergence_limit = 0.1
    # crop MMI to land
    land_crop = False

# station plotting
class STATION:
    dpi = 300
    width = '6i'
    topo_file = topo_file_high
    out_dir = os.path.abspath('PNG_stations')
    # override velocity model region default = None
    # eg: region = (x_min, x_max, y_min, y_max)
    region = None
    # major, minor tick increment on map edges. None for automatic
    tick_major = None
    tick_minor = None
    # 'major' only include major sites, 'all' include all sites,
    # 'auto' to choose between major and all, None for [] (no sites listed)
    # or specify list of sites manualy ['Kaikoura', 'Wellington']...
    sites = 'auto'

# observed / simulated seismogram plotting
class SEISMO:
    title = 'Observed Ground Motions'
    wd = os.path.abspath('GMT_SEIS_WD')
    # output filename excluding file extension
    name = 'ObservedMap'
    dpi = 300
    width = '7i'
    # override velocity model region default = None
    # eg: region = (x_min, x_max, y_min, y_max)
    region = None
    # list of stations that should always be plotted
    wanted_stats = []
    # minimum distance stations should be appart (km)
    min_dist = None
    # GMT seismo files
    # set obs_src or sim_src to None to only plot one
    obs_ts = 'obsVelBB/%s.090'
    obs_src = 'gmt-seismo_obs.xy'
    obs_colour = 'black'
    sim_ts = 'simVelBB/%s.090'
    sim_src = 'gmt-seismo_sim.xy'
    sim_colour = 'red'
    seis_width = 0.3
    # timestep cutoff: with dt = 0.005, 200 per second
    max_ts = 20000
    # timeseries x-azimuth, x length and max y length in degrees
    ts_xaz = 90
    ts_xlen = None
    ts_ymax = None

class SRF:
    # default is the filename without extention
    title = None
    # PNG output, 'srfdir' for same location as SRF
    out_dir = 'srfdir'
    gmt_temp = os.path.abspath('GMT_WD_SRF2D')
    # use average rake (True) or midpoint rake (False)
    rake_average = False
    # length of the longest rake arrow (based on slip)
    rake_length = 0.4
    # place rake arrows every X subfaults, less than 1 for automatic
    rake_decimation = 0

