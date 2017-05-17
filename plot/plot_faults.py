#!/usr/bin/env python2
"""
Plot faults that have been studied by QuakeCoRE
"""

import os

import qcore_path
import gmt
print gmt.sites_major
exit()

title = 'QuakeCoRE Simulations'

faults_folder = '/home/nesi00213/PlottingData/Earthquakes'
faults_historic = ['13Juneb2011FaultPlane.xy', \
         \
        '22Feb2011FaultPlane.xy', \
        '23Deca2011FaultPlane.xy', '23Decb2011FaultPlane.xy', \
        '20161114 Kaikoura Ian02_s103245.txt', \
        'Cook Strait 012017.txt', 'DarfieldFaultPlane.xy', \
        'Wilberforce 122016.txt', '2014Jan20_Eketahuna_m6p3.txt']
faults_future = ['Alpine Fault m7.90-411.0x17.3_s1129570.txt', \
        'Porters Pass M7.20_82.0by14.0_s83747.txt']
srf_files = ['20161114 Kaikoura Ian02_s103245.srf', \
        '20110222 m6p2.srf', '20100904 m7p1.srf',
        '20110613b m6.00-13.0x9.0_s19126.srf', \
        'Alpine Fault m7.90-411.0x17.3_s1129570.srf', \
        'Porters Pass M7.20_82.0by14.0_s83747.srf', \
        'Wilberforce standard_m5.60-5.1x5.1_s103245.srf', \
        'Cook Strait 012017.srf', \
        '2014Jan20_Eketahuna_m6p3.srf']
faults = '/nesi/projects/nesi00213/PlottingData/Paths/faults/FAULTS_20161219.ll'
# convert to absolute paths
rel2abs = lambda fault_file : os.path.join(faults_folder, fault_file)
faults_historic = map(rel2abs, faults_historic)
faults_future = map(rel2abs, faults_future)
srf_files = map(rel2abs, srf_files)

labels = [ \
        [167.4, -44.7, 'RB', 'Alpine Fault'], \
        [173.6, -42.6, 'LB', 'Kaikoura'], \
        [174.3, -41.7, 'LB', 'Cook Strait'], \
        [171.4, -43, 'LB', 'Wilberforce'], \
        [175.95, -40.6, 'LB', 'Eketahuna'], \
        [171.5, -43.4, 'RB', 'Porters Pass'], \
        [172, -43.7, 'RB', 'Darfield'], \
        [172.9, -43.5, 'LB', 'Christchurch'] \
        ]

dpi = 600
script_dir = os.path.abspath(os.path.dirname(__file__))
gmt_temp = os.path.abspath('GMT_WD_FAULTS')
if not os.path.exists(gmt_temp):
    os.makedirs(gmt_temp)
srf_data = []
for i, f in enumerate(srf_files):
    srf_data.append(gmt.srf2map(f, gmt_temp, prefix = i))
plot_file = '%s/fault_plot.ps' % (gmt_temp)
out_dir = os.path.abspath('Png')

si_region = (166.2, 174.6, -47.4, -40.4)
ni_region = (172.4, 178.8, -41.8, -34.2)
# width of south island map, calculate height given projection, region, size
si_width = 6
si_height = gmt.mapproject(si_region[0], si_region[3], \
        wd = gmt_temp, projection = 'M%s' % (si_width), \
        region = si_region)[1]
# space between south and north island map (excludes tick marks/labels)
i_gap = 1
# calculate width and height of north island map to match height
ni_width, ni_height = gmt.map_width('M', si_height, ni_region, \
        wd = gmt_temp, accuracy = 0.001)
# functions to go back and forth between plots
def teleport(plot, location):
    if location == 'si':
        plot.spacial('M', si_region, sizing = si_width, \
                x_shift = -si_width - i_gap)
    elif location == 'ni':
        plot.spacial('M', ni_region, sizing = ni_width, \
                x_shift = si_width + i_gap)

# prepare topography colour palette
cpt_land = '%s/land.cpt' % (gmt_temp)
gmt.makecpt('%s/cpt/palm_springs_1.cpt' % (script_dir), cpt_land, \
        -250, 9000, inc = 10, invert = True)

left_margin = 1
right_margin = 1
bottom_margin = 1
top_margin = 1.5
p = gmt.GMTPlot(plot_file)
p.background(left_margin + si_width + i_gap + ni_width + right_margin, \
        bottom_margin + max(si_height, ni_height) + top_margin, \
        colour = 'gray')
# add title
p.text(left_margin + (si_width + i_gap + ni_width) / 2, \
        bottom_margin + max(si_height, ni_height) + 0.5, \
        title, size = 28, align = 'CB')
# apply initial offset to leave space on page for tickmarks etc..
p.spacial('X', (0, 1, 0, 1), x_shift = left_margin, y_shift = bottom_margin)
# draw blank maps
for i in ['ni', 'si']:
    teleport(p, i)
    if i == 'ni':
        region = ni_region
    elif i == 'si':
        region = si_region
    p.land(fill = 'darkgreen')
    p.topo('/home/nesi00213/PlottingData/Topo/srtm_all_filt_nz.grd', \
            cpt = cpt_land)
    p.water()
    p.coastlines(width = 0.2)
    # active faults in NZ as thin red lines
    p.path(faults, is_file = True, close = False, \
            width = '0.4p', colour = 'red')
    # fault planes before beachballs which are smaller
    for f in faults_historic:
        p.fault(f, is_srf = False, hyp_size = 0.2, \
                plane_width = 0.8, top_width = 1.2, hyp_width = 0.8, \
                plane_colour = 'black', top_colour = 'black', \
                hyp_colour = 'black')
    for f in faults_future:
        p.fault(f, is_srf = False, hyp_size = 0.2, \
                plane_width = 0.8, top_width = 1.2, hyp_width = 0.8, \
                plane_colour = 'blue', top_colour = 'blue', hyp_colour = 'blue')
    # beachballs on top as they are smaller than fault planes
    p.beachballs('/nesi/projects/nesi00213/PlottingData/Earthquakes/CMTData_Mw3p5_5_20170118_Cant_yesFtp.meca', \
            is_file = True, scale = 0.05, colour = 'black')
    p.beachballs('/nesi/projects/nesi00213/PlottingData/Earthquakes/Hoby.meca', \
            is_file = True, scale = 0.05, colour = 'black')
    p.beachballs('/nesi/projects/nesi00213/PlottingData/Earthquakes/Ahsan.meca', \
            is_file = True, scale = 0.05, colour = 'blue')
    if i == 'si' or i == 'ni':
        # loop through srfs and planes
        for s in xrange(len(srf_data)):
            for plane in xrange(len(srf_data[s][1])):
                p.overlay('%s/%d_%d_slip.bin' % (gmt_temp, s, plane), \
                        '%s/%d.cpt' % (gmt_temp, s), \
                        dx = srf_data[s][0][0], dy = srf_data[s][0][1], \
                        climit = 2, \
                        crop_grd = '%s/%d_%d_mask.grd' % (gmt_temp, s, plane), \
                        land_crop = False, transparency = 35, \
                        custom_region = srf_data[s][1][plane])
    p.ticks(major = 2, minor = 0.2)

    # display labels, allow outside region if they started inside the region
    for label in labels:
        if region[0] <= label[0] < region[1] \
                and region[2] <= label[1] <= region[3]:
            p.text(label[0], label[1], label[3], align = label[2], \
                    clip = False, box_fill = 'white@30')
# add QuakeCoRE logo
p.image(si_region[0], si_region[3], '%s/quakecore-logo.png' % (script_dir), \
        width = '3i', align = 'LT')
# ideal command but doesn't work with older GMT 5.1
#p.image('L', 'T', '%s/quakecore-logo.png' % (script_dir), \
#        width = '3i', pos = 'rel')

p.finalise()
p.png(dpi = dpi, clip = True)
