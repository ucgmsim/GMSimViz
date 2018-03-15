#!/usr/bin/env python2
"""
Plot faults that have been studied by QuakeCoRE
"""

from glob import glob
import os

import qcore_path
import geo
import gmt
import srf

#title = 'NZ ERF (Stirling et al., 2012) - Shallow crustal ruptures'
#title = 'NZ ERF (Stirling et al., 2012) - Subduction interface ruptures'
title = ''

faults_folder = '/home/nesi00213/PlottingData/Earthquakes'
faults_historic = ['13Juneb2011FaultPlane.xy', \
         \
        '22Feb2011FaultPlane.xy', \
        '23Deca2011FaultPlane.xy', '23Decb2011FaultPlane.xy', \
        '20161114 Kaikoura Ian02_s103245.txt', \
        'Cook Strait 012017.txt', 'DarfieldFaultPlane.xy', \
        'Wilberforce 122016.txt', '2014Jan20_Eketahuna_m6p3.txt']
faults_historic = []
faults_future = ['Alpine Fault m7.90-411.0x17.3_s1129570.txt', \
        'Porters Pass M7.20_82.0by14.0_s83747.txt']
faults_future = []
srf_files = ['20161114 Kaikoura Ian02_s103245.srf', \
        '20110222 m6p2.srf', '20100904 m7p1.srf',
        '20110613b m6.00-13.0x9.0_s19126.srf', \
        'Alpine Fault m7.90-411.0x17.3_s1129570.srf', \
        'Porters Pass M7.20_82.0by14.0_s83747.srf', \
        'Wilberforce standard_m5.60-5.1x5.1_s103245.srf', \
        'Cook Strait 012017.srf', \
        '2014Jan20_Eketahuna_m6p3.srf']
srf_files = []
faults = '/nesi/projects/nesi00213/PlottingData/Paths/faults/FAULTS_20161219.ll'
# convert to absolute paths
rel2abs = lambda fault_file : os.path.join(faults_folder, fault_file)
faults_historic = map(rel2abs, faults_historic)
faults_future = map(rel2abs, faults_future)
srf_files = map(rel2abs, srf_files)
import sys
srf_files = glob('/home/vap30/scratch/combined_sth_srf/*/Srf/*_HYP01-*1244.srf')
corner_files = glob('/home/vap30/scratch/combined_sth_vm/*/VeloModCorners.txt')

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
labels = []

dpi = 600
script_dir = os.path.abspath(os.path.dirname(__file__))
gmt_temp = os.path.abspath('GMT_WD_FAULTS2')
if not os.path.exists(gmt_temp):
    os.makedirs(gmt_temp)
srf_data = []
good_srf = 0

WANTED = ['ACTIVE_SHALLOW', 'VOLCANIC']
#WANTED = ['SUBDUCTION_INTERFACE']
nhm = '/home/vap30/ucgmsim/Pre-processing/SrfGen/NHM/NZ_FLTmodel_2010.txt'
def nhmtype(name):
    with open(nhm, 'r') as nf:
        db = nf.readlines()
        dbi = 15
        dbl = len(db)
        while dbi < dbl:
            if db[dbi].strip() == name:
                return db[dbi + 1].split()[0]
            dbi += 13 + int(db[dbi + 11])
        return 'NOTFOUND'

# VM bounds
corner_paths = []
if not os.path.isdir('%s/vm_corners' % (gmt_temp)):
    os.makedirs('%s/vm_corners' % (gmt_temp))
for c in corner_files:
    with open(c, 'r') as cf:
        cnrs = [map(float, ll) for ll in map(str.split, cf.readlines()[2:])]
    hr_cnrs = geo.path_from_corners(corners = cnrs, output = None, close = True)
    gmt_cnrs = '\n'.join([' '.join(map(str, ll)) for ll in hr_cnrs])
    name = os.path.basename(os.path.abspath(os.path.dirname(c)))
    with open('%s/vm_corners/%s.txt' % (gmt_temp, name), 'w') as out:
        out.write(gmt_cnrs)
    if nhmtype(name) in WANTED:
        corner_paths.append('%s/vm_corners/%s.txt' % (gmt_temp, name))

#bad_cnrs = []
#good_cnrs = []
for i, f in enumerate(srf_files):
    name = os.path.basename(f).split('_')[0]
    if nhmtype(name) not in WANTED:
            #or (name[-1] in ['2', '3', '4', '5', '6', '7', '8', '9'] and name not in ['NukWaitot1to6', 'Hope1888']) \
            #or (name[-2] == '0' and name[-1] != '1') \
            #or (name[-2] in ['1', '2', '3'] and name[-1] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'] and name not in ['RaukumaraF13']):
    #    print('SKIPPING %s' % (name))
        #cnrs = '%s/bad_%d.cnrs' % (gmt_temp, len(bad_cnrs))
        try:
            print 'bad corners'
            srf.srf2corners(f, cnrs = '%s/%s.bad_cnrs' % (gmt_temp, i))
        except ValueError:
            pass
        continue
        #bad_cnrs.append(cnrs)
    #    continue
    try:
        srf_data.append(gmt.srf2map(f, gmt_temp, prefix = good_srf, \
                wd = gmt_temp))
        cnrs = '%s/%d.cnrs' % (gmt_temp, good_srf)
        srf.srf2corners(f, cnrs = cnrs)
        #good_cnrs.append(cnrs)
        good_srf += 1
    except ValueError:
        pass
plot_file = '%s/fault_plot.ps' % (gmt_temp)
out_dir = os.path.abspath('Png')

si_region = (166.2, 176, -47.4, -40)
si_region = (165, 176, -47.4, -40)
ni_region = (172.4, 178.8, -41.8, -34.2)
ni_region = (171.4, 179.8, -41.8, -34.2)
ni_region = (166.2, 176, -47.4, -40)
# width of south island map, calculate height given projection, region, size
si_width = 0
si_height = 0#gmt.mapproject(si_region[0], si_region[3], \
        #wd = gmt_temp, projection = 'M%s' % (si_width), \
        #region = si_region)[1]
# space between south and north island map (excludes tick marks/labels)
i_gap = 0
# calculate width and height of north island map to match height
ni_width = 9
ni_height = gmt.mapproject(ni_region[0], ni_region[3], \
        wd = gmt_temp, projection = 'M%s' % (ni_width), \
        region = ni_region)[1]
#ni_width, ni_height = 0,0#gmt.map_width('M', si_height, ni_region, \
#        wd = gmt_temp, accuracy = 0.001)
# functions to go back and forth between plots
def teleport(plot, location):
    #if location == 'si':
    #    plot.spacial('M', si_region, sizing = si_width, \
    #            x_shift = -0 - i_gap)
    #            #x_shift = -si_width - i_gap)
    if location == 'ni':
        plot.spacial('M', ni_region, sizing = ni_width, \
                x_shift = 0)#si_width + i_gap)

left_margin = 1
right_margin = 1
bottom_margin = 1.5
top_margin = 1.2
gmt.gmt_defaults(wd = gmt_temp, ps_media = 'Custom_%six%si' \
        % (left_margin + si_width + i_gap + ni_width + right_margin, \
        bottom_margin + max(si_height, ni_height) + top_margin))
p = gmt.GMTPlot(plot_file, reset = False)
p.background(left_margin + si_width + i_gap + ni_width + right_margin, \
        bottom_margin + max(si_height, ni_height) + top_margin, \
        colour = 'white')
# add title
p.text(left_margin + (si_width + i_gap + ni_width) / 2, \
        bottom_margin + max(si_height, ni_height) + 0.5, \
        title, size = 28, align = 'CB')
# apply initial offset to leave space on page for tickmarks etc..
p.spacial('X', (0, 1, 0, 1), x_shift = left_margin, y_shift = bottom_margin)
# draw blank maps
for i in ['ni']:
    teleport(p, i)
    if i == 'ni':
        region = ni_region
    elif i == 'si':
        region = si_region
    p.basemap(topo_cpt = 'grey1')
    ###
    ### VM CORNERS
    ###
    for path in corner_paths:
        p.path(path, is_file = True, close = True, split = '-')
    # non-uniform
    #p.clip(path = corner_paths, is_file = True)
    #p.points('/home/vap30/ucgmsim/Pre-processing/NonUniformGrid/17.06.2017/non_uniform_whole_nz-hh400.ll', shape = 'c', size = '0.8p', fill = 'black@30', line = None)
    #p.clip()
    # active faults in NZ as thin red lines
    #p.path(faults, is_file = True, close = False, \
    #        width = '0.2p', colour = 'red')
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
    #p.beachballs('/nesi/projects/nesi00213/PlottingData/Earthquakes/CMTData_Mw3p5_5_20170118_Cant_yesFtp.meca', \
    #        is_file = True, scale = 0.05, colour = 'black')
    #p.beachballs('/nesi/projects/nesi00213/PlottingData/Earthquakes/Hoby.meca', \
    #        is_file = True, scale = 0.05, colour = 'black')
    #p.beachballs('/nesi/projects/nesi00213/PlottingData/Earthquakes/Ahsan.meca', \
    #        is_file = True, scale = 0.05, colour = 'blue')
    if i == 'si' or i == 'ni':
        # loop through srfs and planes
        for s in xrange(len(srf_data)):
            for plane in xrange(len(srf_data[s][1])):
                p.overlay('%s/%d_%d_slip.bin' % (gmt_temp, s, plane), \
                        'ginc0.cpt', \
                        dx = srf_data[s][0][0], dy = srf_data[s][0][1], \
                        climit = 2, \
                        crop_grd = '%s/%d_%d_mask.grd' % (gmt_temp, s, plane), \
                        land_crop = False, transparency = 35, \
                        custom_region = srf_data[s][1][plane])
            p.fault('%s/%d.cnrs' % (gmt_temp, s), is_srf = False, \
                    hyp_size = 0, plane_width = 0.2, top_width = 0.4, \
                    hyp_width = 0.2, plane_colour = 'black', \
                    top_colour = 'black', hyp_colour = 'black')
            
        #for unwanted in glob('%s/*.bad_cnrs' % (gmt_temp)):
        #    p.fault(unwanted, is_srf = False, \
        #            hyp_size = 0, plane_width = 0.2, top_width = 0.4, \
        #            hyp_width = 0.2, plane_colour = 'red', \
        #            top_colour = 'red', hyp_colour = 'red')
        #for x in good_cnrs:
        #    p.fault(x, is_srf = False, \
        #            hyp_size = 0, plane_width = 0.2, top_width = 0.4, \
        #            hyp_width = 0.2, plane_colour = 'red', \
        #            top_colour = 'red', hyp_colour = 'red')
        #for x in bad_cnrs:
        #    p.fault(x, is_srf = False, \
        #            hyp_size = 0, plane_width = 0.2, top_width = 0.4, \
        #            hyp_width = 0.2, plane_colour = 'black', \
        #            top_colour = 'black', hyp_colour = 'black')
    p.ticks(major = 2, minor = 0.2)

    # display labels, allow outside region if they started inside the region
    for label in labels:
        if region[0] <= label[0] < region[1] \
                and region[2] <= label[1] <= region[3]:
            p.text(label[0], label[1], label[3], align = label[2], \
                    clip = False, box_fill = 'white@30')
# add QuakeCoRE logo
p.image(ni_region[0], ni_region[3], '%s/quakecore-logo.png' % (script_dir), \
        width = '3i', align = 'LT')
# ideal command but doesn't work with older GMT 5.1
#p.image('L', 'T', '%s/quakecore-logo.png' % (script_dir), \
#        width = '3i', pos = 'rel')
# scale
p.cpt_scale((ni_width + si_width + i_gap) / 2.0, -0.5, 'ginc0.cpt', \
        label = 'Slip (cm)', length = si_width + ni_width + i_gap, log = False)

p.finalise()
p.png(dpi = dpi, clip = False)
