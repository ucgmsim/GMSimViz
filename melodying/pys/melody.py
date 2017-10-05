#!/usr/bin/env python
import subprocess
import gmt
import os
import shutil
import glob
COlOR_DICT = {0.0:'red',0.1:'red',1.0:'green',0.8:'yellow',0.9:'orange',0.95:'limegreen',1.1:"red",0.81:"navy",0.91:"red",1.5:"magenta",0.5:"pink",0.75:"blue"}
REGION = (170, 177, -45, -40)
CHCH = (172.639847,-43.525650)
# region = (170, 177, -45, -40)
#
# p = gmt.GMTPlot('melody.ps')
#
# p.spacial('M', region, sizing = 7)
# p.basemap(res = 'f', topo = None)
#
# p.sites(['Wellington,LM'])
#
# p.points('171 -41\n172 -40.5', is_file = False, shape = 'c', size = 0.5)
# p.path('171 -41\n172 -40.5', is_file = False, width = '1c', colour = 'red@30')
#
#
# p.finalise()
# p.png(background = 'white')
# ROAD:  ogr2ogr -f "GMT" 3_15_Nov_2016.gmt 3_15_Nov_2016.shx -s_srs EPSG:2193 -t_srs EPSG:4326
# FLIGHT:  ogr2ogr -f "GMT" 3_15_Nov_2016.gmt 3_15_Nov_2016.shx -s_srs EPSG:2193


# p.path('171 -41\n172 -40.5', is_file = False, width = '1c', colour = 'red@30')

def cp_and_rename(pattern, dir):
    """pattern:common pattern of the images you want to copy. eg 'pic*.png'
       dir: path to images with the specified pattern, default current folder
    """
    for item in glob.glob1(dir,pattern):
        print("item",item)
        s = os.path.join(dir, item)
        print("s",s)
        shutil.copy(s, os.path.join(dir, "cp{}".format(item)))


def custom_make_movie(input_pattern,output_name, fps=1 ,pix_fmt="yuv420p", dir='.'):
    """input_pattern:common pattern of the images you want to copy. eg 'pic*.png'
       output_name: name of output movie
       fps:frame per second
       pix_fmt:resolution of the movie
       dir: path to images with the specified pattern, default current folder
    """
    cp_and_rename(input_pattern,dir)
    input_pattern = "*"+input_pattern
    subprocess.call(
        ["ffmpeg", "-framerate", str(fps), "-pattern_type", "glob", "-i", dir+'/'+input_pattern, "-r", "30", "-pix_fmt", pix_fmt,
         output_name])


def check_width(z):
    """local roads(ends with floating digit 1 eg 0.91, 0.1, 1.1) are skinnier
        z: Los value
    """
    width = '0.14c'
    remainder = str(2-z)
    if '9' in remainder:
        width = '0.08c'
    return width


def set_start(even):
    """set start inedex
       even: to draw a even pos or not
    """
    start = 0 if even else -1
    return start


def set_offset(even):
    offset = 0 if even else 50
    return offset


def draw_base_path(pts,z,p):
    """draw base road of different colors for gmt png
        pts: point list of a road
        z: Los value
        p: GMTPlot object
    """
    strpts = ''.join(pts)  # gmt only takes a string of pts or a file
    p.path(strpts, is_file=False, width=check_width(z), colour=COlOR_DICT[z])
    p.path(strpts, is_file=False, width='0.2p', colour='white', straight=True, split='.-')


def place_vehicle(z, pts, start,offset,p):
    """show selected points(vehicle) on open roads.
       z: Los value
       pts: a list of all the points extracted from a given gmt file
       p: GMTPlot object
    """
    # if road open, show some points(vehicles) on top of the (open) roads
    if COlOR_DICT[z] != 'red':
        length = len(pts)
        print("length is ",length)
        if 610<length < 660:  #road appears longest has 628 points.  Road with most points is 748
            dots = pts[start] + pts[int(length / 8 * (start + 2))] + pts[int(length / 8 * (start + 4))] + pts[400+offset] + pts[int(length / 8 * (start + 6))] + pts[int(length / 8 * (start + 8))]
        elif length > 600:
            dots = pts[start] + pts[int(length / 8 * (start + 2))] + pts[int(length / 8 * (start + 4))] + \
                   pts[int(length / 8 * (start + 6))] + pts[int(length / 8 * (start + 8))]
        elif length > 300:
            dots = pts[start] + pts[length // 8 * (start + 4)] + pts[length // 8 * (start + 8)]
        elif length > 200:
            dots = pts[start] + pts[length // 8 * (start + 8)]
        else:
            dots = pts[start]

        # draw vehicles
        p.points(dots, is_file=False, shape='c',size='0.1', fill='black')
        print("draw points")



def read_road(file,p,even=True):
    """preliminary verison of drawing base path and vehicles, two chunks of functions needs to be merged"""
    # open a gmt file for reading
    input_file = open(file,'r')
    lines = input_file.read()
    splitted_lines = lines.splitlines()

    # first, read all the points and draw the base path.  Base path has to be completely drawn before drawing any vehicle otherwise it might cover some vehicles.
    flag = False
    i = 0
    pts = []
    start = set_start(even)  # one-off check whether to draw start points or not
    offset = set_offset(even)
    for line in splitted_lines:
        i += 1
        if line.startswith("#"):
            if line.startswith("# @D"):  # if a new route
                # print("current line is {}".format(line))
                if flag:
                    # strpts = ''.join(pts)
                    draw_base_path(pts, z, p)  # draw base road
                    pts = []

                line = line.split("|")  # get single info out of line
                z = float(line[2])
                flag = True

        elif line.startswith(">"):  # still the same route
            continue

        else:
            x, y = line.split(" ")
            pt = x + ' ' + y + '\n'
            pts.append(pt)
            # strpts = ''.join(pts)
            if i == len(splitted_lines):
                draw_base_path(pts, z, p)

    # now we've finished drawing the base path, start drawing some vechicles.
    flag = False
    i = 0
    pts = []
    start = set_start(even)  # one-off check whether to draw start points or not
    offset=set_offset(even)
    for line in splitted_lines:
        i+=1
        if line.startswith("#"):
            if line.startswith("# @D"):   # if a new route
                # print("current line is {}".format(line))
                if flag:
                    # strpts = ''.join(pts)
                    # draw_base_path(pts, z, p)  # draw base roads
                    place_vehicle(z, pts, start, offset,p)
                    pts = []

                line=line.split("|")  # get single info out of line
                z=float(line[2])
                flag=True

        elif line.startswith(">"):   # still the same route
            continue

        else:
            x,y=line.split(" ")
            pt = x + ' ' + y + '\n'
            pts.append(pt)
            # strpts = ''.join(pts)
            if i == len(splitted_lines):
                # draw_base_path(pts, z, p)
                place_vehicle(z, pts, start, offset, p)


def draw_road(file, p, even=True):
    """Second verison of drawing base path and vehicles, can be potentially optimized?"""

    # open a gmt file for reading
    input_file = open(file, 'r')
    lines = input_file.read()
    splitted_lines = lines.splitlines()
    start = set_start(even)  # one-off check whether to draw start points or not
    offset = set_offset(even)  # one-off check to set the offset from the centre point in the visually longest road.
    completed=False
    base_drawn=False

    while not completed:
        flag = False  # flag for drawing the first road
        i = 0
        pts = []
        for line in splitted_lines:
            i += 1
            if line.startswith("#"):
                if line.startswith("# @D"):  # if a new route
                    if flag:
                        if not base_drawn:
                            draw_base_path(pts, z, p)  # draw base road
                        else:   # if base path completely drawn, can draw vehicles.
                            place_vehicle(z, pts, start, offset, p)
                        pts = []    # empty point list for the next road

                    line = line.split("|")  # get single info out of line
                    z = float(line[2])
                    flag = True

            elif line.startswith(">"):  # still the same route
                continue

            else:
                x, y = line.split(" ")
                pt = x + ' ' + y + '\n'
                pts.append(pt)
                if i == len(splitted_lines):
                    if not base_drawn:
                        draw_base_path(pts, z, p)
                        base_drawn=True   # reach end of file, base path finished
                    else:
                        place_vehicle(z, pts, start, offset, p)
                        completed = True  # reach end of file, vehicles finished. Both drawings finished.

# def read_road2(file,p,even=False):
#     input_file = open(file,'r')
#     flag=False
#     #points, vertices, zvalues = create_poly()
#     i = 0
#     lines = input_file.read()
#     pts=[]
#     #print(len(lines.splitlines()))
#     splitted_lines = lines.splitlines()
#     for line in splitted_lines:
#         i+=1
#         if line.startswith("#"):
#             if line.startswith("# @D"):   #if a new route
#                 # print("current line is {}".format(line))
#                 if flag:
#                     #poly_data=set_poly(points,vertices,zvalues)
#                     # print("points",pts)
#                     # strpts = ''.join(pts)
#                     # print("sgrprssssssssss",strpts)
#                     # p.path(strpts, is_file=False, width=check_width(z), colour='darkslategray')#'darkslategray'
#                     # p.path(strpts, is_file=False, width = '0.2p',colour='white',straight = True,split='.-')
#                     # if COlOR_DICT[z] != 'red':
#                     #     length = len(pts)
#                     #     if length>600:
#                     #         dots = pts[length // 8] + pts[length // 8 * 3] + pts[length // 8 * 5] + pts[length // 8 * 7] + pts[-1]
#                     #         print("dots",dots)
#                     #     elif 200 < length <= 600:
#                     #         dots = pts[length // 5]+pts[-1]
#                     #     else:
#                     #         dots = pts[-1]
#                     #
#                     #     p.points(dots, is_file=False, size='0.1', fill='black')
#                     place_vehicle(strpts, pts, p, z, even)
#                     pts = []
#
#                 line=line.split("|")
#                 z=float(line[2])
#                 flag=True
#
#         elif line.startswith(">"):   #still the same route
#             continue
#
#         else:
#             x,y=line.split(" ")
#             pt = x + ' ' + y + '\n'
#             pts.append(pt)
#             # strpts = ''.join(pts)
#             if i == len(splitted_lines):
#                 place_vehicle(strpts, pts, p, z, even)


def read_flight(file,p):
    input_file = open(file, 'r')
    flag = False
    i = 0
    lines = input_file.read()
    pts = ''

    splitted_lines = lines.splitlines()
    for line in splitted_lines:
        i += 1
        if line.startswith("#"):
            if line.startswith("# @D"):  # if a new route
                print("current line is {}".format(line))
                if flag:
                    p.path(pts, is_file=False, width=check_width(z), colour=COlOR_DICT[z],split='-')
                    pts = ''
                line = line.split("|")
                if line[1] == "":
                    filename = " ".join(line[3:])
                else:
                    filename = line[1]
                z = float(line[-1])
                print("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz", z)
                flag = True

        elif line.startswith(">"):  # still the same route
            continue

        else:
            # if i % 10 == 2:
            x, y = line.split(" ")
            pts += x + ' ' + y + '\n'
            if i == len(splitted_lines):
                p.path(pts, is_file=False, width=check_width(z), colour=COlOR_DICT[z],split='-')


def read_ferry(file,p):
    input_file = open(file, 'r')
    flag = False
    i = 0
    lines = input_file.read()
    pts = ''

    splitted_lines = lines.splitlines()
    for line in splitted_lines:
        i += 1
        if line.startswith("#"):
            if line.startswith("# @D"):  # if a new route
                print("current line is {}".format(line))
                if flag:
                    p.path(pts, is_file=False, width=check_width(z), colour=COlOR_DICT[z],split='.')
                    pts = ''
                line = line.split("|")
                if line[1] == "":
                    filename = " ".join(line[3:])
                else:
                    filename = line[1]
                z = float(line[-1])
                print("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz", z)
                flag = True

        elif line.startswith(">"):  # still the same route
            continue

        else:

            x, y = line.split(" ")
            pts += x + ' ' + y + '\n'
            if i == len(splitted_lines):
                p.path(pts, is_file=False, width=check_width(z), colour=COlOR_DICT[z],split='.')


def read_rail(file,p):
    input_file = open(file, 'r')
    flag = False
    i = 0
    lines = input_file.read()
    pts = ''

    splitted_lines = lines.splitlines()
    for line in splitted_lines:
        i += 1
        if line.startswith("#"):
            if line.startswith("# @D"):  # if a new route
                print("current line is {}".format(line))
                if flag:
                    p.path(pts, is_file=False, width=check_width(z), colour=COlOR_DICT[z])
                    pts = ''
                line = line.split("|")
                if line[1] == "":
                    filename = " ".join(line[3:])
                else:
                    filename = line[1]
                z = float(line[-1])
                print("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz", z)
                flag = True

        elif line.startswith(">"):  # still the same route
            continue

        else:
            x, y = line.split(" ")
            pts += x + ' ' + y + '\n'
            if i == len(splitted_lines):
                p.path(pts, is_file=False, width=check_width(z), colour=COlOR_DICT[z])

def read_all(folderpath):
    """input a path directing to a folder containing sub folders.
       Each sub folder contains gmt txt files of different types (road,flight,ferry,rail) with the same date.
       Then plots the desired paths"""

    for folder in os.listdir(folderpath):
        print("now,folder",folder)
        # create a new ps for each day
        p = gmt.GMTPlot(pspath='{}.ps'.format(folder))
        p.spacial('M', REGION, sizing=7)
        p.basemap(water='black', res='f', topo=None)
        p.sites(['Wellington,RM', 'Christchurch'])
        for file in os.listdir(folderpath+'/'+folder):
            print("now file",file)
            input_file = open(folderpath+'/'+folder+'/'+file, 'r')
            flag = False
            i = 0
            lines = input_file.read()
            pts = ''

            # set path style
            if 'flight' in file:
                split = '-'
            elif 'ferry' in file:
                split = '.'
            else:
                split = None

            splitted_lines = lines.splitlines()
            for line in splitted_lines:
                i += 1
                if line.startswith("#"):
                    if line.startswith("# @D"):  # if a new route
                        print("current line is {}".format(line))
                        if flag:
                            p.path(pts, is_file=False, width=check_width(z), colour=COlOR_DICT[z], split=split)
                            pts = ''
                        line = line.split("|")
                        if 'road' in file:  # position of Los/z is different in road gmt from other gmts
                            z = float(line[2])
                        else:
                            z = float(line[-1])
                        print("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz", z)
                        flag = True

                elif line.startswith(">"):  # still the same route
                    continue

                else:
                    x, y = line.split(" ")
                    pts += x + ' ' + y + '\n'
                    if i == len(splitted_lines):  # end of file
                        p.path(pts, is_file=False, width=check_width(z), colour=COlOR_DICT[z], split=split)
        # finalise a day
        p.finalise()
        p.png(background='darkblue')




def main():
    # print(CHCH[0])
    p = gmt.GMTPlot(pspath='draw_road_3_15_base.ps')
    # # p.image((REGION[0]+REGION[1])/2.0,(REGION[2]+REGION[3])/2.0,'tiny-pigs-002-08092014.jpg', align='C')
    # # print("adfasfdafdsafdafdssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss")
    p.spacial('M', REGION, sizing = 7)
    p.basemap(water='black',res = 'f', topo = None)
    p.sites(['Wellington,RM','Christchurch'])
    #
    file = '3_15_nov_2016.gmt'
    # # # flight_file = 'flight_14_Nov_2016.gmt'
    # # # # ferry_file = '14_Nov_2016_ferry.gmt'
    # # # # rail_file = 'rail_14_Nov_2016.gmt'

    #draw_road(file, p,even=True)
    # # # # read_ferry(ferry_file,p)
    # # # # read_flight(flight_file, p)
    # # # # read_rail(rail_file, p)
    # # # # read_all('/usr/local/qcore_data/yzh231/qcore/data')
    # # # # p.image(CHCH[0], CHCH[1], 'tiny-pigs-002-08092014.jpg', align='C')
    p.finalise()
    p.png(background='darkblue')
    # cp_and_rename('.', 'pic*.png')
    # custom_make_movie("pic*.png","cp.mp4")


if __name__ == '__main__':
    main()