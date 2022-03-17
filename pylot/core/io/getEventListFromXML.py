#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
   Script to get event parameters from PyLoT-xml file to write
   them into eventlist.
   LK, igem, 03/2021
   Edited for use in PyLoT
   JG, igem, 01/2022
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from obspy.core.event import read_events
from pyproj import Proj
import glob

"""
Creates an eventlist file summarizing all events found in a certain folder. Only called by pressing UI Button eventlis_xml_action

:rtype:
:param path: Path to root folder where single Event folder are to found
"""

def geteventlistfromxml(path, outpath):
    p = Proj(proj='utm', zone=32, ellps='WGS84')


    # open eventlist file and write header
    evlist = outpath + '/eventlist'
    evlistobj = open(evlist, 'w')
    evlistobj.write('EventID      Date      To    Lat    Lon     EAST    NORTH  Dep  Ml  NoP  NoS RMS  errH  errZ  Gap \n')

    # data path
    dp = path + "/e*/*.xml"
    # list of all available xml-files
    xmlnames = glob.glob(dp)

    # read all onset weights
    for names in xmlnames:
        print("Getting location parameters from {}".format(names))
        cat = read_events(names)
        try:
            st = cat.events[0].origins[0].time
            Lat = cat.events[0].origins[0].latitude
            Lon = cat.events[0].origins[0].longitude
            EAST, NORTH = p(Lon, Lat)
            Dep = cat.events[0].origins[0].depth / 1000
            Ml = cat.events[0].magnitudes[1].mag
            NoP = []
            NoS = []
        except IndexError:
            print ('Insufficient data found for event (not localised): ' + names.split('/')[-1].split('_')[-1][:-4] + ' Skipping event for eventlist.' )
            continue
 
        for i in range(len(cat.events[0].origins[0].arrivals)):
            if cat.events[0].origins[0].arrivals[i].phase == 'P':
                NoP.append(cat.events[0].origins[0].arrivals[i].phase)
            elif cat.events[0].origins[0].arrivals[i].phase == 'S':
                NoS.append(cat.events[0].origins[0].arrivals[i].phase)
        #NoP = cat.events[0].origins[0].quality.used_station_count
        errH = cat.events[0].origins[0].origin_uncertainty.max_horizontal_uncertainty
        errZ = cat.events[0].origins[0].depth_errors.uncertainty
        Gap = cat.events[0].origins[0].quality.azimuthal_gap
        #evID = names.split('/')[6]
        evID = names.split('/')[-1].split('_')[-1][:-4]
        Date = str(st.year) + str('%02d' % st.month) + str('%02d' % st.day)
        To = str('%02d' % st.hour) + str('%02d' % st.minute) + str('%02d' % st.second) + \
                  '.' + str('%06d' % st.microsecond)

        # write into eventlist
        evlistobj.write('%s %s %s %9.6f %9.6f %13.6f %13.6f %8.6f %3.1f %d %d NaN %d %d %d\n' %(evID, \
                         Date, To, Lat, Lon, EAST, NORTH, Dep, Ml, len(NoP), len(NoS), errH, errZ, Gap))
        print ('Adding Event ' + names.split('/')[-1].split('_')[-1][:-4]  + ' to eventlist')
    print('Eventlist created and saved in: ' + outpath)
    evlistobj.close()

