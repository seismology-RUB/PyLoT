#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
   Script to get onset uncertainties from Quakeml.xml files created by PyLoT.
   Uncertainties are tranformed into quality classes and visualized via histogram if desired.
   Ludger Küperkoch, BESTEC GmbH, 07/2017
   rev.: Ludger Küperkoch, igem, 10/2020
   Edited for usage in PyLoT: Jeldrik Gaal, igem, 01/2022
"""

import glob

import matplotlib.pyplot as plt
import numpy as np
from obspy.core.event import read_events


def getQualitiesfromxml(path):
    # uncertainties
    ErrorsP = [0.02, 0.04, 0.08, 0.16]
    ErrorsS = [0.04, 0.08, 0.16, 0.32]

    Pw0 = []
    Pw1 = []
    Pw2 = []
    Pw3 = []
    Pw4 = []
    Sw0 = []
    Sw1 = []
    Sw2 = []
    Sw3 = []
    Sw4 = []

    # data path
    dp = path + '/e*/*.xml'
    # list of all available xml-files
    xmlnames = glob.glob(dp)

    # read all onset weights
    for names in xmlnames:
        print("Getting onset weights from {}".format(names))
        cat = read_events(names)
        arrivals = cat.events[0].picks
        for Pick in arrivals:
            if Pick.phase_hint[0] == 'P':
                if Pick.time_errors.uncertainty <= ErrorsP[0]:
                    Pw0.append(Pick.time_errors.uncertainty)
                elif Pick.time_errors.uncertainty > ErrorsP[0] and \
                        Pick.time_errors.uncertainty <= ErrorsP[1]:
                    Pw1.append(Pick.time_errors.uncertainty)
                elif Pick.time_errors.uncertainty > ErrorsP[1] and \
                        Pick.time_errors.uncertainty <= ErrorsP[2]:
                    Pw2.append(Pick.time_errors.uncertainty)
                elif Pick.time_errors.uncertainty > ErrorsP[2] and \
                        Pick.time_errors.uncertainty <= ErrorsP[3]:
                    Pw3.append(Pick.time_errors.uncertainty)
                elif Pick.time_errors.uncertainty > ErrorsP[3]:
                    Pw4.append(Pick.time_errors.uncertainty)
                else:
                    pass
            elif Pick.phase_hint[0] == 'S':
                if Pick.time_errors.uncertainty <= ErrorsS[0]:
                    Sw0.append(Pick.time_errors.uncertainty)
                elif Pick.time_errors.uncertainty > ErrorsS[0] and \
                        Pick.time_errors.uncertainty <= ErrorsS[1]:
                    Sw1.append(Pick.time_errors.uncertainty)
                elif Pick.time_errors.uncertainty > ErrorsS[1] and \
                        Pick.time_errors.uncertainty <= ErrorsS[2]:
                    Sw2.append(Pick.time_errors.uncertainty)
                elif Pick.time_errors.uncertainty > ErrorsS[2] and \
                        Pick.time_errors.uncertainty <= ErrorsS[3]:
                    Sw3.append(Pick.time_errors.uncertainty)
                elif Pick.time_errors.uncertainty > ErrorsS[3]:
                    Sw4.append(Pick.time_errors.uncertainty)
                else:
                    pass
            else:
                print("Phase hint not defined for picking!")
                pass
    # get percentage of weights
    numPweights = np.sum([len(Pw0), len(Pw1), len(Pw2), len(Pw3), len(Pw4)])
    numSweights = np.sum([len(Sw0), len(Sw1), len(Sw2), len(Sw3), len(Sw4)])
    try:
        P0perc = 100.0 / numPweights * len(Pw0)
    except:
        P0perc = 0
    try:
        P1perc = 100.0 / numPweights * len(Pw1)
    except:
        P1perc = 0
    try:
        P2perc = 100.0 / numPweights * len(Pw2)
    except:
        P2perc = 0
    try:
        P3perc = 100.0 / numPweights * len(Pw3)
    except:
        P3perc = 0
    try:
        P4perc = 100.0 / numPweights * len(Pw4)
    except:
        P4perc = 0
    try:
        S0perc = 100.0 / numSweights * len(Sw0)
    except:
        Soperc = 0
    try:
        S1perc = 100.0 / numSweights * len(Sw1)
    except:
        S1perc = 0
    try:
        S2perc = 100.0 / numSweights * len(Sw2)
    except:
        S2perc = 0
    try:
        S3perc = 100.0 / numSweights * len(Sw3)
    except:
        S3perc = 0
    try:
        S4perc = 100.0 / numSweights * len(Sw4)
    except:
        S4perc = 0

    weights = ('0', '1', '2', '3', '4')
    y_pos = np.arange(len(weights))
    width = 0.34
    p1, = plt.bar(0 - width, P0perc, width, color='black')
    p2, = plt.bar(0, S0perc, width, color='red')
    plt.bar(y_pos - width, [P0perc, P1perc, P2perc, P3perc, P4perc], width, color='black')
    plt.bar(y_pos, [S0perc, S1perc, S2perc, S3perc, S4perc], width, color='red')
    plt.ylabel('%')
    plt.xticks(y_pos, weights)
    plt.xlim([-0.5, 4.5])
    plt.xlabel('Qualities')
    plt.title('{0} P-Qualities, {1} S-Qualities'.format(numPweights, numSweights))
    plt.legend([p1, p2], ['P-Weights', 'S-Weights'])
    plt.show()
