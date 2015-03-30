#!/usr/bin/env python
# -*- coding: utf-8 -*-

from obspy.signal.trigger import recSTALTA, triggerOnset


def createSingleTriggerlist(st, station='ZV01', trigcomp='Z', stalta=(1, 10),
                            trigonoff=(6, 1)):
    '''
    uses a single-station trigger to create a triggerlist for this station
    :param st:
    :param station:
    :param trigcomp:
    :param stalta:
    :param trigonoff:
    :return:
    '''
    tr = st.copy().select(component=trigcomp, station=station)[0]
    df = tr.stats.sampling_rate

    cft = recSTALTA(tr.data, int(stalta[0] * df), int(stalta[1] * df))
    triggers = triggerOnset(cft, trigonoff[0], trigonoff[1])
    trigg = []
    for time in triggers:
        trigg.append(tr.stats.starttime + time[0] / df)
    return trigg


def createSubCoincTriggerlist(trig, station='ZV01'):
    '''
    makes a triggerlist with the events, that are triggered by the
    coincidence trigger and are seen at the demanded station
    :param trig: list containing triggers from coincidence trigger
    :type trig: list
    :param station: station name to get triggers for
    :type station: str
    :return: list of triggertimes
    :rtype: list
    '''
    trigg = []
    for tri in trig:
        if station in tri['stations']:
            trigg.append(tri['time'])
    return trigg
