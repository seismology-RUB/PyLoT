#!/usr/bin/env python
# -*- coding: utf-8 -*-

from obspy.signal.trigger import recursive_sta_lta, trigger_onset


def createSingleTriggerlist(st, station='ZV01', trigcomp='Z', stalta=(1, 10),
                            trigonoff=(6, 1)):
    '''
    uses a single-station trigger to create a triggerlist for this station
    :param st: obspy stream
    :type  st:
    :param station: station name to get triggers for  (optional, default = ZV01)
    :type  station: str
    :param trigcomp: (optional, default = Z)
    :type  trigcomp: str
    :param stalta: (optional, default = (1,10))
    :type  stalta: tuple
    :param trigonoff: (optional, default = (6,1))
    :type  trigonoff: tuple
    :return: list of triggtimes
    :rtype: list
    '''
    tr = st.copy().select(component=trigcomp, station=station)[0]
    df = tr.stats.sampling_rate

    cft = recursive_sta_lta(tr.data, int(stalta[0] * df), int(stalta[1] * df))
    triggers = trigger_onset(cft, trigonoff[0], trigonoff[1])
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
    :param station: station name to get triggers for (optional, default = ZV01)
    :type station: str
    :return: list of triggertimes
    :rtype: list
    '''
    trigg = []
    for tri in trig:
        if station in tri['stations']:
            trigg.append(tri['time'])
    return trigg
