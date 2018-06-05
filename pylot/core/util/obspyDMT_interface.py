#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from obspy import UTCDateTime

def check_obspydmt_structure(path):
    '''
    Check path for obspyDMT event structure.
    :param path:
    :return:
    '''
    ev_info = os.path.join(path, 'EVENTS-INFO')
    if os.path.isdir(ev_info):
        if os.path.isfile(os.path.join(ev_info, 'logger_command.txt')):
            return True
    return False

def check_obspydmt_eventfolder(folder):
    try:
        time = folder.split('.')[0]
        time = time.replace('_', 'T')
        time = UTCDateTime(time)
        return True, time
    except Exception as e:
        return False, e

def qml_from_obspyDMT(path):
    import pickle
    from obspy.core.event import Event, Magnitude, Origin

    if not os.path.exists(path):
        return IOError('Could not find Event at {}'.format(path))
    infile = open(path, 'rb')
    event_dmt = pickle.load(infile)
    ev = Event(resource_id=event_dmt['event_id'])
    origin = Origin(resource_id=event_dmt['origin_id'], time=event_dmt['datetime'], longitude=event_dmt['longitude'],
                    latitude=event_dmt['latitude'], depth=event_dmt['depth'])
    mag = Magnitude(mag=event_dmt['magnitude'], magnitude_type=event_dmt['magnitude_type'],
                    origin_id=event_dmt['origin_id'])
    ev.magnitudes.append(mag)
    ev.origins.append(origin)
    return ev

