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

check_obspydmt_eventfolder('20110311_054623.a')