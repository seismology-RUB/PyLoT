#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
from obspy import UTCDateTime

def gse_time_header(lines):
    '''
    takes a path FILE to a GSE data file and returns the time header of the file
    :param file: path to GSE data file
    :type file: str
    :return: time header from FILE
    :rtype: str

    >>> gse_time_header('test.gse')
    "WID2 2005/10/09 20:17:25.000 ATH   SHZ  NSP CM6     9000   50.000000   0.10E+01   1.000    NSP  -1.0  0.0"
    '''

    return lines[1]

def time_from_header(header):
    timeline = header.split(' ')
    time = timeline[1].split('/') + timeline[2].split(':')
    time = time[:-1] + time[-1].split('.')
    time[-1] += '000'
    return [int(t) for t in time]

def check_time(time):
    try:
        UTCDateTime(time)
        return True
    except ValueError:
        return False