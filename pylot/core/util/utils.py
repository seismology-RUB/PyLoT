#!/usr/bin/env python
#
# -*- coding: utf-8 -*-

import re
from obspy.core.event import *

def fnConstructor(s):

    s = s.split('/')[-1]

    badchars = re.compile(r'[^A-Za-z0-9_. ]+|^\.|\.$|^ | $|^$')
    badsuffix = re.compile(r'(aux|com[1-9]|con|lpt[1-9]|prn)(\.|$)')

    fn = badchars.sub('_', s)

    if badsuffix.match(fn):
        fn = '_' + fn
    return fn

def createEvent(origintime, latitude, longitude, depth, **kwargs):
    evt = Event()


