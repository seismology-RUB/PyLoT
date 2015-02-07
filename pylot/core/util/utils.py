#!/usr/bin/env python
#
# -*- coding: utf-8 -*-

import os
import pwd
import re
import obspy.core.event as ope

def fnConstructor(s):

    s = s.split('/')[-1]

    badchars = re.compile(r'[^A-Za-z0-9_. ]+|^\.|\.$|^ | $|^$')
    badsuffix = re.compile(r'(aux|com[1-9]|con|lpt[1-9]|prn)(\.|$)')

    fn = badchars.sub('_', s)

    if badsuffix.match(fn):
        fn = '_' + fn
    return fn

def createEvent(origintime, latitude, longitude, depth, **kwargs):
    evt = ope.Event()

def createArrival(picknum, picktime, eventnum, station, cinfo, phase, wfname,
                  authority_id):
    pickID = 'pick/' + eventnum + '/' + station + '/{0:3d}'.format(picknum)
    pickresID = ope.ResourceIdentifier(id=pickID)
    pickresID.convertIDToQuakeMLURI(authority_id=authority_id)
    pick = ope.Pick()
    pick.resource_id = pickresID
    pick.time = picktime
    pick.creation_info = cinfo
    pick.phase_hint = phase
    pick.waveform_id = ope.ResourceIdentifier(id=wfname, prefix='file:/')

    arriID = 'arrival/' + eventnum + '/' + station + '/{0}'.format(phase)
    arriresID = ope.ResourceIdentifier(id=arriID)
    arriresID.convertIDToQuakeMLURI(authority_id=authority_id)
    arrival = ope.Arrival()
    arrival.resource_id = arriresID
    arrival.creation_info = cinfo
    arrival.pick_id = pickresID
    arrival.phase = pick.phase_hint
    azi = self.location[eventid]['Backazimuth'] - 180
    arrival.azimuth = azi if azi > -180 else azi + 360
    arrival.distance = self.location[eventid]['Distance']['deg']

def getOwner(fn):
    return pwd.getpwuid(os.stat(fn).st_uid).pw_name



