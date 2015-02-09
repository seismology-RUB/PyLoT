#!/usr/bin/env python
#
# -*- coding: utf-8 -*-

import os
import pwd
import re
import hashlib
from obspy.core import UTCDateTime
import obspy.core.event as ope

def fnConstructor(s):

    s = s.split('/')[-1]

    badchars = re.compile(r'[^A-Za-z0-9_. ]+|^\.|\.$|^ | $|^$')
    badsuffix = re.compile(r'(aux|com[1-9]|con|lpt[1-9]|prn)(\.|$)')

    fn = badchars.sub('_', s)

    if badsuffix.match(fn):
        fn = '_' + fn
    return fn

def getHash(origintime):
    '''

    :param origintime:
    :return:
    '''
    hg = hashlib.sha1()
    hg.update(origintime.strftime('%Y-%m-%d %H:%M:%S.%f'))
    return hg.hexdigest()

def createResourceID(timetohash, restype, authority_id=None):
    '''

    :param timetohash:
    :param restype:
    :param authority_id:
    :return:
    '''
    assert isinstance(timetohash, UTCDateTime), "'timetohash' is not an ObsPy" \
                                                "UTCDateTime object"
    hid = getHash(timetohash)
    resID = ope.ResourceIdentifier(restype + '/' + hid[0:6])
    if authority_id is not None:
        resID.convertIDToQuakeMLURI(authority_id=authority_id)


def createOrigin(origintime, latitude, longitude, depth, resID=None,
                 authority_id=None):
    '''
    createOrigin - function to create an ObsPy Origin
    :param origintime: the origins time of occurence
    :type origintime: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :param latitude: latitude in decimal degree of the origins location
    :type latitude: float
    :param longitude: longitude in decimal degree of the origins location
    :type longitude: float
    :param depth: hypocentral depth of the origin
    :type depth: float
    :return: An ObsPy :class: `~obspy.core.event.Origin` object
    '''
    pass


def createEvent(origintime, cinfo, etype, resID=None,
                authority_id=None, **kwargs):
    '''
    createEvent - funtion to create an ObsPy Event
    :param origintime: the events origintime
    :type origintime: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :param cinfo:
    :param etype:
    :param resID:
    :param authority_id:
    :param kwargs:
    :return: An ObsPy :class: `~obspy.core.event.Event` object
    '''
    etype = ope.EventType(etype)
    if resID is None:
        resID = createResourceID(origintime, etype, authority_id)
    event = ope.Event(resource_id=resID)
    event.creation_info = cinfo
    event.event_type = etype
    return event

def createPick(picknum, picktime, eventnum, cinfo, phase, station, wfseedstr,
               authority_id):
    '''
    createPick - function to create an ObsPy Pick

    :param picknum: number of the created pick
    :type picknum: int
    :param eventnum: human-readable event identifier
    :type eventnum: str
    :param cinfo: An ObsPy :class: `~obspy.core.event.CreationInfo` object
        holding information on the creation of the returned object
    :type cinfo: :class: `~obspy.core.event.CreationInfo` object
    :param phase: name of the arrivals seismic phase
    :type phase: str
    :param station: name of the station at which the seismic phase has been
        picked
    :type station: str
    :param wfseedstr: A SEED formatted string of the form
        network.station.location.channel in order to set a referenced waveform
    :type wfseedstr: str, SEED formatted
    :param authority_id: name of the institution carrying out the processing
    :type authority_id: str
    :return: An ObsPy :class: `~obspy.core.event.Pick` object
    '''
    pickID = 'pick/' + eventnum + '/' + station + '/{0:3d}'.format(picknum)
    pickresID = ope.ResourceIdentifier(id=pickID)
    pickresID.convertIDToQuakeMLURI(authority_id=authority_id)
    pick = ope.Pick()
    pick.resource_id = pickresID
    pick.time = picktime
    pick.creation_info = cinfo
    pick.phase_hint = phase
    pick.waveform_id = ope.ResourceIdentifier(id=wfseedstr, prefix='file:/')
    return pick

def createArrival(pickresID, eventnum, cinfo, phase, station, authority_id,
                  azimuth=None, dist=None):
    '''
    createArrival - function to create an Obspy Arrival
    :param pickresID: Resource identifier of the created pick
    :type pickresID: :class: `~obspy.core.event.ResourceIdentifier` object
    :param eventnum: human-readable event identifier
    :type eventnum: str
    :param cinfo: An ObsPy :class: `~obspy.core.event.CreationInfo` object
    holding information on the creation of the returned object
    :type cinfo: :class: `~obspy.core.event.CreationInfo` object
    :param phase: name of the arrivals seismic phase
    :type phase: str
    :param station: name of the station at which the seismic phase has been
    picked
    :type station: str
    :param authority_id: name of the institution carrying out the processing
    :type authority_id: str
    :param azimuth: azimuth between source and receiver
    :type azimuth: float or int, optional
    :param dist: distance between source and receiver
    :type dist: float or int, optional
    :return: An ObsPy :class: `~obspy.core.event.Arrival` object
    '''
    arriID = 'arrival/' + eventnum + '/' + station + '/{0}'.format(phase)
    arriresID = ope.ResourceIdentifier(id=arriID)
    arriresID.convertIDToQuakeMLURI(authority_id=authority_id)
    arrival = ope.Arrival()
    arrival.resource_id = arriresID
    arrival.creation_info = cinfo
    arrival.pick_id = pickresID
    arrival.phase = phase
    if azimuth is not None:
        arrival.azimuth = float(azimuth) if azimuth > -180 else azimuth + 360
    else:
        arrival.azimuth = azimuth
    arrival.distance = None
    return arrival

def getOwner(fn):
    return pwd.getpwuid(os.stat(fn).st_uid).pw_name



