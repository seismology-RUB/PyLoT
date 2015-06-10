#!/usr/bin/env python
#
# -*- coding: utf-8 -*-

import os
import pwd
import re
import hashlib
import numpy as np
from obspy.core import UTCDateTime
import obspy.core.event as ope


def fnConstructor(s):
    if type(s) is str:
        s = s.split('/')[-1]
    else:
        s = getHash(UTCDateTime())

    badchars = re.compile(r'[^A-Za-z0-9_. ]+|^\.|\.$|^ | $|^$')
    badsuffix = re.compile(r'(aux|com[1-9]|con|lpt[1-9]|prn)(\.|$)')

    fn = badchars.sub('_', s)

    if badsuffix.match(fn):
        fn = '_' + fn
    return fn


def getLogin():
    return pwd.getpwuid(os.getuid())[0]


def getHash(time):
    '''
    :param time: time object for which a hash should be calculated
    :type time: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :return: str
    '''
    hg = hashlib.sha1()
    hg.update(time.strftime('%Y-%m-%d %H:%M:%S.%f'))
    return hg.hexdigest()


def getOwner(fn):
    return pwd.getpwuid(os.stat(fn).st_uid).pw_name


def prepTimeAxis(stime, trace):
    nsamp = trace.stats.npts
    srate = trace.stats.sampling_rate
    tincr = trace.stats.delta
    etime = stime + nsamp / srate
    time_ax = np.arange(stime, etime, tincr)
    if len(time_ax) < nsamp:
        print 'elongate time axes by one datum'
        time_ax = np.arange(stime, etime + tincr, tincr)
    elif len(time_ax) > nsamp:
        print 'shorten time axes by one datum'
        time_ax = np.arange(stime, etime - tincr, tincr)
    if len(time_ax) != nsamp:
        raise ValueError('{0} samples of data \n '
                         '{1} length of time vector \n'
                         'delta: {2}'.format(nsamp, len(time_ax), tincr))
    return time_ax


def getGlobalTimes(stream):
    min_start = UTCDateTime()
    max_end = None
    for trace in stream:
        if trace.stats.starttime < min_start:
            min_start = trace.stats.starttime
            if max_end is None or trace.stats.endtime > max_end:
                max_end = trace.stats.endtime
    return [min_start, max_end]


def createCreationInfo(agency_id=None, creation_time=None, author=None):
    if author is None:
        author = getLogin()
    if creation_time is None:
        creation_time = UTCDateTime()
    return ope.CreationInfo(agency_id=agency_id, author=author,
                            creation_time=creation_time)


def createResourceID(timetohash, restype, authority_id=None, hrstr=None):
    '''

    :param timetohash:
    :param restype: type of the resource, e.g. 'orig', 'earthquake' ...
    :type restype: str
    :param authority_id: name of the institution carrying out the processing
    :type authority_id: str, optional
    :return:
    '''
    assert isinstance(timetohash, UTCDateTime), "'timetohash' is not an ObsPy" \
                                                "UTCDateTime object"
    hid = getHash(timetohash)
    if hrstr is None:
        resID = ope.ResourceIdentifier(restype + '/' + hid[0:6])
    else:
        resID = ope.ResourceIdentifier(restype + '/' + hrstr + '_' + hid[0:6])
    if authority_id is not None:
        resID.convertIDToQuakeMLURI(authority_id=authority_id)
    return resID


def createOrigin(origintime, cinfo, latitude, longitude, depth, resID=None,
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
    if resID is None:
        resID = createResourceID(origintime, 'orig', authority_id=authority_id)
    elif isinstance(resID, str):
        resID = createResourceID(origintime, 'orig', authority_id=authority_id,
                                 hrstr=resID)

    origin = ope.Origin()
    origin.resource_id = resID
    origin.time = UTCDateTime(origintime)
    origin.creation_info = cinfo
    origin.latitude = latitude
    origin.longitude = longitude
    origin.depth = depth
    return origin


def createEvent(origintime, cinfo, etype, resID=None, authority_id=None):
    '''
    createEvent - funtion to create an ObsPy Event
    :param origintime: the events origintime
    :type origintime: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :param cinfo: An ObsPy :class: `~obspy.core.event.CreationInfo` object
        holding information on the creation of the returned object
    :type cinfo: :class: `~obspy.core.event.CreationInfo` object
    :param etype: Event type str object. converted via ObsPy to a valid event
        type string.
    :type etype: str
    :param resID: Resource identifier of the created event
    :type resID: :class: `~obspy.core.event.ResourceIdentifier` object
    :param authority_id: name of the institution carrying out the processing
    :type authority_id: str
    :return: An ObsPy :class: `~obspy.core.event.Event` object
    '''
    etype = ope.EventType(etype)
    if etype is None:
        etype = ope.EventType('earthquake')  # defaults to 'earthquake'
    if resID is None:
        resID = createResourceID(origintime, etype, authority_id)
    elif isinstance(resID, str):
        resID = createResourceID(origintime, etype, authority_id, resID)
    event = ope.Event(resource_id=resID)
    event.creation_info = cinfo
    event.event_type = etype
    return event


def createPick(origintime, picknum, picktime, eventnum, cinfo, phase, station,
               wfseedstr, authority_id):
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
    pickID = eventnum + '_' + station + '/{0:3d}'.format(picknum)
    pickresID = createResourceID(origintime, 'pick', authority_id, pickID)
    pick = ope.Pick()
    pick.resource_id = pickresID
    pick.time = picktime
    pick.creation_info = cinfo
    pick.phase_hint = phase
    pick.waveform_id = ope.ResourceIdentifier(id=wfseedstr, prefix='file:/')
    return pick


def createArrival(origintime, pickresID, eventnum, cinfo, phase, station,
                  authority_id, azimuth=None, dist=None):
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
    arriresID = createResourceID(origintime, 'arrival', authority_id, eventnum)
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


def createMagnitude(originID, origintime, cinfo, authority_id=None):
    '''
    createMagnitude - function to create an ObsPy Magnitude object
    :param originID:
    :param origintime:
    :param cinfo:
    :param authority_id:
    :return:
    '''
    magnresID = createResourceID(origintime, 'mag', authority_id)
    magnitude = ope.Magnitude()
    magnitude.resource_id = magnresID
    magnitude.creation_info = cinfo
    magnitude.origin_id = originID
    return magnitude


def createAmplitude(pickID, amp, unit, category, origintime, cinfo,
                    authority_id=None):
    amplresID = createResourceID(origintime, 'ampl', authority_id)
    amplitude = ope.Amplitude()
    amplitude.resource_id = amplresID
    amplitude.creation_info = cinfo
    amplitude.generic_amplitude = amp
    amplitude.unit = ope.AmplitudeUnit(unit)
    amplitude.type = ope.AmplitudeCategory(category)
    amplitude.pick_id = pickID
    return amplitude
