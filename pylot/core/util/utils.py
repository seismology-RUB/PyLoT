#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import pwd
import re
import hashlib
import numpy as np
from obspy.core import UTCDateTime
import obspy.core.event as ope

def createAmplitude(pickID, amp, unit, category, cinfo):
    '''

    :param pickID:
    :param amp:
    :param unit:
    :param category:
    :param cinfo:
    :return:
    '''
    amplitude = ope.Amplitude()
    amplitude.creation_info = cinfo
    amplitude.generic_amplitude = amp
    amplitude.unit = ope.AmplitudeUnit(unit)
    amplitude.type = ope.AmplitudeCategory(category)
    amplitude.pick_id = pickID
    return amplitude

def createArrival(pickresID, cinfo, phase, azimuth=None, dist=None):
    '''
    createArrival - function to create an Obspy Arrival

    :param pickresID: Resource identifier of the created pick
    :type pickresID: :class: `~obspy.core.event.ResourceIdentifier` object
    :param cinfo: An ObsPy :class: `~obspy.core.event.CreationInfo` object
    holding information on the creation of the returned object
    :type cinfo: :class: `~obspy.core.event.CreationInfo` object
    :param phase: name of the arrivals seismic phase
    :type phase: str
    :param azimuth: azimuth between source and receiver
    :type azimuth: float or int, optional
    :param dist: distance between source and receiver
    :type dist: float or int, optional
    :return: An ObsPy :class: `~obspy.core.event.Arrival` object
    '''
    arrival = ope.Arrival()
    arrival.creation_info = cinfo
    arrival.pick_id = pickresID
    arrival.phase = phase
    if azimuth is not None:
        arrival.azimuth = float(azimuth) if azimuth > -180 else azimuth + 360.
    else:
        arrival.azimuth = azimuth
    arrival.distance = dist
    return arrival

def createCreationInfo(agency_id=None, creation_time=None, author=None):
    '''

    :param agency_id:
    :param creation_time:
    :param author:
    :return:
    '''
    if author is None:
        author = getLogin()
    if creation_time is None:
        creation_time = UTCDateTime()
    return ope.CreationInfo(agency_id=agency_id, author=author,
                            creation_time=creation_time)

def createEvent(origintime, cinfo, originloc=None, etype=None, resID=None,
                authority_id=None):
    '''
    createEvent - funtion to create an ObsPy Event

    :param origintime: the events origintime
    :type origintime: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :param cinfo: An ObsPy :class: `~obspy.core.event.CreationInfo` object
        holding information on the creation of the returned object
    :type cinfo: :class: `~obspy.core.event.CreationInfo` object
    :param originloc: tuple containing the location of the origin
        (LAT, LON, DEP) affiliated with the event which is created
    :type originloc: tuple, list
    :param etype: Event type str object. converted via ObsPy to a valid event
        type string.
    :type etype: str
    :param resID: Resource identifier of the created event
    :type resID: :class: `~obspy.core.event.ResourceIdentifier` object, str
    :param authority_id: name of the institution carrying out the processing
    :type authority_id: str
    :return: An ObsPy :class: `~obspy.core.event.Event` object
    '''
    etype = ope.EventType(etype)
    if originloc is not None:
        o = createOrigin(origintime, cinfo,
                         originloc[0], originloc[1], originloc[2])
    else:
        o = None
    if etype is None:
        etype = ope.EventType('earthquake')  # defaults to 'earthquake'
    if not resID:
        resID = createResourceID(origintime, etype, authority_id)
    elif isinstance(resID, str):
        resID = createResourceID(origintime, etype, authority_id, resID)
    elif not isinstance(resID, ope.ResourceIdentifier):
        raise TypeError("unsupported type(resID) for resource identifier "
                        "generation: %s" % type(resID))
    event = ope.Event(resource_id=resID)
    event.creation_info = cinfo
    event.event_type = etype
    if o:
        event.origins = [o]
    return event

def createMagnitude(originID, cinfo):
    '''
    createMagnitude - function to create an ObsPy Magnitude object
    :param originID:
    :type originID:
    :param cinfo:
    :type cinfo:
    :return:
    '''
    magnitude = ope.Magnitude()
    magnitude.creation_info = cinfo
    magnitude.origin_id = originID
    return magnitude

def createOrigin(origintime, cinfo, latitude, longitude, depth):
    '''
    createOrigin - function to create an ObsPy Origin
    :param origintime: the origins time of occurence
    :type origintime: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :param cinfo:
    :type cinfo:
    :param latitude: latitude in decimal degree of the origins location
    :type latitude: float
    :param longitude: longitude in decimal degree of the origins location
    :type longitude: float
    :param depth: hypocentral depth of the origin
    :type depth: float
    :return: An ObsPy :class: `~obspy.core.event.Origin` object
    '''

    assert isinstance(origintime, UTCDateTime), "origintime has to be " \
                                                "a UTCDateTime object, but " \
                                                "actually is of type " \
                                                "'%s'" % type(origintime)

    origin = ope.Origin()
    origin.time = origintime
    origin.creation_info = cinfo
    origin.latitude = latitude
    origin.longitude = longitude
    origin.depth = depth
    return origin

def createPick(origintime, picknum, picktime, eventnum, cinfo, phase, station,
               wfseedstr, authority_id):
    '''
    createPick - function to create an ObsPy Pick

    :param origintime:
    :type origintime:
    :param picknum: number of the created pick
    :type picknum: int
    :param picktime:
    :type picktime:
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
    pickID = eventnum + '_' + station.strip() + '/{0:03d}'.format(picknum)
    pickresID = createResourceID(origintime, 'pick', authority_id, pickID)
    pick = ope.Pick()
    pick.resource_id = pickresID
    pick.time = picktime
    pick.creation_info = cinfo
    pick.phase_hint = phase
    pick.waveform_id = ope.ResourceIdentifier(id=wfseedstr, prefix='file:/')
    return pick

def createResourceID(timetohash, restype, authority_id=None, hrstr=None):
    '''

    :param timetohash:
    :type timetohash
    :param restype: type of the resource, e.g. 'orig', 'earthquake' ...
    :type restype: str
    :param authority_id: name of the institution carrying out the processing
    :type authority_id: str, optional
    :param hrstr:
    :type hrstr:
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

def demeanTrace(trace, window):
    """
    returns the DATA where each trace is demean by the average value within
    WINDOW
    :param trace: waveform trace object
    :type trace: `~obspy.core.stream.Trace`
    :param window:
    :type window: tuple
    :return: trace
    :rtype: `~obspy.core.stream.Trace`
    """
    trace.data -= trace.data[window].mean()
    return trace

def findComboBoxIndex(combo_box, val):
    """
    Function findComboBoxIndex takes a QComboBox object and a string and
    returns either 0 or the index throughout all QComboBox items.
    :param combo_box: Combo box object.
    :type combo_box: QComboBox
    :param val: Name of a combo box to search for.
    :type val:
    :return: index value of item with name val or 0
    """
    return combo_box.findText(val) if combo_box.findText(val) is not -1 else 0

def fnConstructor(s):
    '''

    :param s:
    :type s:
    :return:
    '''
    if type(s) is str:
        s = s.split(':')[-1]
    else:
        s = getHash(UTCDateTime())

    badchars = re.compile(r'[^A-Za-z0-9_. ]+|^\.|\.$|^ | $|^$')
    badsuffix = re.compile(r'(aux|com[1-9]|con|lpt[1-9]|prn)(\.|$)')

    fn = badchars.sub('_', s)

    if badsuffix.match(fn):
        fn = '_' + fn
    return fn

def getGlobalTimes(stream):
    '''

    :param stream:
    :type stream
    :return:
    '''
    min_start = UTCDateTime()
    max_end = None
    for trace in stream:
        if trace.stats.starttime < min_start:
            min_start = trace.stats.starttime
        if max_end is None or trace.stats.endtime > max_end:
            max_end = trace.stats.endtime
    return min_start, max_end

def getHash(time):
    '''
    :param time: time object for which a hash should be calculated
    :type time: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :return: str
    '''
    hg = hashlib.sha1()
    hg.update(time.strftime('%Y-%m-%d %H:%M:%S.%f'))
    return hg.hexdigest()

def getLogin():
    '''

    :return:
    '''
    return pwd.getpwuid(os.getuid())[0]

def getOwner(fn):
    '''

    :param fn:
    :type fn:
    :return:
    '''
    return pwd.getpwuid(os.stat(fn).st_uid).pw_name

def getPatternLine(fn, pattern):
    """
    Takes a file name and a pattern string to search for in the file and
    returns the first line which contains the pattern string otherwise None.

    :param fn: file name
    :type fn: str
    :param pattern: pattern string to search for
    :type pattern: str
    :return: the complete line containing pattern or None

    >>> getPatternLine('utils.py', 'python')
    '#!/usr/bin/env python\\n'
    >>> print(getPatternLine('version.py', 'palindrome'))
    None
    """
    fobj = open(fn, 'r')
    for line in fobj.readlines():
        if pattern in line:
            fobj.close()
            return line

    return None

def isSorted(iterable):
    '''

    :param iterable:
    :type iterable:
    :return:
    '''
    return sorted(iterable) == iterable

def prepTimeAxis(stime, trace):
    '''

    :param stime:
    :type stime:
    :param trace:
     :type trace:
    :return:
    '''
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

def scaleWFData(data, factor=None, components='all'):
    """
    produce scaled waveforms from given waveform data and a scaling factor,
    waveform may be selected by their components name
    :param data: waveform data to be scaled
    :type data: `~obspy.core.stream.Stream` object
    :param factor: scaling factor
    :type factor: float
    :param components: components labels for the traces in data to be scaled by
     the scaling factor (optional, default: 'all')
    :type components: tuple
    :return:  scaled waveform data
    :rtype: `~obspy.core.stream.Stream` object
    """
    if components is not 'all':
        for comp in components:
            if factor is None:
                max_val = np.max(np.abs(data.select(component=comp)[0].data))
                data.select(component=comp)[0].data /= 2 * max_val
            else:
                data.select(component=comp)[0].data /= 2 * factor
    else:
        for tr in data:
            if factor is None:
                max_val = float(np.max(np.abs(tr.data)))
                tr.data /= 2 * max_val
            else:
                tr.data /= 2 * factor

    return data

def runProgram(cmd, parameter=None):
    """
    run an external program specified by cmd with parameters input returning the
    stdout output

    :param cmd: name of the command to run
    :type cmd: str
    :param parameter: filename of parameter file  or parameter string
    :type parameter: str
    :return: stdout output
    :rtype: str
    """

    if parameter:
        cmd.strip()
        cmd += ' %s 2>&1' % parameter

    output = subprocess.check_output('{} | tee /dev/stderr'.format(cmd),
                                     shell = True)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
