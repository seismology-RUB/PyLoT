#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import scipy.io as sio
import obspy.core.event as ope
from obspy.core import UTCDateTime

from pylot.core.util import getOwner, createPick, createArrival, createEvent, \
    createOrigin, createMagnitude


def readPILOTEvent(phasfn=None, locfn=None, authority_id=None, **kwargs):
    """
    readPILOTEvent - function

    Reads Matlab PHASES and LOC files written by Matlab versions of PILOT and
    converts the data into an ObsPy Event object which is returned to the
    calling program.

    :rtype : ~obspy.core.event.Event
    :param eventID:
    :param authority:
    :param kwargs:
    :param phasfn: filename of the old PILOT Matlab PHASES file
    :param locfn: filename of the old PILOT Matlab LOC file
    :return event:  event object containing event and phase information
    """
    sdir = os.path.split(phasfn)[0]
    if phasfn is not None and os.path.isfile(phasfn):
        phases = sio.loadmat(phasfn)
        phasctime = UTCDateTime(os.path.getmtime(phasfn))
        phasauthor = getOwner(phasfn)
    else:
        phases = None
        phasctime = None
        phasauthor = None
    if locfn is not None and os.path.isfile(locfn):
        loc = sio.loadmat(locfn)
        locctime = UTCDateTime(os.path.getmtime(locfn))
        locauthor = getOwner(locfn)
    else:
        loc = None
        locctime = None
        locauthor = None
    pickcinfo = ope.CreationInfo(agency_id=authority_id,
                                 author=phasauthor,
                                 creation_time=phasctime)
    loccinfo = ope.CreationInfo(agency_id=authority_id,
                                author=locauthor,
                                creation_time=locctime)
    np = 0
    try:
        eventNum = loc['ID'][0]

        # retrieve eventID for the actual database
        idsplit = eventNum.split('.')

        # retrieve date information
        julday = int(idsplit[1])
        year = int(idsplit[2])
        hour = int(loc['hh'])
        minute = int(loc['mm'])
        second = int(loc['ss'])

        if UTCDateTime(year=year + 2000) < UTCDateTime.utcnow():
            year += 2000
        else:
            year += 1900

        eventDate = UTCDateTime(year=year, julday=julday, hour=hour,
                                minute=minute, second=second)

        stations = [stat for stat in phases['stat'][0:-1:3]]

        event = createEvent(eventDate, loccinfo, 'earthquake', eventNum,
                            authority_id)

        lat = float(loc['LAT'])
        lon = float(loc['LON'])
        dep = float(loc['DEP'])

        origin = createOrigin(eventDate, loccinfo, lat, lon, dep, eventNum)
        for n, pick in enumerate(phases['Ptime']):
            kwargs = {'year': int(pick[0]),
                      'month': int(pick[1]),
                      'day': int(pick[2]),
                      'hour': int(pick[3]),
                      'minute': int(pick[4]),
                      'second': int(str(pick[5]).split('.')[0]),
                      'microsecond': int(str(pick[5]).split('.')[1][0:6])}
            spick = phases['Stime'][n]
            if spick[0] > 0:
                skwargs = {'year': int(spick[0]),
                           'month': int(spick[1]),
                           'day': int(spick[2]),
                           'hour': int(spick[3]),
                           'minute': int(spick[4]),
                           'second': int(str(spick[5]).split('.')[0]),
                           'microsecond': int(str(spick[5]).split('.')[1][0:6])}
                spicktime = UTCDateTime(**skwargs)
            else:
                spicktime = None
            ppicktime = UTCDateTime(**kwargs)

            for picktime, phase in [(ppicktime, 'P'), (spicktime, 'S')]:
                if phase == 'P':
                    wffn = os.path.join([sdir, '{0}*{1}*'.format(stations[n],
                                                                 'z')])
                else:
                    wffn = os.path.join([sdir, '{0}*{1}*'.format(stations[n],
                                                                 '[ne]')])
                pick = createPick(eventDate, np, picktime, eventNum, pickcinfo,
                                  phase,
                                  stat[n], wffn, authority_id)
                event.picks.append(pick)
                pickID = pick.get('resource_id')
                arrival = createArrival(pickID, eventNum, pickcinfo, phase,
                                        stat[n], authority_id)
                origin.arrivals.append(arrival)
                np += 1

        magnitude = createMagnitude(origin.get('resource_id'), eventDate,
                                    loccinfo,
                                    authority_id)
        magnitude.mag = float(loc['Mnet'])
        magnitude.magnitude_type = 'Ml'

        event.picks.append(pick)
        event.origins.append(origin)
        event.magnitudes.append(magnitude)

    except AttributeError, e:
        raise AttributeError('{0} - Matlab LOC files {1} and {2} contains \
                              insufficient data!'.format(e, phasfn, locfn))



