#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import scipy.io as sio
import obspy.core.event as ope
from obspy.core import UTCDateTime

from pylot.core.util.utils import getOwner, createPick, createArrival, \
    createEvent, createOrigin, createMagnitude


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

        if year + 2000 < UTCDateTime.utcnow().year:
            year += 2000
        else:
            year += 1900

        eventDate = UTCDateTime(year=year, julday=julday, hour=hour,
                                minute=minute, second=second)

        stations = [stat for stat in phases['stat'][0:-1:3]]

        event = createEvent(eventDate, loccinfo, etype='earthquake', resID=eventNum,
                            authority_id=authority_id)

        lat = float(loc['LAT'])
        lon = float(loc['LON'])
        dep = float(loc['DEP'])

        origin = createOrigin(eventDate, loccinfo, lat, lon, dep)
        for n, pick in enumerate(phases['Ptime']):
            if pick[0] > 0:
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
                if picktime is not None:
                    if phase == 'P':
                        wffn = os.path.join(sdir, '{0}*{1}*'.format(
                            stations[n].strip(), 'z'))
                    else:
                        wffn = os.path.join(sdir, '{0}*{1}*'.format(
                            stations[n].strip(), '[ne]'))
                print wffn
                pick = createPick(eventDate, np, picktime, eventNum, pickcinfo,
                                  phase, stations[n], wffn, authority_id)
                event.picks.append(pick)
                pickID = pick.get('id')
                arrival = createArrival(pickID, pickcinfo, phase)
                origin.arrivals.append(arrival)
                np += 1

        magnitude = createMagnitude(origin.get('id'), loccinfo)
        magnitude.mag = float(loc['Mnet'])
        magnitude.magnitude_type = 'Ml'

        event.picks.append(pick)
        event.origins.append(origin)
        event.magnitudes.append(magnitude)
        return event

    except AttributeError as e:
        raise AttributeError('{0} - Matlab LOC files {1} and {2} contains \
                              insufficient data!'.format(e, phasfn, locfn))

def picks_from_obs(fn):
    picks = dict()
    station_name = str()
    for line in open(fn, 'r'):
        if line.startswith('#'):
            continue
        else:
            phase_line = line.split()
            if not station_name == phase_line[0]:
                phase = dict()
            station_name = phase_line[0]
            phase_name = phase_line[4].upper()
            pick = UTCDateTime(phase_line[6] + phase_line[7] + phase_line[8])
            phase[phase_name] = dict(mpp=pick, fm=phase_line[5])
            picks[station_name] = phase
    return picks


def picks_from_evt(evt):
    '''
    Takes an Event object and return the pick dictionary commonly used within
    PyLoT
    :param evt: Event object contain all available information
    :type evt: `~obspy.core.event.Event`
    :return: pick dictionary
    '''
    picks = {}
    for pick in evt.picks:
        phase = {}
        station = pick.waveform_id.station_code
        try:
            onsets = picks[station]
        except KeyError as e:
            print(e)
            onsets = {}
        mpp = pick.time
        lpp = mpp + pick.time_errors.upper_uncertainty
        epp = mpp - pick.time_errors.lower_uncertainty
        spe = pick.time_errors.uncertainty
        phase['mpp'] = mpp
        phase['epp'] = epp
        phase['lpp'] = lpp
        phase['spe'] = spe
        try:
            picker = str(pick.method_id)
            if picker.startswith('smi:local/'):
                picker = picker.split('smi:local/')[1]
            phase['picker'] = picker
        except IndexError:
            pass

        onsets[pick.phase_hint] = phase.copy()
        picks[station] = onsets.copy()
    return picks
