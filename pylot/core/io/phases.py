#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import warnings
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
                print(wffn)
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


def picks_from_pilot(fn):
    picks = dict()
    phases_pilot = sio.loadmat(fn)
    stations = stations_from_pilot(phases_pilot['stat'])
    for n, station in enumerate(stations):
        phases = dict()
        for onset_name in 'PS':
            onset_label = '{0}time'.format(onset_name)
            pick = phases_pilot[onset_label][n]
            if not pick[0]:
                continue
            pick = convert_pilot_times(pick)
            phases[onset_name] = dict(mpp=pick)
        picks[station] = phases

    return picks


def stations_from_pilot(stat_array):
    stations = list()
    cur_stat = None
    for stat in stat_array:
        if stat == cur_stat:
            continue
        cur_stat = stat
        stations.append(stat.strip())

    return stations


def convert_pilot_times(time_array):
    times = [int(time) for time in time_array]
    microseconds = int((time_array[-1] - times[-1]) * 1e6)
    times.append(microseconds)
    return UTCDateTime(*times)


def picksdict_from_obs(fn):
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


def picks_to_dict(evt):
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

def picks_from_dict(picks):
    firstonset = None
    for station, onsets in picks.items():
        print('Reading picks on station %s' % station)
        for label, phase in onsets.items():
            if not isinstance(phase, dict):
                continue
            onset = phase['mpp']
            epp = phase['epp']
            lpp = phase['lpp']
            error = phase['spe']
            try:
                picker = phase['picker']
            except KeyError as e:
                warnings.warn(str(e), Warning)
                picker = 'Unknown'
            pick = ope.Pick()
            pick.time = onset
            pick.time_errors.lower_uncertainty = onset - epp
            pick.time_errors.upper_uncertainty = lpp - onset
            pick.time_errors.uncertainty = error
            pick.phase_hint = label
            pick.method_id = ope.ResourceIdentifier(id=picker)
            pick.waveform_id = ope.WaveformStreamID(station_code=station)
            try:
                polarity = phase['fm']
                if polarity == 'U' or '+':
                    pick.polarity = 'positive'
                elif polarity == 'D' or '-':
                    pick.polarity = 'negative'
                else:
                    pick.polarity = 'undecidable'
            except KeyError as e:
                print('No polarity information found for %s' % phase)
            if firstonset is None or firstonset > onset:
                firstonset = onset


def reassess_pilot_event(root_dir, event_id):
    from obspy import read
    from pylot.core.util.defaults import AUTOMATIC_DEFAULTS
    from pylot.core.io.inputs import AutoPickParameter
    from pylot.core.pick.utils import earllatepicker

    default = AutoPickParameter(AUTOMATIC_DEFAULTS)

    search_base = os.path.join(root_dir, event_id)
    phases_file = glob.glob(os.path.join(search_base, 'PHASES.mat'))
    picks_dict = picks_from_pilot(phases_file)
    for station in picks_dict.keys():
        fn_pattern = os.path.join(search_base, '{0}*'.format(station))
        try:
            st = read(fn_pattern)
        except TypeError as e:
            print(e.message)
            st = read(fn_pattern, format='GSE2')
        if not st:
            raise RuntimeError('no waveform data found for station {station}'.format(station=station))
        for phase in picks_dict[station].keys():
            try:
                mpp = picks_dict[station][phase]['mpp']
            except KeyError as e:
                print(e.message, station)
                continue
            epp, lpp, spe = earllatepicker(st,
                                           default.get('nfac{0}'.format(phase)),
                                           default.get('tsnrz' if phase == 'P' else 'tsnrh'),
                                           mpp
                                           )
            picks_dict[station][phase] = dict(epp=epp, mpp=mpp, lpp=lpp, spe=spe)
    # create Event object for export
    evt = ope.Event(resource_id=event_id)
    evt.picks = picks_from_dict(picks_dict)
    # write phase information to file
    evt.write('{0}.xml'.format(event_id), format='QUAKEML')


def writephases(arrivals, fformat, filename):
    '''
    Function of methods to write phases to the following standard file
    formats used for locating earthquakes:

    HYPO71, NLLoc, VELEST, HYPOSAT, and hypoDD

    :param: arrivals
    :type: dictionary containing all phase information including
           station ID, phase, first motion, weight (uncertainty),
           ....

    :param: fformat
    :type:  string, chosen file format (location routine),
            choose between NLLoc, HYPO71, HYPOSAT, VELEST,
            HYPOINVERSE, and hypoDD

    :param: filename, full path and name of phase file
    :type: string
    '''

    if fformat == 'NLLoc':
        print ("Writing phases to %s for NLLoc" % filename)
        fid = open("%s" % filename, 'w')
        # write header
        fid.write('# EQEVENT:  Label: EQ001  Loc:  X 0.00  Y 0.00  Z 10.00  OT 0.00 \n')
        for key in arrivals:
            # P onsets
            if arrivals[key]['P']:
                fm = arrivals[key]['P']['fm']
                if fm == None:
                    fm = '?'
                onset = arrivals[key]['P']['mpp']
                year = onset.year
                month = onset.month
                day = onset.day
                hh = onset.hour
                mm = onset.minute
                ss = onset.second
                ms = onset.microsecond
                ss_ms = ss + ms / 1000000.0
                if arrivals[key]['P']['weight'] < 4:
                    pweight = 1  # use pick
                else:
                    pweight = 0  # do not use pick
                fid.write('%s ? ? ? P   %s %d%02d%02d %02d%02d %7.4f GAU 0 0 0 0 %d \n' % (key,
                                                                                           fm,
                                                                                           year,
                                                                                           month,
                                                                                           day,
                                                                                           hh,
                                                                                           mm,
                                                                                           ss_ms,
                                                                                           pweight))
            # S onsets
            if arrivals[key]['S']:
                fm = '?'
                onset = arrivals[key]['S']['mpp']
                year = onset.year
                month = onset.month
                day = onset.day
                hh = onset.hour
                mm = onset.minute
                ss = onset.second
                ms = onset.microsecond
                ss_ms = ss + ms / 1000000.0
                if arrivals[key]['S']['weight'] < 4:
                    sweight = 1  # use pick
                else:
                    sweight = 0  # do not use pick
                fid.write('%s ? ? ? S   %s %d%02d%02d %02d%02d %7.4f GAU 0 0 0 0 %d \n' % (key,
                                                                                           fm,
                                                                                           year,
                                                                                           month,
                                                                                           day,
                                                                                           hh,
                                                                                           mm,
                                                                                           ss_ms,
                                                                                           sweight))

        fid.close()

    elif fformat == 'HYPO71':
        print ("Writing phases to %s for HYPO71" % filename)
        fid = open("%s" % filename, 'w')
        # write header
        fid.write('                                                              EQ001\n')
        for key in arrivals:
            if arrivals[key]['P']['weight'] < 4:
                Ponset = arrivals[key]['P']['mpp']
                Sonset = arrivals[key]['S']['mpp']
                pweight = arrivals[key]['P']['weight']
                sweight = arrivals[key]['S']['weight']
                fm = arrivals[key]['P']['fm']
                if fm is None:
                    fm = '-'
                Ao = arrivals[key]['S']['Ao']
                if Ao is None:
                    Ao = ''
                else:
                    Ao = str('%7.2f' % Ao)
                year = Ponset.year
                if year >= 2000:
                    year = year - 2000
                else:
                    year = year - 1900
                month = Ponset.month
                day = Ponset.day
                hh = Ponset.hour
                mm = Ponset.minute
                ss = Ponset.second
                ms = Ponset.microsecond
                ss_ms = ss + ms / 1000000.0
                if pweight < 2:
                    pstr = 'I'
                elif pweight >= 2:
                    pstr = 'E'
                if arrivals[key]['S']['weight'] < 4:
                    Sss = Sonset.second
                    Sms = Sonset.microsecond
                    Sss_ms = Sss + Sms / 1000000.0
                    Sss_ms = str('%5.02f' % Sss_ms)
                    if sweight < 2:
                        sstr = 'I'
                    elif sweight >= 2:
                        sstr = 'E'
                    fid.write('%s%sP%s%d %02d%02d%02d%02d%02d%5.2f       %s%sS %d   %s\n' % (key,
                                                                                             pstr,
                                                                                             fm,
                                                                                             pweight,
                                                                                             year,
                                                                                             month,
                                                                                             day,
                                                                                             hh,
                                                                                             mm,
                                                                                             ss_ms,
                                                                                             Sss_ms,
                                                                                             sstr,
                                                                                             sweight,
                                                                                             Ao))
                else:
                    fid.write('%s%sP%s%d %02d%02d%02d%02d%02d%5.2f                  %s\n' % (key,
                                                                                             pstr,
                                                                                             fm,
                                                                                             pweight,
                                                                                             year,
                                                                                             month,
                                                                                             day,
                                                                                             hh,
                                                                                             mm,
                                                                                             ss_ms,
                                                                                             Ao))

        fid.close()
