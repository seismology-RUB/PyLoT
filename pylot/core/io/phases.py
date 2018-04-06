#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import obspy.core.event as ope
import scipy.io as sio
from obspy.core import UTCDateTime
from obspy.core.event import read_events
from obspy.core.util import AttribDict
from pylot.core.io.inputs import PylotParameter
from pylot.core.io.location import create_event, \
    create_magnitude
from pylot.core.pick.utils import select_for_phase
from pylot.core.util.utils import getOwner, full_range, four_digits, transformFilteroptions2String, \
    transformFilterString4Export, backtransformFilterString


def add_amplitudes(event, amplitudes):
    amplitude_list = []
    for pick in event.picks:
        try:
            a0 = amplitudes[pick.waveform_id.station_code]
            amplitude = ope.Amplitude(generic_amplitude=a0 * 1e-3)
            amplitude.unit = 'm'
            amplitude.category = 'point'
            amplitude.waveform_id = pick.waveform_id
            amplitude.magnitude_hint = 'ML'
            amplitude.pick_id = pick.resource_id
            amplitude.type = 'AML'
            amplitude_list.append(amplitude)
        except KeyError:
            continue
    event.amplitudes = amplitude_list
    return event


def readPILOTEvent(phasfn=None, locfn=None, authority_id='RUB', **kwargs):
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

    eventNum = str(loc['ID'][0])

    # retrieve eventID for the actual database
    idsplit = eventNum.split('.')

    # retrieve date information
    julday = int(idsplit[1])
    year = int(idsplit[2])
    hour = int(loc['hh'])
    minute = int(loc['mm'])
    second = int(loc['ss'])

    year = four_digits(year)

    eventDate = UTCDateTime(year=year, julday=julday, hour=hour,
                            minute=minute, second=second)

    stations = [stat for stat in phases['stat'][0:-1:3]]

    lat = float(loc['LAT'])
    lon = float(loc['LON'])
    dep = float(loc['DEP'])

    event = create_event(eventDate, loccinfo, originloc=(lat, lon, dep),
                         etype='earthquake', resID=eventNum,
                         authority_id=authority_id)

    picks = picksdict_from_pilot(phasfn)

    event.picks = picks_from_picksdict(picks, creation_info=pickcinfo)

    if event.origins:
        origin = event.origins[0]
        magnitude = create_magnitude(origin.get('id'), loccinfo)
        magnitude.mag = float(loc['Mnet'])
        magnitude.magnitude_type = 'Ml'
        event.magnitudes.append(magnitude)
    return event


def picksdict_from_pilot(fn):
    """
    Create pick dictionary from matlab file
    :param fn: matlab file
    :type fn:
    :return: pick dictionary
    :rtype: dict
    """
    from pylot.core.util.defaults import TIMEERROR_DEFAULTS
    picks = dict()
    phases_pilot = sio.loadmat(fn)
    stations = stations_from_pilot(phases_pilot['stat'])
    params = PylotParameter(TIMEERROR_DEFAULTS)
    timeerrors = dict(P=params.get('timeerrorsP'),
                      S=params.get('timeerrorsS'))
    for n, station in enumerate(stations):
        phases = dict()
        for onset_name in 'PS':
            onset_label = '{0}time'.format(onset_name)
            pick = phases_pilot[onset_label][n]
            if not pick[0]:
                continue
            pick = convert_pilot_times(pick)
            uncertainty_label = '{0}weight'.format(onset_name.lower())
            ierror = phases_pilot[uncertainty_label][0, n]
            try:
                spe = timeerrors[onset_name][ierror]
            except IndexError as e:
                print(e.message + '\ntake two times the largest default error value')
                spe = timeerrors[onset_name][-1] * 2
            phases[onset_name] = dict(mpp=pick, spe=spe, weight=ierror)
        picks[station] = phases

    return picks


def stations_from_pilot(stat_array):
    """
    Create stations list from pilot station array
    :param stat_array:
    :type stat_array:
    :return:
    :rtype: list
    """
    stations = list()
    cur_stat = None
    for stat in stat_array:
        stat = stat.strip()
        if stat == cur_stat:
            continue
        cur_stat = stat
        if stat not in stations:
            stations.append(stat)
        else:
            warnings.warn('station {0} listed at least twice, might corrupt '
                          'phase times', RuntimeWarning)

    return stations


def convert_pilot_times(time_array):
    """
    Convert pilot times to UTCDateTimes
    :param time_array: pilot times
    :type time_array:
    :return:
    :rtype:
    """
    times = [int(time) for time in time_array]
    microseconds = int((time_array[-1] - times[-1]) * 1e6)
    times.append(microseconds)
    return UTCDateTime(*times)


def picksdict_from_obs(fn):
    """
    create pick dictionary from obs file
    :param fn: filename
    :type fn:
    :return:
    :rtype:
    """
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


def picksdict_from_picks(evt):
    """
    Takes an Event object and return the pick dictionary commonly used within
    PyLoT
    :param evt: Event object contain all available information
    :type evt: `~obspy.core.event.Event`
    :return: pick dictionary (auto and manual)
    """
    picksdict = {
        'manual': {},
        'auto': {}
    }
    for pick in evt.picks:
        phase = {}
        station = pick.waveform_id.station_code
        channel = pick.waveform_id.channel_code
        network = pick.waveform_id.network_code
        mpp = pick.time
        spe = pick.time_errors.uncertainty
        if pick.filter_id:
            filter_id = backtransformFilterString(str(pick.filter_id.id))
        else:
            filter_id = None
        try:
            picker = str(pick.method_id)
            if picker.startswith('smi:local/'):
                picker = picker.split('smi:local/')[1]
        except IndexError:
            picker = 'manual' # MP MP TODO maybe improve statement
        try:
            onsets = picksdict[picker][station]
        except KeyError as e:
            # print(e)
            onsets = {}
        try:
            lpp = mpp + pick.time_errors.upper_uncertainty
            epp = mpp - pick.time_errors.lower_uncertainty
        except TypeError as e:
            if not spe:
                msg = 'No uncertainties found for pick: {}. Uncertainty set to 0'.format(pick)
                lpp = mpp
                epp = mpp
            else:
                msg = str(e) + ',\n falling back to symmetric uncertainties'
                lpp = mpp + spe
                epp = mpp - spe
            warnings.warn(msg)
        phase['mpp'] = mpp
        phase['epp'] = epp
        phase['lpp'] = lpp
        phase['spe'] = spe
        phase['channel'] = channel
        phase['network'] = network
        phase['picker'] = picker
        phase['filter_id'] = filter_id if filter_id is not None else ''

        onsets[pick.phase_hint] = phase.copy()
        picksdict[picker][station] = onsets.copy()
    return picksdict


def picks_from_picksdict(picks, creation_info=None):
    """
    Create a list of picks out of a pick dictionary
    :param picks: pick dictionary
    :type picks: dict
    :param creation_info: obspy creation information to apply to picks
    :type creation_info:
    :param creation_info: obspy creation information to apply to picks
    :return: list of picks
    :rtype: list
    """
    picks_list = list()
    for station, onsets in picks.items():
        for label, phase in onsets.items():
            if not isinstance(phase, dict) and not isinstance(phase, AttribDict):
                continue
            onset = phase['mpp']
            try:
                ccode = phase['channel']
                ncode = phase['network']
            except:
                continue
            pick = ope.Pick()
            if creation_info:
                pick.creation_info = creation_info
            pick.time = onset
            error = phase['spe']
            pick.time_errors.uncertainty = error
            try:
                epp = phase['epp']
                lpp = phase['lpp']
                pick.time_errors.lower_uncertainty = onset - epp
                pick.time_errors.upper_uncertainty = lpp - onset
            except (KeyError, TypeError) as e:
                warnings.warn(str(e), RuntimeWarning)
            try:
                picker = phase['picker']
            except KeyError as e:
                warnings.warn(e.message, RuntimeWarning)
                picker = 'Unknown'
            pick.phase_hint = label
            pick.method_id = ope.ResourceIdentifier(id=picker)
            pick.waveform_id = ope.WaveformStreamID(station_code=station,
                                                    channel_code=ccode,
                                                    network_code=ncode)
            try:
                filter_id = phase['filteroptions']
                filter_id = transformFilterString4Export(filter_id)
            except KeyError as e:
                warnings.warn(e.message, RuntimeWarning)
                filter_id = ''
            pick.filter_id = filter_id
            try:
                polarity = phase['fm']
                if polarity == 'U' or '+':
                    pick.polarity = 'positive'
                elif polarity == 'D' or '-':
                    pick.polarity = 'negative'
                else:
                    pick.polarity = 'undecidable'
            except KeyError as e:
                if 'fm' in str(e):  # no polarity information found for this phase
                    pass
                else:
                    raise e
            picks_list.append(pick)
    return picks_list

def reassess_pilot_db(root_dir, db_dir, out_dir=None, fn_param=None, verbosity=0):
    import glob

    db_root = os.path.join(root_dir, db_dir)
    evt_list = glob.glob1(db_root, 'e????.???.??')

    for evt in evt_list:
        if verbosity > 0:
            print('Reassessing event {0}'.format(evt))
        reassess_pilot_event(root_dir, db_dir, evt, out_dir, fn_param, verbosity)


def reassess_pilot_event(root_dir, db_dir, event_id, out_dir=None, fn_param=None, verbosity=0):
    from obspy import read

    from pylot.core.io.inputs import PylotParameter
    from pylot.core.pick.utils import earllatepicker

    if fn_param is None:
        fn_param = defaults.AUTOMATIC_DEFAULTS

    default = PylotParameter(fn_param, verbosity)

    search_base = os.path.join(root_dir, db_dir, event_id)
    phases_file = glob.glob(os.path.join(search_base, 'PHASES.mat'))
    if not phases_file:
        return
    if verbosity > 1:
        print('Opening PILOT phases file: {fn}'.format(fn=phases_file[0]))
    picks_dict = picksdict_from_pilot(phases_file[0])
    if verbosity > 0:
        print('Dictionary read from PHASES.mat:\n{0}'.format(picks_dict))
    datacheck = list()
    info = None
    for station in picks_dict.keys():
        fn_pattern = os.path.join(search_base, '{0}*'.format(station))
        try:
            st = read(fn_pattern)
        except TypeError as e:
            if 'Unknown format for file' in e.message:
                try:
                    st = read(fn_pattern, format='GSE2')
                except ValueError as e:
                    if e.message == 'second must be in 0..59':
                        info = 'A known Error was raised. Please find the list of corrupted files and double-check these files.'
                        datacheck.append(fn_pattern + ' (time info)\n')
                        continue
                    else:
                        raise ValueError(e.message)
                except Exception as e:
                    if 'No file matching file pattern:' in e.message:
                        if verbosity > 0:
                            warnings.warn('no waveform data found for station {station}'.format(station=station),
                                          RuntimeWarning)
                        datacheck.append(fn_pattern + ' (no data)\n')
                        continue
                    else:
                        raise e
            else:
                raise e
        for phase in picks_dict[station].keys():
            try:
                mpp = picks_dict[station][phase]['mpp']
            except KeyError as e:
                print(e.message, station)
                continue
            sel_st = select_for_phase(st, phase)
            if not sel_st:
                msg = 'no waveform data found for station {station}'.format(station=station)
                warnings.warn(msg, RuntimeWarning)
                continue
            stime, etime = full_range(sel_st)
            rel_pick = mpp - stime
            epp, lpp, spe = earllatepicker(sel_st,
                                           default.get('nfac{0}'.format(phase)),
                                           default.get('tsnrz' if phase == 'P' else 'tsnrh'),
                                           Pick1=rel_pick,
                                           iplot=0,
                                           verbosity=0)
            if epp is None or lpp is None:
                continue
            epp = stime + epp
            lpp = stime + lpp
            min_diff = 3 * st[0].stats.delta
            if lpp - mpp < min_diff:
                lpp = mpp + min_diff
            if mpp - epp < min_diff:
                epp = mpp - min_diff
            picks_dict[station][phase] = dict(epp=epp, mpp=mpp, lpp=lpp, spe=spe)
    if datacheck:
        if info:
            if verbosity > 0:
                print(info + ': {0}'.format(search_base))
        fncheck = open(os.path.join(search_base, 'datacheck_list'), 'w')
        fncheck.writelines(datacheck)
        fncheck.close()
        del datacheck
    # create Event object for export
    evt = ope.Event(resource_id=event_id)
    evt.picks = picks_from_picksdict(picks_dict)
    # write phase information to file
    if not out_dir:
        fnout_prefix = os.path.join(root_dir, db_dir, event_id, 'PyLoT_{0}.'.format(event_id))
    else:
        out_dir = os.path.join(out_dir, db_dir)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        fnout_prefix = os.path.join(out_dir, 'PyLoT_{0}.'.format(event_id))
    evt.write(fnout_prefix + 'xml', format='QUAKEML')
    # evt.write(fnout_prefix + 'cnv', format='VELEST')


def writephases(arrivals, fformat, filename, parameter=None, eventinfo=None):
    """
    Function of methods to write phases to the following standard file
    formats used for locating earthquakes:

    HYPO71, NLLoc, VELEST, HYPOSAT, and hypoDD

    :param arrivals:dictionary containing all phase information including
     station ID, phase, first motion, weight (uncertainty), ...
    :type arrivals: dict

    :param fformat: chosen file format (location routine),
    choose between NLLoc, HYPO71, HYPOSAT, VELEST,
    HYPOINVERSE, and hypoDD
    :type fformat: str

    :param filename: full path and name of phase file
    :type filename:  string

    :param parameter: all input information
    :type parameter:  object

    :param eventinfo: optional, needed for VELEST-cnv file
            and FOCMEC- and HASH-input files 
    :type eventinfo: `obspy.core.event.Event` object
    """

    if fformat == 'NLLoc':
        print("Writing phases to %s for NLLoc" % filename)
        fid = open("%s" % filename, 'w')
        # write header
        fid.write('# EQEVENT: %s Label: EQ%s  Loc:  X 0.00  Y 0.00  Z 10.00  OT 0.00 \n' %
                  (parameter.get('database'), parameter.get('eventID')))
        for key in arrivals:
            # P onsets
            if arrivals[key].has_key('P'):
                try:
                    fm = arrivals[key]['P']['fm']
                except KeyError as e:
                    print(e)
                    fm = None
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
                pweight = 1  # use pick
                try:
                    if arrivals[key]['P']['weight'] >= 4:
                        pweight = 0  # do not use pick
                except KeyError as e:
                    print(e.message + '; no weight set during processing')
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
            if arrivals[key].has_key('S') and arrivals[key]['S']:
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
                sweight = 1  # use pick
                try:
                    if arrivals[key]['S']['weight'] >= 4:
                        sweight = 0  # do not use pick
                except KeyError as e:
                    print(str(e) + '; no weight set during processing')
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
        print("Writing phases to %s for HYPO71" % filename)
        fid = open("%s" % filename, 'w')
        # write header
        fid.write('                                                                %s\n' %
                  parameter.get('eventID'))
        for key in arrivals:
            if arrivals[key]['P']['weight'] < 4:
                stat = key
                if len(stat) > 4:  # HYPO71 handles only 4-string station IDs
                    stat = stat[1:5]
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
                    fid.write('%-4s%sP%s%d %02d%02d%02d%02d%02d%5.2f       %s%sS %d   %s\n' % (stat,
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
                    fid.write('%-4s%sP%s%d %02d%02d%02d%02d%02d%5.2f                  %s\n' % (stat,
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

    elif fformat == 'HYPOSAT':
        print("Writing phases to %s for HYPOSAT" % filename)
        fid = open("%s" % filename, 'w')
        # write header
        fid.write('%s, event %s \n' % (parameter.get('database'), parameter.get('eventID')))
        for key in arrivals:
            # P onsets
            if arrivals[key].has_key('P'):
                if arrivals[key]['P']['weight'] < 4:
                    Ponset = arrivals[key]['P']['mpp']
                    pyear = Ponset.year
                    pmonth = Ponset.month
                    pday = Ponset.day
                    phh = Ponset.hour
                    pmm = Ponset.minute
                    pss = Ponset.second
                    pms = Ponset.microsecond
                    Pss = pss + pms / 1000000.0
                    # use symmetrized picking error as std
                    # (read the HYPOSAT manual)
                    pstd = arrivals[key]['P']['spe']
                    fid.write('%-5s P1       %4.0f %02d %02d %02d %02d %05.02f   %5.3f -999.   0.00 -999.  0.00\n'
                              % (key, pyear, pmonth, pday, phh, pmm, Pss, pstd))
            # S onsets
            if arrivals[key].has_key('S') and arrivals[key]['S']:
                if arrivals[key]['S']['weight'] < 4:
                    Sonset = arrivals[key]['S']['mpp']
                    syear = Sonset.year
                    smonth = Sonset.month
                    sday = Sonset.day
                    shh = Sonset.hour
                    smm = Sonset.minute
                    sss = Sonset.second
                    sms = Sonset.microsecond
                    Sss = sss + sms / 1000000.0
                    sstd = arrivals[key]['S']['spe']
                    fid.write('%-5s S1       %4.0f %02d %02d %02d %02d %05.02f   %5.3f -999.   0.00 -999.  0.00\n'
                              % (key, syear, smonth, sday, shh, smm, Sss, sstd))
        fid.close()

    elif fformat == 'VELEST':
        print("Writing phases to %s for VELEST" % filename)
        fid = open("%s" % filename, 'w')
        # get informations needed in cnv-file
        # check, whether latitude is N or S and longitude is E or W
        try:
            eventsource = eventinfo.origins[0]
        except:
            print("No source origin calculated yet, thus no cnv-file creation possible!")
            return
        if eventsource['latitude'] < 0:
            cns = 'S'
        else:
            cns = 'N'
        if eventsource['longitude'] < 0:
            cew = 'W'
        else:
            cew = 'E'
        # get last two integers of origin year
        stime = eventsource['time']
        if stime.year - 2000 >= 0:
            syear = stime.year - 2000
        else:
            syear = stime.year - 1900
        ifx = 0  # default value, see VELEST manual, pp. 22-23
        # write header
        fid.write('%s%02d%02d %02d%02d %05.2f %7.4f%c %8.4f%c %7.2f %6.2f     %02.0f  0.0 0.03  1.0  1.0\n' % (
            syear, stime.month, stime.day, stime.hour, stime.minute, stime.second, eventsource['latitude'],
            cns, eventsource['longitude'], cew, eventsource['depth'], eventinfo.magnitudes[0]['mag'], ifx))
        n = 0
        for key in arrivals:
            # P onsets
            if arrivals[key].has_key('P'):
                if arrivals[key]['P']['weight'] < 4:
                    n += 1
                    stat = key
                    if len(stat) > 4:  # VELEST handles only 4-string station IDs
                        stat = stat[1:5]
                    Ponset = arrivals[key]['P']['mpp']
                    Pweight = arrivals[key]['P']['weight']
                    Prt = Ponset - stime  # onset time relative to source time
                    if n % 6 is not 0:
                        fid.write('%-4sP%d%6.2f' % (stat, Pweight, Prt))
                    else:
                        fid.write('%-4sP%d%6.2f\n' % (stat, Pweight, Prt))
                        # S onsets
            if arrivals[key].has_key('S'):
                if arrivals[key]['S']['weight'] < 4:
                    n += 1
                    stat = key
                    if len(stat) > 4:  # VELEST handles only 4-string station IDs
                        stat = stat[1:5]
                    Sonset = arrivals[key]['S']['mpp']
                    Sweight = arrivals[key]['S']['weight']
                    Srt = Ponset - stime  # onset time relative to source time
                    if n % 6 is not 0:
                        fid.write('%-4sS%d%6.2f' % (stat, Sweight, Srt))
                    else:
                        fid.write('%-4sS%d%6.2f\n' % (stat, Sweight, Srt))
        fid.close()

    elif fformat == 'hypoDD':
        print("Writing phases to %s for hypoDD" % filename)
        fid = open("%s" % filename, 'w')
        # get event information needed for hypoDD-phase file
        eventsource = eventinfo.origins[0]
        stime = eventsource['time']
        event = parameter.get('eventID')
        hddID = event.split('.')[0][1:5]
        # write header
        fid.write('# %d  %d %d %d %d %5.2f %7.4f +%6.4f %7.4f %4.2f 0.1 0.5 %4.2f      %s\n' % (
            stime.year, stime.month, stime.day, stime.hour, stime.minute, stime.second,
            eventsource['latitude'], eventsource['longitude'], eventsource['depth'] / 1000,
            eventinfo.magnitudes[0]['mag'], eventsource['quality']['standard_error'], hddID))
        for key in arrivals:
            if arrivals[key].has_key('P'):
                # P onsets
                if arrivals[key]['P']['weight'] < 4:
                    Ponset = arrivals[key]['P']['mpp']
                    Prt = Ponset - stime  # onset time relative to source time
                    fid.write('%s    %6.3f  1  P\n' % (key, Prt))
                    # S onsets
                if arrivals[key]['S']['weight'] < 4:
                    Sonset = arrivals[key]['S']['mpp']
                    Srt = Sonset - stime  # onset time relative to source time
                    fid.write('%-5s    %6.3f  1  S\n' % (key, Srt))

        fid.close()

    elif fformat == 'FOCMEC':
        print("Writing phases to %s for FOCMEC" % filename)
        fid = open("%s" % filename, 'w')
        # get event information needed for FOCMEC-input file
        eventsource = eventinfo.origins[0]
        stime = eventsource['time']
        # write header line including event information
        fid.write('%s %d%02d%02d%02d%02d%02.0f %7.4f %6.4f %3.1f %3.1f\n' % (parameter.get('eventID'),
                                                                             stime.year, stime.month, stime.day,
                                                                             stime.hour, stime.minute, stime.second,
                                                                             eventsource['latitude'],
                                                                             eventsource['longitude'],
                                                                             eventsource['depth'] / 1000,
                                                                             eventinfo.magnitudes[0]['mag']))
        picks = eventinfo.picks
        for key in arrivals:
            if arrivals[key].has_key('P'):
                if arrivals[key]['P']['weight'] < 4 and arrivals[key]['P']['fm'] is not None:
                    stat = key
                    for i in range(len(picks)):
                        station = picks[i].waveform_id.station_code
                        if station == stat:
                            # get resource ID
                            resid_picks = picks[i].get('resource_id')
                            # find same ID in eventinfo
                            # there it is the pick_id!!
                            for j in range(len(eventinfo.origins[0].arrivals)):
                                resid_eventinfo = eventinfo.origins[0].arrivals[j].get('pick_id')
                                if resid_eventinfo == resid_picks and eventinfo.origins[0].arrivals[j].phase == 'P':
                                    if len(stat) > 4:  # FOCMEC handles only 4-string station IDs
                                        stat = stat[1:5]
                                    az = eventinfo.origins[0].arrivals[j].get('azimuth')
                                    inz = eventinfo.origins[0].arrivals[j].get('takeoff_angle')
                                    fid.write('%-4s  %6.2f  %6.2f%s \n' % (stat,
                                                                           az,
                                                                           inz,
                                                                           arrivals[key]['P']['fm']))
                                    break

        fid.close()

    elif fformat == 'HASH':
        # two different input files for 
        # HASH-driver 1 and 2 (see HASH manual!)
        filename1 = filename + 'drv1' + '.phase'
        filename2 = filename + 'drv2' + '.phase'
        print("Writing phases to %s for HASH for HASH-driver 1" % filename1)
        fid1 = open("%s" % filename1, 'w')
        print("Writing phases to %s for HASH for HASH-driver 2" % filename2)
        fid2 = open("%s" % filename2, 'w')
        # get event information needed for HASH-input file
        eventsource = eventinfo.origins[0]
        event = parameter.get('eventID')
        hashID = event.split('.')[0][1:5]
        latdeg = eventsource['latitude']
        latmin = eventsource['latitude'] * 60 / 10000
        londeg = eventsource['longitude']
        lonmin = eventsource['longitude'] * 60 / 10000
        erh = 1 / 2 * (eventsource.origin_uncertainty['min_horizontal_uncertainty'] +
                       eventsource.origin_uncertainty['max_horizontal_uncertainty']) / 1000
        erz = eventsource.depth_errors['uncertainty']
        stime = eventsource['time']
        if stime.year - 2000 >= 0:
            syear = stime.year - 2000
        else:
            syear = stime.year - 1900
        picks = eventinfo.picks
        # write header line including event information
        # for HASH-driver 1
        fid1.write('%s%02d%02d%02d%02d%5.2f%2dN%5.2f%3dE%5.2f%6.3f%4.2f%5.2f%5.2f%s\n' % (syear,
                                                                                          stime.month, stime.day,
                                                                                          stime.hour, stime.minute,
                                                                                          stime.second,
                                                                                          latdeg, latmin, londeg,
                                                                                          lonmin, eventsource['depth'],
                                                                                          eventinfo.magnitudes[0][
                                                                                              'mag'], erh, erz,
                                                                                          hashID))
        # write header line including event information
        # for HASH-driver 2
        fid2.write(
            '%d%02d%02d%02d%02d%5.2f%dN%5.2f%3dE%6.2f%5.2f    %d                                          %5.2f %5.2f                                        %4.2f      %s \n' % (
                syear, stime.month, stime.day,
                stime.hour, stime.minute, stime.second,
                latdeg, latmin, londeg, lonmin,
                eventsource['depth'],
                eventsource['quality']['used_phase_count'],
                erh, erz, eventinfo.magnitudes[0]['mag'],
                hashID))

        # write phase lines
        for key in arrivals:
            if arrivals[key].has_key('P'):
                if arrivals[key]['P']['weight'] < 4 and arrivals[key]['P']['fm'] is not None:
                    stat = key
                    ccode = arrivals[key]['P']['channel']
                    ncode = arrivals[key]['P']['network']

                    if arrivals[key]['P']['weight'] < 2:
                        Pqual = 'I'
                    else:
                        Pqual = 'E'

                    for i in range(len(picks)):
                        station = picks[i].waveform_id.station_code
                        if station == stat:
                            # get resource ID
                            resid_picks = picks[i].get('resource_id')
                            # find same ID in eventinfo
                            # there it is the pick_id!!
                            for j in range(len(eventinfo.origins[0].arrivals)):
                                resid_eventinfo = eventinfo.origins[0].arrivals[j].get('pick_id')
                                if resid_eventinfo == resid_picks and eventinfo.origins[0].arrivals[j].phase == 'P':
                                    if len(stat) > 4:  # HASH handles only 4-string station IDs
                                        stat = stat[1:5]
                                    az = eventinfo.origins[0].arrivals[j].get('azimuth')
                                    inz = eventinfo.origins[0].arrivals[j].get('takeoff_angle')
                                    dist = eventinfo.origins[0].arrivals[j].get('distance')
                                    # write phase line for HASH-driver 1
                                    fid1.write(
                                        '%-4s%sP%s%d                   0                                   %3.1f          %03d %03d   2     1   %s\n' % (
                                            stat, Pqual, arrivals[key]['P']['fm'], arrivals[key]['P']['weight'],
                                            dist, inz, az, ccode))
                                    # write phase line for HASH-driver 2
                                    fid2.write('%-4s %s   %s %s %s                    \n' % (
                                        stat,
                                        ncode,
                                        ccode,
                                        Pqual,
                                        arrivals[key]['P']['fm']))
                                    break

        fid1.write('                                    %s' % hashID)
        fid1.close()
        fid2.close()


def merge_picks(event, picks):
    """
    takes an event object and a list of picks and searches for matching
    entries by comparing station name and phase_hint and overwrites the time
    and time_errors value of the event picks' with those from the picks
    without changing the resource identifiers
    :param event: `obspy.core.event.Event` object (e.g. from NLLoc output)
    :param picks: list of `obspy.core.event.Pick` objects containing the
    original time and time_errors values
    :return: merged `obspy.core.event.Event` object
    """
    for pick in picks:
        time = pick.time
        err = pick.time_errors
        phase = pick.phase_hint
        station = pick.waveform_id.station_code
        network = pick.waveform_id.network_code
        method = pick.method_id
        for p in event.picks:
            if p.waveform_id.station_code == station\
                    and p.waveform_id.network_code == network\
                    and p.phase_hint == phase\
                    and (str(p.method_id) in str(method)
                         or str(method) in str(p.method_id)):
                p.time, p.time_errors, p.waveform_id.network_code, p.method_id = time, err, network, method
        del time, err, phase, station, network, method
    return event


def getQualitiesfromxml(xmlnames, ErrorsP, ErrorsS, plotflag=1):
    """
    Script to get onset uncertainties from Quakeml.xml files created by PyLoT.
   Uncertainties are tranformed into quality classes and visualized via histogram if desired.
   Ludger Küperkoch, BESTEC GmbH, 07/2017
    :param xmlnames: list of xml obspy event files containing picks
    :type xmlnames: list
    :param ErrorsP: time errors of P waves for the four discrete quality classes
    :type ErrorsP:
    :param ErrorsS: time errors of S waves for the four discrete quality classes
    :type ErrorsS:
    :param plotflag:
    :type plotflag:
    :return:
    :rtype:
    """

    from pylot.core.pick.utils import getQualityFromUncertainty
    from pylot.core.util.utils import loopIdentifyPhase, identifyPhase

    # read all onset weights
    Pw0 = []
    Pw1 = []
    Pw2 = []
    Pw3 = []
    Pw4 = []
    Sw0 = []
    Sw1 = []
    Sw2 = []
    Sw3 = []
    Sw4 = []
    for names in xmlnames:
        print("Getting onset weights from {}".format(names))
        cat = read_events(names)
        cat_copy = cat.copy()
        arrivals = cat.events[0].picks
        arrivals_copy = cat_copy.events[0].picks
        # Prefere manual picks if qualities are sufficient!
        for Pick in arrivals:
            if (Pick.method_id.id).split('/')[1] == 'manual':
                mstation = Pick.waveform_id.station_code
                mstation_ext = mstation + '_'
                for mpick in arrivals_copy:
                    phase = identifyPhase(loopIdentifyPhase(Pick.phase_hint))
                    if phase == 'P':
                        if ((mpick.waveform_id.station_code == mstation) or
                                (mpick.waveform_id.station_code == mstation_ext)) and \
                                ((mpick.method_id).split('/')[1] == 'auto') and \
                                (mpick.time_errors['uncertainty'] <= ErrorsP[3]):
                            del mpick
                            break
                    elif phase == 'S':
                        if ((mpick.waveform_id.station_code == mstation) or
                                (mpick.waveform_id.station_code == mstation_ext)) and \
                                ((mpick.method_id).split('/')[1] == 'auto') and \
                                (mpick.time_errors['uncertainty'] <= ErrorsS[3]):
                            del mpick
                            break
        lendiff = len(arrivals) - len(arrivals_copy)
        if lendiff is not 0:
            print("Found manual as well as automatic picks, prefered the {} manual ones!".format(lendiff))

        for Pick in arrivals_copy:
            phase = identifyPhase(loopIdentifyPhase(Pick.phase_hint))
            if phase == 'P':
                Pqual = getQualityFromUncertainty(Pick.time_errors.uncertainty, ErrorsP)
                if Pqual == 0:
                    Pw0.append(Pick.time_errors.uncertainty)
                elif Pqual == 1:
                    Pw1.append(Pick.time_errors.uncertainty)
                elif Pqual == 2:
                    Pw2.append(Pick.time_errors.uncertainty)
                elif Pqual == 3:
                    Pw3.append(Pick.time_errors.uncertainty)
                elif Pqual == 4:
                    Pw4.append(Pick.time_errors.uncertainty)
            elif phase == 'S':
                Squal = getQualityFromUncertainty(Pick.time_errors.uncertainty, ErrorsS)
                if Squal == 0:
                    Sw0.append(Pick.time_errors.uncertainty)
                elif Squal == 1:
                    Sw1.append(Pick.time_errors.uncertainty)
                elif Squal == 2:
                    Sw2.append(Pick.time_errors.uncertainty)
                elif Squal == 3:
                    Sw3.append(Pick.time_errors.uncertainty)
                elif Squal == 4:
                    Sw4.append(Pick.time_errors.uncertainty)
            else:
                print("Phase hint not defined for picking!")
                pass

    if plotflag == 0:
        Punc = [Pw0, Pw1, Pw2, Pw3, Pw4]
        Sunc = [Sw0, Sw1, Sw2, Sw3, Sw4]
        return Punc, Sunc
    else:
        # get percentage of weights
        numPweights = np.sum([len(Pw0), len(Pw1), len(Pw2), len(Pw3), len(Pw4)])
        numSweights = np.sum([len(Sw0), len(Sw1), len(Sw2), len(Sw3), len(Sw4)])
        if len(Pw0) > 0:
            P0perc = 100 / numPweights * len(Pw0)
        else:
            P0perc = 0
        if len(Pw1) > 0:
            P1perc = 100 / numPweights * len(Pw1)
        else:  
            P1perc = 0
        if len(Pw2) > 0:
            P2perc = 100 / numPweights * len(Pw2)
        else:  
            P2perc = 0
        if len(Pw3) > 0:
            P3perc = 100 / numPweights * len(Pw3)
        else:  
            P3perc = 0
        if len(Pw4) > 0:
            P4perc = 100 / numPweights * len(Pw4)
        else:  
            P4perc = 0
        if len(Sw0) > 0:
            S0perc = 100 / numSweights * len(Sw0)
        else:
            S0perc = 0
        if len(Sw1) > 0:
            S1perc = 100 / numSweights * len(Sw1)
        else:  
            S1perc = 0
        if len(Sw2) > 0:
            S2perc = 100 / numSweights * len(Sw2)
        else:  
            S2perc = 0
        if len(Sw3) > 0:
            S3perc = 100 / numSweights * len(Sw3)
        else:  
            S3perc = 0
        if len(Sw4) > 0:
            S4perc = 100 / numSweights * len(Sw4)
        else:  
            S4perc = 0

        weights = ('0', '1', '2', '3', '4')
        y_pos = np.arange(len(weights))
        width = 0.34
        plt.bar(y_pos - width, [P0perc, P1perc, P2perc, P3perc, P4perc], width, color='black')
        plt.bar(y_pos, [S0perc, S1perc, S2perc, S3perc, S4perc], width, color='red')
        plt.ylabel('%')
        plt.xticks(y_pos, weights)
        plt.xlim([-0.5, 4.5])
        plt.xlabel('Qualities')
        plt.title('{0} P-Qualities, {1} S-Qualities'.format(numPweights, numSweights))
        plt.show()
