#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import obspy.core.event as ope
import os
import scipy.io as sio
import warnings
from obspy.core import UTCDateTime

from pylot.core.io.inputs import AutoPickParameter
from pylot.core.io.location import create_arrival, create_event, \
    create_magnitude, create_origin, create_pick
from pylot.core.pick.utils import select_for_phase
from pylot.core.util.utils import getOwner, getGlobalTimes, four_digits


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
    from pylot.core.util.defaults import TIMEERROR_DEFAULTS
    picks = dict()
    phases_pilot = sio.loadmat(fn)
    stations = stations_from_pilot(phases_pilot['stat'])
    params = AutoPickParameter(TIMEERROR_DEFAULTS)
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
            phases[onset_name] = dict(mpp=pick, spe=spe)
        picks[station] = phases

    return picks


def stations_from_pilot(stat_array):
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


def picksdict_from_picks(evt):
    """
    Takes an Event object and return the pick dictionary commonly used within
    PyLoT
    :param evt: Event object contain all available information
    :type evt: `~obspy.core.event.Event`
    :return: pick dictionary
    """
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
        spe = pick.time_errors.uncertainty
        try:
            lpp = mpp + pick.time_errors.upper_uncertainty
            epp = mpp - pick.time_errors.lower_uncertainty
        except TypeError as e:
            msg = e.message + ',\n falling back to symmetric uncertainties'
            warnings.warn(msg)
            lpp = mpp + spe
            epp = mpp - spe
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

def picks_from_picksdict(picks, creation_info=None):
    picks_list = list()
    for station, onsets in picks.items():
        for label, phase in onsets.items():
            if not isinstance(phase, dict):
                continue
            onset = phase['mpp']
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
            except KeyError as e:
                warnings.warn(e.message, RuntimeWarning)
            try:
                picker = phase['picker']
            except KeyError as e:
                warnings.warn(e.message, RuntimeWarning)
                picker = 'Unknown'
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
                if 'fm' in e.message: # no polarity information found for this phase
                    pass
                else:
                    raise e
            picks_list.append(pick)
    return picks_list


def reassess_pilot_db(root_dir, db_dir, out_dir=None, fn_param=None, verbosity=0):
    import glob

    db_root = os.path.join(root_dir, db_dir)
    evt_list = glob.glob1(db_root,'e????.???.??')

    for evt in evt_list:
        if verbosity > 0:
            print('Reassessing event {0}'.format(evt))
        reassess_pilot_event(root_dir, db_dir, evt, out_dir, fn_param, verbosity)



def reassess_pilot_event(root_dir, db_dir, event_id, out_dir=None, fn_param=None, verbosity=0):
    from obspy import read

    from pylot.core.io.inputs import AutoPickParameter
    from pylot.core.pick.utils import earllatepicker

    if fn_param is None:
        import pylot.core.util.defaults as defaults
        fn_param = defaults.AUTOMATIC_DEFAULTS

    default = AutoPickParameter(fn_param, verbosity)

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
                            warnings.warn('no waveform data found for station {station}'.format(station=station), RuntimeWarning)
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
            stime, etime = getGlobalTimes(sel_st)
            rel_pick = mpp - stime
            epp, lpp, spe = earllatepicker(sel_st,
                                           default.get('nfac{0}'.format(phase)),
                                           default.get('tsnrz' if phase == 'P' else 'tsnrh'),
                                           Pick1=rel_pick,
                                           iplot=None,
                                           stealth_mode=True)
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
        fnout_prefix = os.path.join(root_dir, db_dir, event_id, '{0}.'.format(event_id))
    else:
        out_dir = os.path.join(out_dir, db_dir)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        fnout_prefix = os.path.join(out_dir, '{0}.'.format(event_id))
    evt.write(fnout_prefix + 'xml', format='QUAKEML')
    #evt.write(fnout_prefix + 'cnv', format='VELEST')


def writephases(arrivals, fformat, filename):
    """
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
    """

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
