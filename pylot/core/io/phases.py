#!/usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import logging
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
from pylot.core.pick.utils import select_for_phase, get_quality_class
from pylot.core.util.utils import get_owner, full_range, four_digits, transformFilterString4Export, \
    backtransformFilterString, loopIdentifyPhase, identifyPhase


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
        phasauthor = get_owner(phasfn)
    else:
        phases = None
        phasctime = None
        phasauthor = None
    if locfn is not None and os.path.isfile(locfn):
        loc = sio.loadmat(locfn)
        locctime = UTCDateTime(os.path.getmtime(locfn))
        locauthor = get_owner(locfn)
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


def picksdict_from_picks(evt, parameter=None):
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
        errors = None
        phase = {}
        station = pick.waveform_id.station_code
        if pick.waveform_id.channel_code is None:
            channel = ''
        else:
            channel = pick.waveform_id.channel_code
        network = pick.waveform_id.network_code
        mpp = pick.time
        spe = pick.time_errors.uncertainty
        if pick.filter_id:
            filter_id = backtransformFilterString(str(pick.filter_id.id))
        else:
            filter_id = None
        try:
            pick_method = str(pick.method_id)
            if pick_method.startswith('smi:local/'):
                pick_method = pick_method.split('smi:local/')[1]
        except IndexError:
            pick_method = 'manual'  # MP MP TODO maybe improve statement
        if pick_method == 'None':
            pick_method = 'manual'
        try:
            onsets = picksdict[pick_method][station]
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
        weight = phase.get('weight')
        if not weight:
            if not parameter:
                logging.warning('Using ')
                logging.warning('Using default input parameter')
                parameter = PylotParameter()
            pick.phase_hint = identifyPhase(pick.phase_hint)
            if pick.phase_hint == 'P':
                errors = parameter['timeerrorsP']
            elif pick.phase_hint == 'S':
                errors = parameter['timeerrorsS']
            if errors:
                weight = get_quality_class(spe, errors)
                phase['weight'] = weight
        phase['channel'] = channel
        phase['network'] = network
        phase['picker'] = pick_method
        if pick.polarity == 'positive':
            phase['fm'] = 'U'
        elif pick.polarity == 'negative':
            phase['fm'] = 'D'
        else:
            phase['fm'] = 'N'
        phase['filter_id'] = filter_id if filter_id is not None else ''

        onsets[pick.phase_hint] = phase.copy()
        picksdict[pick_method][station] = onsets.copy()
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
                warnings.warn(str(e), RuntimeWarning)
                filter_id = ''
            pick.filter_id = filter_id

            try:
                polarity = picks[station][label]['fm']
                if polarity == 'U' or polarity == '+':
                    pick.polarity = 'positive'
                elif polarity == 'D' or polarity == '-':
                    pick.polarity = 'negative'
                else:
                    pick.polarity = 'undecidable'
            except:
                pick.polarity = 'undecidable'
                print("No polarity information available!")
            picks_list.append(pick)
    return picks_list


def writephases(arrivals, fformat, filename, parameter=None, eventinfo=None):
    """
    Writes earthquake phase data to different file formats.

    :param arrivals: Dictionary containing phase information (station ID, phase, first motion, weight, etc.)
    :type arrivals: dict
    :param fformat: File format to write to (e.g., 'NLLoc', 'HYPO71', 'HYPOSAT', 'VELEST', 'HYPODD', 'FOCMEC')
    :type fformat: str
    :param filename: Path and name of the output phase file
    :type filename: str
    :param parameter: Additional parameters for writing the phase data
    :type parameter: object
    :param eventinfo: Event information needed for specific formats like VELEST, FOCMEC, and HASH
    :type eventinfo: obspy.core.event.Event
    """

    def write_nlloc():
        with open(filename, 'w') as fid:
            fid.write('# EQEVENT: {} Label: EQ{}  Loc:  X 0.00  Y 0.00  Z 10.00  OT 0.00 \n'.format(
                parameter.get('database'), parameter.get('eventID')))
            for key, value in arrivals.items():
                for phase in ['P', 'S']:
                    if phase in value:
                        fm = value[phase].get('fm', '?')
                        onset = value[phase]['mpp']
                        ss_ms = onset.second + onset.microsecond / 1000000.0
                        weight = 1 if value[phase].get('weight', 0) < 4 else 0
                        amp = value[phase].get('Ao', 0.0) if phase == 'S' else ''
                        fid.write('{} ? ? ? {}   {}{}{} {}{} {:7.4f} GAU 0 {} 0 0 {}\n'.format(
                            key, phase, fm, onset.year, onset.month, onset.day, onset.hour, onset.minute, ss_ms, amp,
                            weight))

    def write_hypo71():
        with open(filename, 'w') as fid:
            fid.write(
                '                                                                {}\n'.format(parameter.get('eventID')))
            for key, value in arrivals.items():
                if value['P'].get('weight', 0) < 4:
                    stat = key[:4]
                    Ponset = value['P']['mpp']
                    Sonset = value.get('S', {}).get('mpp')
                    pweight = value['P'].get('weight', 0)
                    sweight = value.get('S', {}).get('weight', 0)
                    fm = value['P'].get('fm', '-')
                    Ao = value.get('S', {}).get('Ao', '')
                    year = Ponset.year - 2000 if Ponset.year >= 2000 else Ponset.year - 1900
                    ss_ms = Ponset.second + Ponset.microsecond / 1000000.0
                    if Sonset:
                        Sss_ms = Sonset.second + Sonset.microsecond / 1000000.0
                        fid.write('{}P{}{}{} {}{}{}{}{} {:5.2f}       {}{}S {}   {}\n'.format(
                            stat, 'I' if pweight < 2 else 'E', fm, pweight, year, Ponset.month, Ponset.day,
                            Ponset.hour, Ponset.minute, ss_ms, Sss_ms, 'I' if sweight < 2 else 'E', sweight, Ao))
                    else:
                        fid.write('{}P{}{}{} {}{}{}{}{} {:5.2f}                  {}\n'.format(
                            stat, 'I' if pweight < 2 else 'E', fm, pweight, year, Ponset.month, Ponset.day,
                            Ponset.hour, Ponset.minute, ss_ms, Ao))

    def write_hyposat():
        with open(filename, 'w') as fid:
            fid.write('{}, event {} \n'.format(parameter.get('database'), parameter.get('eventID')))
            for key, value in arrivals.items():
                for phase in ['P', 'S']:
                    if phase in value and value[phase].get('weight', 0) < 4:
                        onset = value[phase]['mpp']
                        ss_ms = onset.second + onset.microsecond / 1000000.0
                        std = value[phase].get('spe', parameter.get('timeerrorsP')[value[phase].get('weight', 0)])
                        fid.write(
                            '{:<5} {}1       {:4} {:02} {:02} {:02} {:02} {:05.02f}   {:5.3f} -999.   0.00 -999.  0.00\n'.format(
                                key, phase, onset.year, onset.month, onset.day, onset.hour, onset.minute, ss_ms, std))

    def write_velest():
        if not eventinfo:
            print("No source origin calculated yet, thus no cnv-file creation possible!")
            return
        with open(filename, 'w') as fid:
            origin = eventinfo.origins[0]
            lat_dir = 'S' if origin.latitude < 0 else 'N'
            lon_dir = 'W' if origin.longitude < 0 else 'E'
            year = origin.time.year - 2000 if origin.time.year >= 2000 else origin.time.year - 1900
            fid.write(
                '{}{}{} {}{} {} {:05.2f} {:7.4f}{} {:8.4f}{} {:7.2f} {:6.2f}     {:02.0f}  0.0 0.03  1.0  1.0\n'.format(
                    year, origin.time.month, origin.time.day, origin.time.hour, origin.time.minute, origin.time.second,
                    origin.latitude, lat_dir, origin.longitude, lon_dir, origin.depth, eventinfo.magnitudes[0].mag, 0))
            for key, value in arrivals.items():
                for phase in ['P', 'S']:
                    if phase in value and value[phase].get('weight', 0) < 4:
                        onset = value[phase]['mpp']
                        rt = (onset - origin.time).total_seconds()
                        fid.write('{:<4}{}{}{:6.2f}\n'.format(key[:4], phase, value[phase].get('weight', 0), rt))

    def write_hypodd():
        if not eventinfo:
            print("No source origin calculated yet, thus no hypoDD-infile creation possible!")
            return
        with open(filename, 'w') as fid:
            origin = eventinfo.origins[0]
            stime = origin.time
            fid.write('# {}  {} {} {} {} {} {:7.4f} +{:6.4f} {:7.4f} {:4.2f} 0.1 0.5 {:4.2f}      {}\n'.format(
                stime.year, stime.month, stime.day, stime.hour, stime.minute, stime.second,
                origin.latitude, origin.longitude, origin.depth / 1000, eventinfo.magnitudes[0].mag,
                origin.quality.standard_error, "00000"))
            for key, value in arrivals.items():
                for phase in ['P', 'S']:
                    if phase in value and value[phase].get('weight', 0) < 4:
                        onset = value[phase]['mpp']
                        rt = (onset - stime).total_seconds()
                        fid.write('{}    {:6.3f}  1  {}\n'.format(key, rt, phase))

    def write_focmec():
        if not eventinfo:
            print("No source origin calculated yet, thus no FOCMEC-infile creation possible!")
            return
        with open(filename, 'w') as fid:
            origin = eventinfo.origins[0]
            stime = origin.time
            fid.write('{} {}{:02d}{:02d}{:02d}{:02d}{:02.0f} {:7.4f} {:6.4f} {:3.1f} {:3.1f}\n'.format(
                parameter.get('eventid', 'e0000'), stime.year, stime.month, stime.day, stime.hour, stime.minute,
                stime.second, origin.latitude, origin.longitude, origin.depth / 1000, eventinfo.magnitudes[0].mag))
            for key, value in arrivals.items():
                if 'P' in value and value['P'].get('weight', 0) < 4 and value['P'].get('fm'):
                    for pick in eventinfo.picks:
                        if pick.waveform_id.station_code == key:
                            for arrival in origin.arrivals:
                                if arrival.pick_id == pick.resource_id and arrival.phase == 'P':
                                    stat = key[:4]
                                    az = arrival.azimuth
                                    inz = arrival.takeoff_angle
                                    fid.write('{:<4}  {:6.2f}  {:6.2f}{}\n'.format(stat, az, inz, value['P']['fm']))
                                    break

    if fformat == 'NLLoc':
        write_nlloc()
    elif fformat == 'HYPO71':
        write_hypo71()
    elif fformat == 'HYPOSAT':
        write_hyposat()
    elif fformat == 'VELEST':
        write_velest()
    elif fformat == 'HYPODD':
        write_hypodd()
    elif fformat == 'FOCMEC':
        write_focmec()



def chooseArrivals(arrivals):
    """
    takes arrivals and returns the manual picks if manual and automatic ones are there
    returns automatic picks if only automatic picks are there
    :param arrivals: 'dictionary' with automatic and or manual arrivals
    :return: arrivals but with the manual picks prefered if possible
    """
    # If len of arrivals is greater than 2 it comes from autopicking so only autopicks are available
    if len(arrivals) > 2:
        return arrivals
    if arrivals['auto'] and arrivals['manual']:
        usedarrivals = arrivals['manual']
    elif arrivals['auto']:
        usedarrivals = arrivals['auto']
    elif arrivals['manual']:
        usedarrivals = arrivals['manual']
    return usedarrivals


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
            if p.waveform_id.station_code == station \
                    and p.waveform_id.network_code == network \
                    and p.phase_hint == phase \
                    and (str(p.method_id) in str(method)
                         or str(method) in str(p.method_id)):
                p.time, p.time_errors, p.waveform_id.network_code, p.method_id = time, err, network, method
        del time, err, phase, station, network, method
    return event


def getQualitiesfromxml(path, errorsP, errorsS, plotflag=1, figure=None, verbosity=0):
    """
    Script to get onset uncertainties from Quakeml.xml files created by PyLoT.
    Uncertainties are tranformed into quality classes and visualized via histogram if desired.
    Ludger Küperkoch, BESTEC GmbH, 07/2017
    :param path: path containing xml files
    :type path: str
    :param errorsP: time errors of P waves for the four discrete quality classes
    :type errorsP:
    :param errorsS: time errors of S waves for the four discrete quality classes
    :type errorsS:
    :param plotflag:
    :type plotflag:
    :return:
    :rtype:
    """

    def calc_perc(uncertainties, ntotal):
        ''' simple function that calculates percentage of number of uncertainties (list length)'''
        if len(uncertainties) == 0:
            return 0
        else:
            return 100. / ntotal * len(uncertainties)

    def calc_weight_perc(psweights, weight_ids):
        ''' calculate percentages of different weights (pick classes!?) of total number of uncertainties of a phase'''
        # count total number of list items for this phase
        numWeights = np.sum([len(weight) for weight in psweights.values()])

        # iterate over all available weights to return a list with percentages for plotting
        plot_list = []
        for weight_id in weight_ids:
            plot_list.append(calc_perc(psweights[weight_id], numWeights))

        return plot_list, numWeights

    # get all xmlfiles in path (maybe this should be changed to one xml file for this function, selectable via GUI?)
    xmlnames = glob.glob(os.path.join(path, '*.xml'))
    if len(xmlnames) == 0:
        print(f'No files found in path {path}.')
        return False

    # first define possible phases here
    phases = ['P', 'S']

    # define possible weights (0-4)
    weight_ids = list(range(5))

    # put both error lists in a dictionary with P/S key so that amount of code can be halfed by simply using P/S as key
    errors = dict(P=errorsP, S=errorsS)

    # create dictionaries for each phase (P/S) with a dictionary of empty list for each weight defined in weights
    # tuple above
    weights = {}
    for phase in phases:
        weights[phase] = {weight_id: [] for weight_id in weight_ids}

    for names in xmlnames:
        print("Getting onset weights from {}".format(names))
        cat = read_events(names)
        cat_copy = cat.copy()
        arrivals = cat.events[0].picks
        arrivals_copy = cat_copy.events[0].picks
        # Prefere manual picks if qualities are sufficient!
        for pick in arrivals:
            if pick.method_id.id.split('/')[1] == 'manual':
                mstation = pick.waveform_id.station_code
                mstation_ext = mstation + '_'
                for mpick in arrivals_copy:
                    phase = identifyPhase(loopIdentifyPhase(pick.phase_hint)) # MP MP catch if this fails?
                    if ((mpick.waveform_id.station_code == mstation) or
                        (mpick.waveform_id.station_code == mstation_ext)) and \
                            (mpick.method_id.id.split('/')[1] == 'auto') and \
                            (mpick.time_errors['uncertainty'] <= errors[phase][3]):
                        del mpick
                        break
        lendiff = len(arrivals) - len(arrivals_copy)
        if lendiff != 0:
            print("Found manual as well as automatic picks, prefered the {} manual ones!".format(lendiff))

        for pick in arrivals_copy:
            phase = identifyPhase(loopIdentifyPhase(pick.phase_hint))
            uncertainty = pick.time_errors.uncertainty
            if not uncertainty:
                if verbosity > 0:
                    print('No uncertainty, pick {} invalid!'.format(pick.method_id.id))
                continue
            # check P/S phase
            if phase not in phases:
                print("Phase hint not defined for picking!")
                continue

            qual = get_quality_class(uncertainty, errors[phase])
            weights[phase][qual].append(uncertainty)

    if plotflag == 0:
        p_unc = [weights['P'][weight_id] for weight_id in weight_ids]
        s_unc = [weights['S'][weight_id] for weight_id in weight_ids]
        return p_unc, s_unc
    else:
        if not figure:
            fig = plt.figure()
        ax = fig.add_subplot(111)
        # get percentage of weights
        listP, numPweights = calc_weight_perc(weights['P'], weight_ids)
        listS, numSweights = calc_weight_perc(weights['S'], weight_ids)

        y_pos = np.arange(len(weight_ids))
        width = 0.34
        ax.bar(y_pos - width, listP, width, color='black')
        ax.bar(y_pos, listS, width, color='red')
        ax.set_ylabel('%')
        ax.set_xticks(y_pos, weight_ids)
        ax.set_xlim([-0.5, 4.5])
        ax.set_xlabel('Qualities')
        ax.set_title('{0} P-Qualities, {1} S-Qualities'.format(numPweights, numSweights))

        if not figure:
            fig.show()

        return listP, listS
