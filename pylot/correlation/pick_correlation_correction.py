#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
import copy
import logging
import random
import traceback

import matplotlib.pyplot as plt
import yaml
from code_base.fmtomo_tools.fmtomo_teleseismic_utils import *
from code_base.util.utils import get_metadata
from joblib import Parallel, delayed
from obspy import read, Stream
from obspy.core.event.base import WaveformStreamID, ResourceIdentifier
from obspy.core.event.origin import Pick
from obspy.signal.cross_correlation import correlate, xcorr_max
from pylot.core.io.inputs import PylotParameter
from pylot.core.io.phases import picks_from_picksdict
from pylot.core.pick.autopick import autopickstation
from pylot.core.util.utils import check4rotated
from pylot.core.util.utils import identifyPhaseID
from pylot.core.util.event import Event as PylotEvent

class CorrelationParameters(object):
    """
    :param filter_type: filter type for data (e.g. bandpass)
    :param filter_options: filter options (min/max_freq etc.)
    :param filter_options_final: filter options for second iteration, final pick (min/max_freq etc.)
    :param freq: resampling frequence (IMPORTANT for correlation)
    :param station_list: list of stations used as possible reference stations
    :param min_corr_stacking: minimum correlation coefficient for stacking (between 0 and 1)
    :param min_corr_export: minimum correlation coefficient for export (between 0 and 1)
    :param min_stack: minimum number of stacks for exporting of picks
    :param t_before: time before initial pick taken into account
    :param t_after: time after initial pick taken into account
    :param initial_pick_outlier_threshold: remove initial picks that differ more than x seconds from median corrected
     taupy onset
     :param export_threshold: remove final correlated picks that differ more than x seconds from (newly calculated)
      median corrected taupy onsets
    :param cc_maxlag: cross correlation shift
    :param plot: plot results
    :param plot_detailed: make detailed plots (of every trace)
    :param save_fig: save figures
    :param use_taupy_onsets: use taupy theoretical onsets instead of initial picks from folders (bool)
    :param ncores: number of processed spawned for correlation (multiprocessing)
    :param use_stacked_trace: use existing trace of previous correlation (e.g. when adding more data to an already
    evaluated event)
    :param data_dir: obspy_dmt data dir (processed/raw)
    """

    def __init__(self, **kwargs):
        self.__parameter = {}
        self.add_parameters(**kwargs)

    # String representation of the object
    def __repr__(self):
        return "CorrelationParameters "

    # Boolean test
    def __nonzero__(self):
        return bool(self.__parameter)

    def __getitem__(self, key):
        return self.__parameter.get(key)

    def __setitem__(self, key, value):
        self.__parameter[key] = value

    def __delitem__(self, key):
        del self.__parameter[key]

    def __iter__(self):
        return iter(self.__parameter)

    def __len__(self):
        return len(self.__parameter.keys())

    def add_parameters(self, **kwargs):
        for key, value in kwargs.items():
            self.__parameter[key] = value


class Xcorr_pick_correction:
    def __init__(self, pick1, trace1, pick2, trace2, t_before, t_after, cc_maxlag, frac_max=0.5):
        """
        MP MP : Modified version of obspy xcorr_pick_correction

        Calculate the correction for the differential pick time determined by cross
        correlation of the waveforms in narrow windows around the pick times.
        For details on the fitting procedure refer to [Deichmann1992]_.

        The parameters depend on the epicentral distance and magnitude range. For
        small local earthquakes (Ml ~0-2, distance ~3-10 km) with consistent manual
        picks the following can be tried::

            t_before=0.05, t_after=0.2, cc_maxlag=0.10,

        The appropriate parameter sets can and should be determined/verified
        visually using the option `plot=True` on a representative set of picks.

        To get the corrected differential pick time calculate: ``((pick2 +
        pick2_corr) - pick1)``. To get a corrected differential travel time using
        origin times for both events calculate: ``((pick2 + pick2_corr - ot2) -
        (pick1 - ot1))``

        :type pick1: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param pick1: Time of pick for `trace1`.
        :type trace1: :class:`~obspy.core.trace.Trace`
        :param trace1: Waveform data for `pick1`. Add some time at front/back.
                The appropriate part of the trace is used automatically.
        :type pick2: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param pick2: Time of pick for `trace2`.
        :type trace2: :class:`~obspy.core.trace.Trace`
        :param trace2: Waveform data for `pick2`. Add some time at front/back.
                The appropriate part of the trace is used automatically.
        :type t_before: float
        :param t_before: Time to start cross correlation window before pick times
                in seconds.
        :type t_after: float
        :param t_after: Time to end cross correlation window after pick times in
                seconds.
        :type cc_maxlag: float
        :param cc_maxlag: Maximum lag/shift time tested during cross correlation in
            seconds.
        :rtype: (float, float)
        :returns: Correction time `pick2_corr` for `pick2` pick time as a float and
                corresponding correlation coefficient.
        """

        self.trace1 = trace1
        self.trace2 = trace2

        self.pick1 = pick1
        self.pick2 = pick2

        self.cc_maxlag = cc_maxlag
        self.t_before = t_before
        self.t_after = t_after
        self.frac_max = frac_max

        self.samp_rate = 0

        # perform some checks on the traces
        self.check_traces()

        # check data and take correct slice of traces
        self.tr1_slice = self.slice_trace(self.trace1, self.pick1)
        self.tr2_slice = self.slice_trace(self.trace2, self.pick2)

    def check_traces(self):
        if self.trace1.stats.sampling_rate != self.trace2.stats.sampling_rate:
            msg = "Sampling rates do not match: %s != %s" % (
            self.trace1.stats.sampling_rate, self.trace2.stats.sampling_rate)
            raise Exception(msg)
        self.samp_rate = self.trace1.stats.sampling_rate

    def slice_trace(self, tr, pick):
        start = pick - self.t_before - (self.cc_maxlag / 2.0)
        end = pick + self.t_after + (self.cc_maxlag / 2.0)
        # check if necessary time spans are present in data
        if tr.stats.starttime > start:
            msg = f"Trace {tr.id} starts too late."
            raise Exception(msg)
        if tr.stats.endtime < end:
            msg = f"Trace {tr.id} ends too early."
            raise Exception(msg)
        # apply signal processing and take correct slice of data
        return tr.slice(start, end)

    def cross_correlation(self, plot, fig_dir, plot_name, min_corr=None):

        def get_timeaxis(trace):
            return np.linspace(0,trace.stats.endtime -trace.stats.starttime, trace.stats.npts)

        def plot_cc(fig_dir, plot_name):
            if fig_dir and os.path.isdir(fig_dir):
                filename = os.path.join(fig_dir, 'corr_{}_{}.svg'.format(self.trace2.id, plot_name))
            else:
                filename = None
            # with MatplotlibBackend(filename and "AGG" or None, sloppy=True):

            fig = plt.figure(figsize=(16, 9))
            ax1 = fig.add_subplot(211)
            tmp_t1 = get_timeaxis(self.tr1_slice)
            tmp_t2 = get_timeaxis(self.tr2_slice)
            # MP MP normalize slices (not only by positive maximum!
            tr1_slice_norm = self.tr1_slice.copy().normalize()
            tr2_slice_norm = self.tr2_slice.copy().normalize()
            ax1.plot(tmp_t1, tr1_slice_norm.data, "b", label="Trace 1 (reference)", lw=0.75)
            ax1.plot(tmp_t2, tr2_slice_norm.data, "g--", label="Trace 2 (pick shifted)", lw=0.75)
            ax1.plot(tmp_t2 - dt, tr2_slice_norm.data, "k", label="Trace 2 (correlation shifted)", lw=1.)
            delta_pick_ref = (self.pick1 - self.tr1_slice.stats.starttime)
            # correct pick time shift in traces for trace1
            ax1.axvline(delta_pick_ref, color='k', linestyle='dashed', label='Pick', lw=0.5)
            ylims = ax1.get_ylim()
            ax1.fill_between([delta_pick_ref - uncert, delta_pick_ref + uncert], ylims[0], ylims[1], alpha=.25,
                             color='g', label='pick uncertainty)'.format(self.frac_max))
            ax1.legend(loc="lower right", prop={'size': "small"})
            ax1.set_title("Correlated {} with {}".format(self.tr2_slice.id, self.tr1_slice.id))
            ax1.set_xlabel("time [s]")
            ax1.set_ylabel("norm. amplitude")
            ax2 = fig.add_subplot(212)
            ax2.plot(cc_t, cc_convex, ls="", marker=".", color="k", label="xcorr (convex)")
            ax2.plot(cc_t, cc_concave, ls="", marker=".", color="0.7", label="xcorr (concave)")
            ax2.plot(cc_t[first_sample:last_sample + 1], cc[first_sample:last_sample + 1], "b.",
                     label="used for fitting")
            tmp_t = np.linspace(cc_t[first_sample], cc_t[last_sample], num_samples * 10)
            ax2.plot(tmp_t, np.polyval(coeffs, tmp_t), "b", label="fit")
            ax2.axvline(-dt, color="g", label="vertex")
            ax2.axhline(coeff, color="g")
            ax2.hlines(self.frac_max * fmax, -dt - 0.5 * fwfm, -dt + 0.5 * fwfm, color='r', label='FWFM')
            fill_color = 'r' if min_corr and coeff < min_corr else 'c'
            ax2.fill_between([-dt - 0.5 * fwfm, -dt + 0.5 * fwfm], coeff, 1, alpha=.25, color=fill_color,
                             label='pick error')
            ax2.set_xlabel("%.2f at %.3f s correction. Pickerror: %.2f s" % (coeff, -dt, uncert))
            ax2.set_ylabel("correlation coefficient")
            ax2.set_ylim(-1, 1)
            ax2.set_xlim(cc_t[0], cc_t[-1])
            ax2.legend(loc="lower right", prop={'size': "x-small"})
            # plt.legend(loc="lower left")
            if filename:
                fig.savefig(filename)
            else:
                plt.show()

            plt.close(fig)

        # start of cross correlation method
        shift_len = int(self.cc_maxlag * self.samp_rate)
        cc = correlate(self.tr1_slice.data, self.tr2_slice.data, shift_len, demean=False)
        _cc_shift, cc_max = xcorr_max(cc)
        cc_curvature = np.concatenate((np.zeros(1), np.diff(cc, 2), np.zeros(1)))
        cc_convex = np.ma.masked_where(np.sign(cc_curvature) >= 0, cc)
        cc_concave = np.ma.masked_where(np.sign(cc_curvature) < 0, cc)

        # check results of cross correlation
        if cc_max < 0:
            msg = "Trace: {} Absolute maximum is negative: {:.3}. Set cc_max to 0"
            logging.debug(msg.format(self.trace2.id, cc_max))
            cc_max = 0

        # make array with time shifts in seconds corresponding to cc function
        cc_t = np.linspace(-self.cc_maxlag, self.cc_maxlag, len(cc))
        # take the subportion of the cross correlation around the maximum that is
        # convex and fit a parabola.
        # use vertex as subsample resolution best cc fit.
        peak_index = cc.argmax()
        first_sample = peak_index
        # XXX this could be improved..
        while first_sample > 0 and cc_curvature[first_sample - 1] <= 0:
            first_sample -= 1
        last_sample = peak_index
        while last_sample < len(cc) - 1 and cc_curvature[last_sample + 1] <= 0:
            last_sample += 1
        if first_sample == 0 or last_sample == len(cc) - 1:
            msg = "Fitting at maximum lag. Maximum lag time should be increased."
            logging.debug(msg)

        # work on subarrays
        num_samples = last_sample - first_sample + 1
        if num_samples < 3:
            msg = "Less than 3 samples selected for fit to cross " + "correlation: %s" % num_samples
            raise Exception(msg)
        if num_samples < 5:
            msg = "Less than 5 samples selected for fit to cross " + "correlation: %s" % num_samples
            logging.debug(msg)

        # quadratic fit for small subwindow
        coeffs, cov_mat = np.polyfit(cc_t[first_sample:last_sample + 1], cc[first_sample:last_sample + 1], deg=2,
            full=False, cov=True)[:2]

        a, b, c = coeffs

        # check results of fit
        if a >= 0:
            msg = "Fitted parabola opens upwards!"
            logging.debug(msg)

        # LON coordinate of vertex of parabola gives time shift to correct
        # differential pick time. LAT coordinate gives maximum correlation
        # coefficient.
        dt = -b / 2.0 / a
        coeff = (4 * a * c - b ** 2) / (4 * a)

        # this is the shift to apply on the time axis of `trace2` to align the
        # traces. Actually we do not want to shift the trace to align it but we
        # want to correct the time of `pick2` so that the traces align without
        # shifting. This is the negative of the cross correlation shift.
        # MP MP calculate full width at half maximum
        fmax = coeff
        fwfm = 2 * (np.sqrt((b / (2 * a)) ** 2 + (self.frac_max * fmax - c) / a))

        # # calculated error using a and ccc
        # if a < 0:
        #     asqrtcc = np.sqrt(-1. / a) * (1. - coeff)
        # else:
        #     # warning already printed above
        #     asqrtcc = np.nan

        # uncertainty is half of two times the fwhm scaled by 1 - maxcc
        uncert = fwfm * (1. - coeff)
        if uncert < 0:
            logging.warning('negative uncertainty. Check Parabola fit!!')
        dt = -dt
        pick2_corr = dt

        # set coeff to 0 if cc_max was negative
        if cc_max == 0.:
            coeff = 0.

        if plot:
            try:
                plot_cc(fig_dir, plot_name)
            except Exception as e:
                logging.error(f'{e}: {traceback.format_exc()}')

        return pick2_corr, coeff, uncert, fwfm


def correlation_main(database_path_dmt, pylot_infile_path, params, channel_config, istart=0, istop=1e9, update=False,
                     event_blacklist=None):
    '''
    Main function of this program, correlates waveforms around theoretical, or other initial (e.g. automatic cf) picks
    on one or more reference stations.
    All stations with a correlation higher "min_corr_stacking" are stacked onto the station with the highest mean
    correlation. The stacked seismogram will be re-picked using parameters set in "pylot_infile_path".
    Finally, all other stations will be correlated with the re-picked, stacked seismogram around their theoretical (or
    initial) onsets and the shifted pick time of the stacked seismogram will be assigned as their new pick time.

    :param database_path_dmt: obspy_dmt databse directory
    :param pylot_infile_path: path containing input file for autoPyLoT (automatic picking parameters)
    :param params: dictionary with parameter class objects

    :return:
    '''
    assert os.path.isdir(database_path_dmt), 'Unrecognized directory {}'.format(database_path_dmt)

    tstart = datetime.now()
    logging.info(50 * '#')
    logging.info('Starting pick correlation script at {}'.format(tstart))
    logging.info('Eventindices selected: {} - {}'.format(istart, istop))
    logging.info(50 * '#')

    for phase_type in params.keys():
        if params[phase_type]['plot_detailed']: params[phase_type]['plot'] = True

    eventdirs = glob.glob(os.path.join(database_path_dmt, '*.?'))

    pylot_parameter = PylotParameter(pylot_infile_path)

    if event_blacklist:
        with open(event_blacklist, 'r') as infile:
            event_blacklist = [line.split('\n')[0] for line in infile.readlines()]

    # iterate over all events in "database_path_dmt"
    for eventindex, eventdir in enumerate(eventdirs):
        if not istart <= eventindex < istop:
            continue

        # MP MP testing +++
        # ids_filter = ['20181229_033914.a']
        # if not os.path.split(eventdir)[-1] in ids_filter:
        #   continue
        # MP MP testing ---

        logging.info('\n' + 100 * '#')
        logging.info('Working on event {} ({}/{})'.format(eventdir, eventindex + 1, len(eventdirs)))
        if event_blacklist and get_event_id(eventdir) in event_blacklist:
            logging.info('Event on blacklist. Continue')

        correlate_event(eventdir, pylot_parameter, params=params, channel_config=channel_config, update=update)

    logging.info('Finished script after {} at {}'.format(datetime.now() - tstart, datetime.now()))


def get_estimated_inclination(stations_dict, origin, phases, model='ak135'):
    ''' calculate a mean inclination angle for all stations for seismometer rotation to spare computation time'''
    model = TauPyModel(model)
    phases = [*phases.split(',')]
    avg_lat = np.median([item['latitude'] for item in stations_dict.values()])
    avg_lon = np.median([item['longitude'] for item in stations_dict.values()])
    arrivals = model.get_ray_paths_geo(origin.depth, origin.latitude, origin.longitude, avg_lat, avg_lon,
                                       phase_list=phases)
    if len(arrivals) > 0:
        return arrivals[0].incident_angle


def modify_horizontal_name(wfdata):
    ''' This function only renames (does not rotate!) numeric channel IDs to be able to rotate them to LQT. It is
     therefore not accurate and only used for picking with correlation. Returns a copy of the original stream. '''
    # wfdata = wfdata.copy()
    stations = np.unique([tr.stats.station for tr in wfdata])
    for station in stations:
        st = wfdata.select(station=station)
        locations = np.unique([tr.stats.location for tr in st])
        for location in locations:
            st = wfdata.select(station=station, location=location)
            channels = [tr.stats.channel for tr in st]
            check_numeric = [channel[-1].isnumeric() for channel in channels]
            check_nonnumeric = [not (cn) for cn in check_numeric]
            if not any(check_numeric):
                continue
            numeric_channels = np.array(channels)[check_numeric].tolist()
            nonnumeric_channels = np.array(channels)[check_nonnumeric].tolist()
            if not 'Z' in [nc[-1] for nc in nonnumeric_channels]:
                logging.warning('Modify_horizontal_name failed: Only implemented for existing Z component! Return original data.')
                return wfdata
            numeric_characters = sorted([nc[-1] for nc in numeric_channels])
            if numeric_characters == ['1', '2']:
                channel_dict = {'1': 'N', '2': 'E'}
            elif numeric_characters == ['2', '3']:
                channel_dict = {'2': 'N', '3': 'E'}
            else:
                logging.warning('Modify_horizontal_name failed: Channel order not implemented/unknown. Return original data.')
                return wfdata
            for tr in st:
                channel = tr.stats.channel
                for ch_key in channel_dict.keys():
                    if channel[-1] == ch_key:
                        channel = channel[:-1] + channel_dict[ch_key]
                        tr.stats.channel = channel

    return wfdata


def cut_stream_to_same_length(wfdata):
    def remove(wfdata, st):
        for tr in st:
            wfdata.remove(tr)

    stations = np.unique([tr.stats.station for tr in wfdata])
    for station in stations:
        st = wfdata.select(station=station)
        tstart = max([tr.stats.starttime for tr in st])
        tend = min([tr.stats.endtime for tr in st])
        if tstart > tend:
            logging.info('Starttime > Endtime, remove stream from dataset:', str(st))
            remove(wfdata, st)
            continue
        st.trim(tstart, tend, pad=True, fill_value=0.)
        # check for stream length
        if len(st) < 3:
            logging.info('Not enough traces in stream, remove it from dataset:', str(st))
            remove(wfdata, st)
            continue
        if st[0].stats.starttime != st[1].stats.starttime or st[1].stats.starttime != st[2].stats.starttime:
            # in this case, tstart and tend have also changed
            tstart = max([tr.stats.starttime for tr in st])
            tend = min([tr.stats.endtime for tr in st])
            sampling_rate = st[0].stats.sampling_rate
            st.interpolate(sampling_rate=sampling_rate, method='linear', starttime=tstart)
            st.trim(tstart, tend, pad=True, fill_value=0.)


def get_unique_phase(picks, rename_phases=None):
    phases = [pick.phase_hint for pick in picks]
    if rename_phases:
        phases = [rename_phases[phase] if phase in rename_phases else phase for phase in phases]

    if len(set(phases)) == 1:
        return phases[0]

    return False


def correlate_event(eventdir, pylot_parameter, params, channel_config, update):
    rename_phases = {'Pdiff': 'P'}

    # create ObsPy event from .pkl file in dmt eventdir
    event = get_event_obspy_dmt(eventdir)
    # catch event id
    event_id = get_event_id(eventdir)
    if len(event.origins) > 1:
        raise Exception('Ambiguous number of origins for event {}'.format(event))
    origin = event.origins[0]
    # get metadata for current event from dmt_database
    metadata = get_metadata(eventdir)
    # get a dictionary containing coordinates for all sources
    stations_dict = metadata.get_all_coordinates()

    # read processed (restituted) data (assume P and S are the same...)
    wfdata_raw = get_data(eventdir, params['P']['data_dir'], headonly=False)

    # create dictionaries for final export of P and S phases together
    # ps_correlation_dict = {}
    # ps_taup_picks = {}
    # ps_picks = {}
    correlations_dict_stacked = {}

    # iterate over P and S to create one model each
    for phase_type in params.keys():
        ncores = params[phase_type]['ncores']
        filter_type = params[phase_type]['filter_type']
        filter_options = params[phase_type]['filter_options']
        filter_options_final = params[phase_type]['filter_options_final']

        logging.info('\n' + 100 * '#')
        logging.info(f'Working on {phase_type} picks.')
        logging.info('RMS check set to: {}'.format(params[phase_type]['check_RMS']))
        logging.info('Use TauPy onsets instead of picks: {}'.format(params[phase_type]['use_taupy_onsets']))
        logging.info('Use previously generated stacked trace: {}'.format(params[phase_type]['use_stacked_trace']))
        logging.info('Data directory: {}'.format(params[phase_type]['data_dir']))
        logging.info('Pick file extension: {}'.format(params[phase_type]['pickfile_extension']))

        taupypicks_orig = get_taupy_picks(origin, stations_dict, pylot_parameter, phase_type=phase_type,
                                          ncores=params[phase_type]['ncores'])
        if not taupypicks_orig:
            logging.info('No tau-p picks for event {} found.'.format(eventdir))
            continue

        # rename phase to actual phase from tau-p calculation if unique
        true_phase_type = get_unique_phase(taupypicks_orig, rename_phases)
        if not true_phase_type:
            logging.info(f'Ambiguous phase information found in taupy-picks for event {event}. Continue with next one')
            continue

        if update:
            pickfile_path = os.path.join(eventdir,
                                         get_pickfile_name_correlated(event_id, filter_options_final, true_phase_type))
            if os.path.isfile(pickfile_path):
                logging.info(
                    f'Option "update" selected. Found file {pickfile_path}. Skipping phase {phase_type} for this event')
                continue

        wfdata = wfdata_raw.copy()
        # resample and filter
        if params[phase_type]['sampfreq']:
            wfdata = resample_parallel(wfdata, params[phase_type]['sampfreq'], ncores=params[phase_type]['ncores'])
        else:
            logging.warning('Resampling deactivated! '
                            'Make sure that the sampling rate of all input waveforms is identical')

        # pick-file-extension (file naming)
        pfe = params[phase_type]['pickfile_extension']
        # get picks or taupy onsets
        if params[phase_type]['use_taupy_onsets']:
            method = 'taupy'
            picks = taupypicks_orig
        else:
            method = 'auto'
            # get all picks from PyLoT *.xml file
            if not pfe.startswith('_'): pfe = '_' + pfe
            picks = get_picks(eventdir, extension=pfe)
            # take picks of selected phase only
            picks = [pick for pick in picks if pick.phase_hint == phase_type]
            # exclude quality 4 picks (invalid, mostly just beginning of trace)
            picks = remove_invalid_picks(picks)
            if not picks:
                logging.info('No picks for event {} found.'.format(eventdir))
                continue
            if len(picks) < params[phase_type]['min_picks_autopylot']:
                logging.info('Not enough automatic picks for correlation. Continue!')
                continue
            # calculate corrected taupy picks and remove strong outliers
            taupypicks_corr_initial, median_diff = get_corrected_taupy_picks(picks, taupypicks_orig)
            picks = remove_outliers(picks, taupypicks_corr_initial,
                                    params[phase_type]['initial_pick_outlier_threshold'])

            # # MP MP TEST  # taupypicks, time_shift = get_corrected_taupy_picks(picks, taupypicks_orig, all_available=True)  #  # picks = copy.deepcopy(taupypicks)  # for pick in picks:  #     pick.method_id.id = 'auto'  # # MP MP

        if phase_type == 'S':
            # check whether rotation to ZNE is possible (to get rid of horizontal channel 1,2,3)
            # the two rotations here (first to ZNE then LQT) could, however, be optimised into one
            wfdata = check4rotated(wfdata, metadata)
            # if not, try to rename horizontal channels to allow rotation to LQT (however, less accurate)
            wfdata = modify_horizontal_name(wfdata)
            # get an average (one ray path calculation only to spare time) inclination angle
            inclination = get_estimated_inclination(stations_dict, origin, phases=pylot_parameter['taup_phases'],
                                                    model='ak135')
            # rotate ZNE -> LQT (also cut to same length)
            wfdata = rotate_stream(wfdata, metadata=metadata, origin=origin, stations_dict=stations_dict,
                                   channels=channels[phase_type], inclination=inclination, ncores=ncores)

        # make copies before filtering
        wfdata_lowf = wfdata.copy()
        wfdata_highf = wfdata.copy()

        if filter_type and filter_options:
            wfdata_lowf.filter(type=filter_type, **filter_options)
            wfdata_highf.filter(type=filter_type, **filter_options_final)

        channels_list = channel_config[phase_type]

        # make directory for correlation output
        correlation_out_dir = create_correlation_output_dir(eventdir, filter_options_final, true_phase_type)
        fig_dir = create_correlation_figure_dir(correlation_out_dir) if params[phase_type]['save_fig'] == True else ''

        if not params[phase_type]['use_stacked_trace']:
            # first stack mastertrace by using stations with high correlation on that station in station list with the
            # highest correlation coefficient
            stack_result = stack_mastertrace(wfdata_lowf, wfdata_highf, wfdata, picks, params=params[phase_type],
                                             channels=channels_list, method=method, fig_dir=fig_dir)
        else:
            stack_result = load_stacked_trace(eventdir, params[phase_type]['min_corr_stacking'])

        if not stack_result:
            logging.info('No stack result. Continue.')
            continue

        #############################################
        # NOW DO THE FINAL CORRELATION
        # extract stack result
        correlations_dict, nwst_id, trace_master, nstack = stack_result

        if params[phase_type]['plot']:
            # plot correlations of traces used to generate stacked trace
            plot_mastertrace_corr(nwst_id, correlations_dict, wfdata_lowf, stations_dict, picks, trace_master, method,
                                  min_corr=params[phase_type]['min_corr_stacking'], title=eventdir + '_before_stacking',
                                  fig_dir=fig_dir)

        # continue if there is not enough traces for stacking
        if nstack < params[phase_type]['min_stack']:
            logging.warning('Not enough traces to stack: {} (min_stack = {}). Skip this event'.format(nstack,
                                                                                                      params[phase_type][
                                                                                                          'min_stack']))
            continue

        # write unfiltered trace
        trace_master.write(os.path.join(correlation_out_dir, '{}_stacked.mseed'.format(trace_master.id)))

        # now pick stacked trace with PyLoT for a more precise pick (use raw trace, gets filtered by autoPyLoT)
        pick_stacked = repick_mastertrace(wfdata_lowf, trace_master, pylot_parameter, event, event_id, metadata,
                                          phase_type, correlation_out_dir)
        if not pick_stacked:
            continue

        # correlate stations with repicked and stacked master trace
        fig_dir_traces = make_figure_dirs(fig_dir, trace_master.id)

        # filter master trace which was stacked unfiltered
        trace_master_highf = trace_master.copy()
        trace_master_lowf = trace_master.copy()
        if filter_type and filter_options:
            trace_master_lowf.filter(filter_type, **filter_options)
            trace_master_highf.filter(filter_type, **filter_options_final)

        input_list = prepare_correlation_input(wfdata_lowf, picks, channels_list, trace_master_lowf, pick_stacked,
                                               params=params[phase_type], plot=params[phase_type]['plot'],
                                               fig_dir=fig_dir_traces, ncorr=2, wfdata_highf=wfdata_highf,
                                               trace_master_highf=trace_master_highf)

        if (params[phase_type]['plot_detailed'] and not fig_dir_traces) or ncores == 1:
            correlations_dict_stacked = correlate_serial(input_list, plot=params[phase_type]['plot_detailed'])
        else:
            correlations_dict_stacked = correlate_parallel(input_list, ncores=ncores)

        # write output to correlation output directory
        write_correlation_output(correlation_out_dir, correlations_dict, correlations_dict_stacked)

        # export picks to obspy_dmt directories
        export_picks(eventdir=eventdir,
                     correlations_dict=get_stations_min_corr(correlations_dict_stacked,
                                                             params[phase_type]['min_corr_export']),
                     params=params[phase_type], picks=picks, taupypicks=taupypicks_orig, phase_type=true_phase_type,
                     pf_extension=pfe)

        # plot results
        if params[phase_type]['plot']:
            # some random stations for comparison
            coords_master = stations_dict[nwst_id]
            # stations2compare = get_random_stations(coords_master, stations_dict, 60)
            # stations2compare = stations_dict
            stations2compare = get_stations_min_corr(correlations_dict_stacked,
                                                     min_corr=params[phase_type]['min_corr_export'])

            if not stations2compare:
                continue

            fig, axes = plot_stations(stations_dict, stations2compare, coords_master, correlations_dict_stacked,
                                      trace_master, pick_stacked, window_title=eventdir)
            if axes is not None:
                # Note: method changed to 'auto' for exporting picks
                plot_section(wfdata_lowf, trace_master, pick_stacked, picks, stations2compare, channels_list,
                             correlations_dict_stacked, method='auto', axes=axes)
            if fig_dir:
                fig.savefig(os.path.join(fig_dir, 'correlation_result.svg'), dpi=400)
                try:
                    fig_stack = plot_stacked_trace_pick(trace_master, pick_stacked, pylot_parameter)
                    fig_stack.savefig(os.path.join(fig_dir, '{}.svg'.format(trace_master.id)), dpi=400)
                    plt.close(fig_stack)
                except Exception as e:
                    logging.error('Could not plot stacked trace: {}'.format(e))
            else:
                plt.show()

            plt.close(fig)

    if len(correlations_dict_stacked) > 0:
        return True


def remove_outliers(picks, corrected_taupy_picks, threshold):
    ''' delete a pick if difference to corrected taupy pick is larger than threshold'''

    n_picks = len(picks)
    n_outl = 0
    n_no_ref = 0
    for index, pick in list(reversed(list(enumerate(picks)))):
        taup_pick = get_pick4station(corrected_taupy_picks, network_code=pick.waveform_id.network_code,
                                     station_code=pick.waveform_id.station_code, method='taupy')
        if not taup_pick:
            deleted = picks.pop(index)
            logging.debug('Deleted pick as there is no taup reference pick found: {}'.format(deleted.waveform_id))
            n_no_ref += 1
            continue
        diff = pick.time - taup_pick.time
        if abs(diff) > threshold:
            deleted = picks.pop(index)
            logging.debug('Deleted pick outlier ({} seconds): {}'.format(diff, deleted.waveform_id))
            n_outl += 1

    if n_outl:
        logging.info(f'Remove_outliers: Removed {n_outl}/{n_picks} picks with over {threshold} '
                     f'seconds difference to corrected theoretical onset.\n')
    if n_no_ref:
        logging.info(f'Remove_outliers: Removed {n_no_ref} picks because of missing reference theoretical onset.')

    return picks


def remove_invalid_picks(picks):
    ''' Remove picks without uncertainty (invalid PyLoT picks)'''
    count = 0
    deleted_picks_ids = []
    for index, pick in list(reversed(list(enumerate(picks)))):
        uncert = pick.time_errors.uncertainty
        if not uncert:
            deleted = picks.pop(index)
            deleted_picks_ids.append(deleted.waveform_id)
            logging.debug('Removed invalid pick: {}'.format(deleted.waveform_id))
            count += 1
    logging.info('Removed {} invalid picks'.format(count))
    return picks


def get_picks_median(picks):
    return UTCDateTime(int(np.median([pick.time.timestamp for pick in picks if pick.time])))


def get_picks_mean(picks):
    return UTCDateTime(np.mean([pick.time.timestamp for pick in picks if pick.time]))


def get_corrected_taupy_picks(picks, taupypicks, all_available=False):
    ''' get mean/median from picks taupy picks, correct latter for the difference '''

    def nwst_id_from_wfid(wfid):
        return '{}.{}'.format(wfid.network_code if wfid.network_code else '',
                              wfid.station_code if wfid.station_code else '')

    taupypicks = copy.deepcopy(taupypicks)
    # get taupypicks from same stations as there are picks to calculate median/mean of differences
    picks_wf_ids = set([nwst_id_from_wfid(pick.waveform_id) for pick in picks])
    taupypicks_new = [pick for pick in taupypicks if nwst_id_from_wfid(pick.waveform_id) in picks_wf_ids]

    # calculate median and mean of differences
    picks_median = get_picks_median(picks)
    taupy_median = get_picks_median(taupypicks_new)
    median_diff = taupy_median - picks_median

    picks_mean = get_picks_mean(picks)
    taupy_mean = get_picks_mean(taupypicks_new)
    mean_diff = taupy_mean - picks_mean

    logging.info(f'Calculated {len(taupypicks_new)} TauPy theoretical picks.')
    logging.info('Calculated median difference from TauPy theoretical picks of {} seconds. (mean: {})'.format(median_diff,
                                                                                                       mean_diff))

    # return all available taupypicks corrected for median difference to autopicks
    if all_available:
        logging.info('Correcting and returning all available taupy-picks!')
        taupypicks_new = copy.deepcopy(taupypicks)

    for pick in taupypicks_new:
        pick.time -= median_diff

    return taupypicks_new, median_diff


def load_stacked_trace(eventdir, min_corr_stack):
    # load stacked stream (miniseed)
    str_stack_fpaths = glob.glob(os.path.join(eventdir, 'correlation', '*_stacked.mseed'))
    if not len(str_stack_fpaths) == 1:
        logging.warning('No stacked trace for event {} found!'.format(eventdir))
        return
    str_stack_fpath = str_stack_fpaths[0]

    # load dictionary with correlations for stacking (json)
    corr_dict_fpath = os.path.join(eventdir, 'correlation', 'correlations_for_stacking.json')
    if not os.path.isfile(corr_dict_fpath):
        logging.warning('No correlations_for_stacking dict found for event {}!'.format(eventdir))
        return
    corr_dict = json.load(corr_dict_fpath)  # TODO: Check this line

    # get stations for stacking and nstack
    stations4stack = get_stations_min_corr(corr_dict, min_corr_stack)
    nstack = len(stations4stack)

    stacked_stream = read(str_stack_fpath)
    if len(stacked_stream) != 1:
        raise Exception('N_traces!=1 in Stream')
    stacked_trace = stacked_stream[0]
    nwst_id = '{}.{}'.format(stacked_trace.stats.network, stacked_trace.stats.station)

    return corr_dict, nwst_id, stacked_trace, nstack


def repick_mastertrace(wfdata, trace_master, pylot_parameter, event, event_id, metadata, phase_type, corr_out_dir):
    rename_LQT = phase_type == 'S'
    # create an 'artificial' stream object which can be used as input for PyLoT
    stream_master = modify_master_trace4pick(trace_master.copy(), wfdata, rename_LQT=rename_LQT)
    stream_master.write(os.path.join(corr_out_dir, 'stacked_master_trace_pylot.mseed'.format(trace_master.id)))
    try:
        picksdict = autopickstation(stream_master, pylot_parameter, verbose=True, metadata=metadata,
                                    origin=event.origins, iplot=0)
    except Exception as e:
        logging.error('Exception in autopickstation for event {}: {}'.format(event_id, str(e)))
        logging.error(traceback.format_exc())
        return

    # modify input for picks_from_picksdict
    picksdict = {trace_master.stats.station: picksdict[0]}
    # transform to obspy pick
    stacked_picks = picks_from_picksdict(picksdict)
    # only this phase
    stacked_picks = [pick for pick in stacked_picks if pick.phase_hint == phase_type]
    # get pick from stacked picks
    pick_stacked = get_pick4station(stacked_picks, network_code=trace_master.stats.network,
                                    station_code=trace_master.stats.station, method='auto')
    if not pick_stacked:
        raise ValueError('Could not pick stacked trace of event {}. Abort'.format(event_id))

    return pick_stacked


def create_correlation_output_dir(eventdir, fopts, phase_type):
    folder = 'correlation'
    folder = add_fpath_extension(folder, fopts, phase_type)
    export_path = os.path.join(eventdir, folder)
    if not os.path.isdir(export_path):
        os.mkdir(export_path)
    return export_path


def create_correlation_figure_dir(correlation_out_dir):
    export_path = os.path.join(correlation_out_dir, 'figures')
    if not os.path.isdir(export_path):
        os.mkdir(export_path)
    return export_path


def add_fpath_extension(fpath, fopts, phase):
    if fopts:
        freqmin, freqmax = fopts['freqmin'], fopts['freqmax']
        if freqmin:
            fpath += f'_{freqmin}'
        if freqmax:
            fpath += f'-{freqmax}'
    if phase:
        fpath += f'_{phase}'
    return fpath


def write_correlation_output(export_path, correlations_dict, correlations_dict_stacked):
    write_json(correlations_dict, os.path.join(export_path, 'correlations_for_stacking.json'))
    write_json(correlations_dict_stacked, os.path.join(export_path, 'correlation_results.json'))


def modify_master_trace4pick(trace_master, wfdata, rename_LQT=True):
    '''
    Create an artificial Stream with master trace instead of the trace of its corresponding channel.
    This is done to find metadata for correct station metadata which were modified when stacking (e.g. loc=ST)
    USE THIS ONLY FOR PICKING!
    '''
    # fake ZNE coordinates for autopylot
    lqt_zne = {'L': 'Z', 'Q': 'N', 'T': 'E'}

    stream_master = Stream()
    # select stream of master-trace station (includes master trace)
    stream = wfdata.select(network=trace_master.stats.network, station=trace_master.stats.station)
    # remove master trace from stream
    for trace in stream.select(location='99'):
        stream.remove(trace)
    # iterate over all traces (except for master trace with location 99:)
    for trace in stream:
        # if trace should be overwritten with master-trace (same channel)
        if trace.stats.channel == trace_master.stats.channel:
            # take location from old trace and overwrite
            trace_master.stats.location = trace.stats.location
            trace = trace_master
        if rename_LQT:
            channel = trace.stats.channel
            component_new = lqt_zne.get(channel[-1])
            if not component_new:
                logging.warning('Could not rename LQT.')
                continue
            trace.stats.channel = channel[:-1] + component_new
        # append to new stream
        stream_master += trace
    stream_master.trim(starttime=trace_master.stats.starttime, endtime=trace_master.stats.endtime)
    stream_master.normalize()
    return stream_master


def export_picks(eventdir, correlations_dict, picks, taupypicks, params, phase_type, pf_extension):
    threshold = params['export_threshold']
    min_picks_export = params['min_picks_export']
    # make copy so that modified picks are not exported
    new_picks = copy.deepcopy(picks)
    # get current event
    event_id = get_event_id(eventdir)

    event = get_event_pylot(eventdir, extension=pf_extension)
    if not event:
        event = PylotEvent(eventdir)

    # make safety copy of old picks
    fname = get_pickfile_name(event_id, 'pre-correlated')
    fpath_copy = write_event(event, eventdir, fname)
    logging.info('Made safety copy of old pickfile: {}'.format(fpath_copy))

    # iterate over picks and remove picks where correlation failed
    for index, pick in reversed(list(enumerate(new_picks))):
        network = pick.waveform_id.network_code
        station = pick.waveform_id.station_code
        nwst_id = '{}.{}'.format(network, station)
        if nwst_id in correlations_dict:
            dpick = correlations_dict[nwst_id]['dpick']
            ccc = correlations_dict[nwst_id]['ccc']
            uncertainty = correlations_dict[nwst_id]['uncertainty']
            # in case dpick and ccc exist, change pick attributes and go to next pick
            if not np.isnan(dpick) and not np.isnan(ccc) and uncertainty > 0:
                pick.time -= dpick
                # modify resource id and phase hint for pylot
                pick.method_id = ResourceIdentifier('auto')
                pick.phase_hint = phase_type
                # just a little error estimation for testing
                pick.time_errors.uncertainty = uncertainty
                pick.time_errors.lower_uncertainty = uncertainty
                pick.time_errors.upper_uncertainty = uncertainty
                continue
        # in case pick is not in correlations dict or np.nan
        new_picks.pop(index)

    if not new_picks:
        logging.info('No picks left to export! Return.')
        return
    if len(new_picks) < min_picks_export:
        logging.info('Not enough picks to export: {} (min: {}). Return'.format(len(new_picks), min_picks_export))
        return

    fopts = params['filter_options_final']

    # correct taupypicks for new, correlated picks
    taupypicks, time_shift = get_corrected_taupy_picks(new_picks, taupypicks, all_available=True)

    # save taupy picks to file
    write_taupy_picks(event, eventdir, taupypicks, time_shift,
                      extension=add_fpath_extension('corrected_taup_times', fopts, phase_type))

    # remove outliers (more than threshold [s] from corrected taupypicks)
    new_picks = remove_outliers(new_picks, taupypicks, threshold)

    if len(new_picks) < min_picks_export:
        logging.info(
            'Not enough picks to export after removing outliers: {} (min: {}). Return'.format(len(new_picks),
                                                                                              min_picks_export))
        return

    event.picks = new_picks
    fname = get_pickfile_name_correlated(event_id, fopts, phase_type)
    fpath = write_event(event, eventdir, fname)

    logging.info('Wrote {} correlated picks to file {}'.format(len(event.picks), fpath))


def write_taupy_picks(event, eventdir, taupypicks, time_shift, extension='corrected_taup_times'):
    # make copies because both objects are being modified
    event = copy.deepcopy(event)
    taupypicks = copy.deepcopy(taupypicks)

    # get eventname
    eventname = os.path.split(eventdir)[-1]

    # save timeshift to file
    with open(os.path.join(eventdir, 'taup_timeshift.txt'), 'w') as outfile:
        outfile.write('{}\n'.format(time_shift))

    # set method id to 'auto' for PyLoT
    for pick in taupypicks:
        pick.method_id = 'auto'
        # artificial uncertainty
        pick.time_errors.uncertainty = 0.0001
    event.picks = taupypicks
    fname = get_pickfile_name(eventname, extension)
    write_event(event, eventdir, fname)


def write_event(event, eventdir, fname):
    fpath = os.path.join(eventdir, fname)
    event.write(fpath, format='QUAKEML')
    return fpath


def get_pickfile_name(event_id, fname_extension):
    return 'PyLoT_{}_{}.xml'.format(event_id, fname_extension)


def get_pickfile_name_correlated(event_id, fopts, phase_type):
    fname_extension = add_fpath_extension('correlated', fopts, phase_type)
    return get_pickfile_name(event_id, fname_extension)


# def prepare_trace(wfdata, network_code, station_code, channels, freq):
#     # get trace of this station
#     for channel in channels:
#         st_this = wfdata.select(network=network_code, station=station_code, channel=channel)
#         if st_this:
#             break
#     if not st_this:
#         return
#     trace_this = st_this[0]
#
#     # resample and filter trace
#     trace_this.resample(freq)
#     trace_this.filter('bandpass', freqmin=0.03, freqmax=0.5)
#     return trace_this


def prepare_correlation_input(wfdata, picks, channels, trace_master, pick_master, params, plot=None, fig_dir=None,
                              ncorr=1, wfdata_highf=None, trace_master_highf=None):
    # prepare input for multiprocessing worker for all 'other' picks to correlate with current master-trace
    input_list = []

    assert (ncorr in [1, 2]), 'ncorr has to be 1 or 2'

    for pick_other in picks:
        trace_other_highf = None
        network_other = pick_other.waveform_id.network_code
        station_other = pick_other.waveform_id.station_code
        nwst_id_other = '{nw}.{st}'.format(nw=network_other, st=station_other)
        if not pick_other.time:
            logging.debug('No pick for station {nwst}!'.format(nwst=nwst_id_other))
            continue
        for channel in channels:
            stream_other = wfdata.select(network=network_other, station=station_other, channel=channel)
            if ncorr == 2:
                stream_other_highf = wfdata_highf.select(network=network_other, station=station_other, channel=channel)
            if stream_other:
                break
        if not stream_other:
            logging.debug('No waveform data for station {nwst}!'.format(nwst=nwst_id_other))
            continue
        trace_other = stream_other[0]
        if ncorr == 2:
            trace_other_highf = stream_other_highf[0]
        if trace_other == stream_other:
            continue

        input_dict = {'nwst_id': nwst_id_other, 'trace1': trace_master, 'pick1': pick_master, 'trace2': trace_other,
                      'pick2': pick_other, 'channel': channel, 't_before': params['t_before'],
                      't_after': params['t_after'], 'cc_maxlag': params['cc_maxlag'],
                      'cc_maxlag2': params['cc_maxlag2'], 'plot': plot, 'fig_dir': fig_dir, 'ncorr': ncorr,
                      'trace1_highf': trace_master_highf, 'trace2_highf': trace_other_highf,
                      'min_corr': params['min_corr_export']}
        input_list.append(input_dict)

    return input_list


def stack_mastertrace(wfdata_lowf, wfdata_highf, wfdata_raw, picks, params, channels, method, fig_dir):
    '''
    Correlate all stations with the first available station given in station_list, a list containing central,
    permanent, long operating, low noise stations with descending priority.
    A master trace will be created by stacking well correlating traces onto this station.
    '''

    def get_best_station4stack(station_results):
        ''' return station with maximum mean_ccc'''
        ccc_means = {nwst_id: value['mean_ccc'] for nwst_id, value in station_results.items() if
                     not np.isnan(value['mean_ccc'])}
        if len(ccc_means) < 1:
            logging.warning('No valid station found for stacking! Return.')
            return
        best_station_id = max(ccc_means, key=ccc_means.get)
        logging.info('Found highest mean correlation for station {} ({})'.format(best_station_id, max(ccc_means.values())))
        return best_station_id

    station_results = iterate_correlation(wfdata_lowf, wfdata_highf, channels, picks, method, params, fig_dir=fig_dir)
    nwst_id_master = get_best_station4stack(station_results)

    # in case no stream with a valid pick is found
    if not nwst_id_master:
        logging.info('No mastertrace found! Will skip this event.')
        return False

    trace_master = station_results[nwst_id_master]['trace']
    stations4stack = station_results[nwst_id_master]['stations4stack']
    correlations_dict = station_results[nwst_id_master]['correlations_dict']

    wfdata_highf += trace_master

    dt_pre, dt_post = params['dt_stacking']
    trace_master, nstack = apply_stacking(trace_master, stations4stack, wfdata_raw, picks, method=method,
                                          check_RMS=params['check_RMS'], plot=params['plot'], fig_dir=fig_dir,
                                          dt_pre=dt_pre, dt_post=dt_post)

    return correlations_dict, nwst_id_master, trace_master, nstack


def iterate_correlation(wfdata_lowf, wfdata_highf, channels, picks, method, params, fig_dir=None):
    '''iterate over possible stations for master-trace store them and return a dictionary'''

    station_results = {nwst_id: {'mean_ccc': np.nan, 'correlations_dict': None, 'stations4stack': None, 'trace': None}
                       for nwst_id in params['station_list']}

    for nwst_id_master in params['station_list']:
        logging.info(20 * '#')
        logging.info('Starting correlation for station: {}'.format(nwst_id_master))
        nw, st = nwst_id_master.split('.')
        for channel in channels:
            stream_master_lowf = wfdata_lowf.select(network=nw, station=st, channel=channel)
            stream_master_highf = wfdata_highf.select(network=nw, station=st, channel=channel)
            # TODO: Account for Q and T
            if stream_master_lowf and stream_master_highf:
                break
        if not stream_master_lowf and not stream_master_highf:
            continue
        pick_master = get_pick4station(picks, nw, st, method)
        if not pick_master or not pick_master.time:
            continue

        # prepare master-trace
        trace_master_lowf = stream_master_lowf[0]
        trace_master_highf = stream_master_highf[0]

        # make new trace with artificial location ID
        trace_master_lowf = trace_master_lowf.copy()
        trace_master_highf = trace_master_highf.copy()
        trace_master_highf.stats.location = '99'

        # trace_master.resample(freq)
        # trace_master.filter(filter_type, **filter_options)

        fig_dir_traces = ''
        if fig_dir and os.path.isdir(fig_dir):
            fig_dir_traces = os.path.join(fig_dir, 'corr_with_{}'.format(trace_master_highf.id))
            if not os.path.isdir(fig_dir_traces):
                os.mkdir(fig_dir_traces)

        input_list = prepare_correlation_input(wfdata_lowf, picks, channels, trace_master_lowf, pick_master,
                                               trace_master_highf=trace_master_highf, wfdata_highf=wfdata_highf,
                                               params=params, plot=params['plot_detailed'], fig_dir=fig_dir_traces,
                                               ncorr=2)

        if params['plot_detailed'] and not fig_dir_traces:
            # serial
            correlations_dict = correlate_serial(input_list, plot=True)
        else:
            # parallel
            correlations_dict = correlate_parallel(input_list, ncores=params['ncores'])

        stations4stack = get_stations_min_corr(correlations_dict, min_corr=params['min_corr_stacking'])

        station_results[nwst_id_master]['mean_ccc'] = np.nanmean([item['ccc'] for item in correlations_dict.values()])
        station_results[nwst_id_master]['correlations_dict'] = correlations_dict
        station_results[nwst_id_master]['stations4stack'] = stations4stack
        station_results[nwst_id_master]['trace'] = trace_master_lowf

    return station_results


def apply_stacking(trace_master, stations4stack, wfdata, picks, method, check_RMS, dt_pre=250., dt_post=250.,
                   plot=False, fig_dir=None):
    def check_trace_length_and_cut(trace, correlated_midtime, dt_pre, dt_post):
        starttime = correlated_midtime - dt_pre
        endtime = correlated_midtime + dt_post

        if trace.stats.starttime > starttime or trace.stats.endtime < endtime:
            raise ValueError('Trace too short for safe cutting. Would create inconsistent stacking. Abort!')

        trace.trim(starttime=starttime, endtime=endtime)

    traces4stack = []

    pick = get_pick4station(picks, trace_master.stats.network, trace_master.stats.station, method)

    trace_master = trace_master.copy()
    trace_master.stats.location = 'ST'
    check_trace_length_and_cut(trace_master, pick.time, dt_pre=dt_pre, dt_post=dt_post)

    # empty trace so that it does not stack twice
    trace_master.data = np.zeros(len(trace_master.data))

    count = 0
    for nwst_id, corr_result in stations4stack.items():
        nw, st = nwst_id.split('.')
        pick_other = get_pick4station(picks, nw, st, method)
        dpick = corr_result['dpick']
        channel = corr_result['channel']
        stream_other = wfdata.select(network=nw, station=st, channel=channel)
        correlated_midtime = pick_other.time - dpick
        trace_other = stream_other[0].copy()
        check_trace_length_and_cut(trace_other, correlated_midtime, dt_pre=dt_pre, dt_post=dt_post)

        if not len(trace_other) == len(trace_master):
            logging.warning('Can not stack trace on master trace because of different lengths: {}. Continue'.format(nwst_id))
            continue

        traces4stack.append(trace_other)

    if check_RMS:
        traces4stack = checkRMS(traces4stack, plot=plot, fig_dir=fig_dir)

    if plot:
        fig = plt.figure(figsize=(16, 9))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

    for trace_other in traces4stack:
        # stack on mastertrace
        trace_normalized = trace_other.copy().normalize()
        trace_master.data += trace_normalized.data
        count += 1
        if plot:
            ax1.plot(trace_master.data, label=count)
            ax2.plot(trace_normalized.data, color='k', alpha=0.1)

    if plot:
        if fig_dir:
            fig.savefig(os.path.join(fig_dir, 'traces_stacked.svg'), dpi=300.)
        else:
            fig.show()

        plt.close(fig)

    logging.info('Successfully stacked {} traces onto master trace.'.format(count))

    return trace_master, count


def checkRMS(traces4stack, plot=False, fig_dir=None, max_it=10, ntimes_std=5.):
    rms_list = []
    trace_names = []

    traces4stack = sorted(traces4stack, key=lambda x: x.id)

    for trace in traces4stack:
        rms_list.append(RMS(trace.data))
        trace_names.append(trace.id)

    # iterative elimination of RMS outliers
    iterate = True
    count = 0
    while iterate:
        count += 1
        if count >= max_it:
            logging.debug('Maximum iterations ({}) of checkRMS function reached.'.format(max_it))
            break
        std = np.std(rms_list)
        mean = np.mean(rms_list)

        iterate = False
        for index, tr_rms in list(reversed(list(enumerate(zip(traces4stack, rms_list))))):
            trace, rms = tr_rms
            if not mean - ntimes_std * std <= rms <= mean + ntimes_std * std:
                logging.debug('Removing trace {}.{} because of strong RMS deviation.'.format(trace.stats.network,
                                                                                             trace.stats.station))
                if plot:
                    plot_rms(rms_list, trace_names, mean, std, fig_dir, count, ntimes_std)

                traces4stack.remove(trace)
                rms_list.pop(index)
                trace_names.pop(index)
                iterate = True

    if plot:
        plot_rms(rms_list, trace_names, mean, std, fig_dir, count, ntimes_std)

    return traces4stack


def plot_rms(rms_list, trace_names, mean, std, fig_dir, count, ntimes_std):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.semilogy(rms_list, 'b.')
    ax.axhline(mean, color='k', linestyle='-')
    ax.axhline(mean + ntimes_std * std, color='k', linestyle='--')
    ax.axhline(mean - ntimes_std * std, color='k', linestyle='--')
    ax.set_ylabel('RMS [m/s]')
    ax.set_xticks(list(range(len(trace_names))))
    ax.set_xticklabels(trace_names, rotation=90, fontsize=2)
    if fig_dir:
        fig.savefig(os.path.join(fig_dir, 'RMS_check_{:02}.svg'.format(count)), dpi=400.)
    else:
        plt.show()

    plt.close(fig)


def RMS(X):
    """
    Returns root mean square of a given array LON
    """
    return np.sqrt(np.sum(np.power(X, 2)) / len(X))


def resample_parallel(stream, freq, ncores):
    input_list = [{'trace': trace, 'freq': freq} for trace in stream.traces]

    parallel = Parallel(n_jobs=ncores)

    logging.info('Resample_parallel: Generated {} parallel jobs.'.format(ncores))
    output_list = parallel(delayed(resample_worker)(item) for item in input_list)

    logging.info('Parallel execution finished.')

    stream.traces = output_list
    return stream


def rotate_stream(stream, metadata, origin, stations_dict, channels, inclination, ncores):
    ''' Rotate stream to LQT. To do this all traces have to be cut to equal lenths'''
    input_list = []
    cut_stream_to_same_length(stream)
    new_stream = Stream()
    # nwst_ids = np.unique([f'{tr.stats.network}.{tr.stats.station}' for tr in stream])
    # for nwst_id in nwst_ids:
    #     if not nwst_id in stations_dict.keys():
    #         logging.info(f'{nwst_id} not found')
    for nwst_id, rec_coords in stations_dict.items():
        if not rec_coords: continue
        nw, sta = nwst_id.split('.')
        st = stream.select(network=nw, station=sta)
        if not st:
            logging.debug(f'No stream for {nwst_id}')
            continue
        locations = np.unique([tr.stats.location for tr in st])
        # TODO: Take other locations into account? Excluded here due to possible ambiguity with usage of nwst_id
        for location in locations:
            st = stream.select(network=nw, station=sta, location=location)
            dic = {'stream': st, 'src_lat': origin.get('latitude'), 'src_lon': origin.get('longitude'),
                   'rec_lat': rec_coords.get('latitude'), 'rec_lon': rec_coords.get('longitude'),
                   'inclination': inclination}

            rotation_worker(dic)
            new_stream += st  # input_list.append(dic)

    logging.info(50 * '#')
    logging.info(f'LENGTH NEW STREAM: {len(new_stream)}, LENGTH OLD STREAM: {len(stream)}')
    logging.info(50 * '#')

    return new_stream

    pool = multiprocessing.Pool(ncores, maxtasksperchild=100)
    logging.info('Resample_parallel: Generated multiprocessing pool with {} cores.'.format(ncores))
    output_list = pool.map(rotation_worker, input_list, chunksize=10)
    pool.close()
    logging.info('Closed multiprocessing pool.')
    pool.join()
    del (pool)
    stream.traces = [tr for tr in output_list if tr is not None]
    return stream


def correlate_parallel(input_list, ncores):
    parallel = Parallel(n_jobs=ncores)

    logging.info('Correlate_parallel: Generated {} parallel jobs.'.format(ncores))
    correlation_result = parallel(delayed(correlation_worker)(item) for item in input_list)

    logging.info('Parallel execution finished.')

    return unpack_result(correlation_result)


def correlate_serial(input_list, plot=False):
    correlation_result = []
    for input_dict in input_list:
        input_dict['plot'] = plot
        correlation_result.append(correlation_worker(input_dict))
    return unpack_result(correlation_result)


def taupy_parallel(input_list, ncores):
    parallel = Parallel(n_jobs=ncores)

    logging.info('Taupy_parallel: Generated {} parallel jobs.'.format(ncores))
    taupy_results = parallel(delayed(taupy_worker)(item) for item in input_list)

    logging.info('Parallel execution finished.')

    return unpack_result(taupy_results)


def unpack_result(result):
    result_dict = {item['nwst_id']: {key: item[key] for key in item.keys()} for item in result}
    nerr = 0
    for item in result_dict.values():
        if item['err']:
            logging.debug(item['err'])
            nerr += 1
    logging.info('Unpack results: Found {} errors after multiprocessing'.format(nerr))
    return result_dict


def get_taupy_picks(origin, stations_dict, pylot_parameter, phase_type, ncores=None):
    input_list = []

    taup_phases = pylot_parameter['taup_phases'].split(',')
    taup_model = pylot_parameter['taup_model']
    taup_phases = [phase for phase in taup_phases if identifyPhaseID(phase) == phase_type]

    model = TauPyModel(taup_model)

    src_depth = origin.depth
    src_lat = origin.latitude
    src_lon = origin.longitude
    for nwst_id, coords in stations_dict.items():
        rec_lat = coords['latitude']
        rec_lon = coords['longitude']
        taupy_input = {'source_depth_in_km': src_depth, 'source_latitude_in_deg': src_lat,
                       'source_longitude_in_deg': src_lon, 'receiver_latitude_in_deg': rec_lat,
                       'receiver_longitude_in_deg': rec_lon, 'phase_list': taup_phases, }
        input_dict = {'nwst_id': nwst_id, 'taupy_input': taupy_input, 'model': model, 'source_time': origin.time, }
        input_list.append(input_dict)

    taupy_result = taupy_parallel(input_list, ncores=ncores)
    check_taupy_phases(taupy_result)

    taupy_picks = create_artificial_picks(taupy_result)
    return taupy_picks


def create_artificial_picks(taupy_result):
    artificial_picks = []
    for nwst_id, taupy_dict in taupy_result.items():
        nw, st = nwst_id.split('.')
        waveform_id = WaveformStreamID(network_code=nw, station_code=st)
        method = ResourceIdentifier('taupy')
        pick = Pick(force_resource_id=True, waveform_id=waveform_id, time=taupy_dict['phase_time'], method_id=method,
                    phase_hint=taupy_dict['phase_name'])
        artificial_picks.append(pick)
    return artificial_picks


def check_taupy_phases(taupy_result):
    test_phase_name = list(taupy_result.values())[0]['phase_name']
    phase_names_equal = [item['phase_name'] == test_phase_name for item in taupy_result.values()]
    if not all(phase_names_equal):
        logging.info('Different first arriving phases detected in TauPy phases for this event.')


def taupy_worker(input_dict):
    model = input_dict['model']
    taupy_input = input_dict['taupy_input']
    source_time = input_dict['source_time']

    try:
        arrivals = model.get_travel_times_geo(**taupy_input)
        first_arrival = arrivals[0]
        output_dict = dict(nwst_id=input_dict['nwst_id'], phase_name=first_arrival.name,
                           phase_time=source_time + first_arrival.time, phase_dist=first_arrival.distance, err=None, )
    except Exception as e:
        output_dict = dict(nwst_id=input_dict['nwst_id'], phase_name=None, phase_time=None, err=str(e))
    return output_dict


def resample_worker(input_dict):
    trace = input_dict['trace']
    freq = input_dict['freq']
    if freq == trace.stats.sampling_rate:
        return trace
    else:
        return trace.resample(freq, no_filter=False)


def rotation_worker(input_dict):
    stream = input_dict['stream']
    tstart = max([tr.stats.starttime for tr in stream])
    tend = min([tr.stats.endtime for tr in stream])
    if not tstart or not tend:
        logging.debug('Could not cut stream {} for rotation'.format(stream[0].id))
        return
    stream.trim(tstart, tend, nearest_sample=False)
    dist, az, baz = gps2dist_azimuth(input_dict.get('src_lat'), input_dict.get('src_lon'), input_dict.get('rec_lat'),
                                     input_dict.get('rec_lon'))
    try:
        stream.rotate('ZNE->LQT', back_azimuth=baz, inclination=input_dict.get('inclination'))
    except ValueError as e:
        logging.error(f'Could not rotate Stream to LQT: {e}')
    return stream


def correlation_worker(input_dict):
    '''worker function for multiprocessing'''

    # unpack input dictionary
    nwst_id = input_dict['nwst_id']
    pick_this = input_dict['pick1']
    other_pick = input_dict['pick2']
    channel = input_dict['channel']
    plot = input_dict['plot']
    fig_dir = input_dict['fig_dir']
    t_before = input_dict['t_before']
    t_after = input_dict['t_after']
    cc_maxlag = input_dict['cc_maxlag']
    cc_maxlag2 = input_dict['cc_maxlag2']

    # apply correlation and pick correction
    # Note: Do correlation twice, once to get rough time shift (lower freq),
    # second time to get more precise error for shifted pick.
    # Therefore, traces have to be cut to same length around pick and correlated again

    try:
        logging.debug(f'Starting Pick correction for {nwst_id}')
        xcpc = Xcorr_pick_correction(pick_this.time, input_dict['trace1'], other_pick.time, input_dict['trace2'],
                                     t_before=t_before, t_after=t_after, cc_maxlag=cc_maxlag)

        dpick, ccc, uncert, fwm = xcpc.cross_correlation(plot, fig_dir, plot_name='dpick')
        logging.debug(f'dpick of first correlation: {dpick}')

        if input_dict['ncorr'] > 1:  # and not ccc <= 0:
            xcpc2 = Xcorr_pick_correction(pick_this.time, input_dict['trace1_highf'], other_pick.time + dpick,
                                          input_dict['trace2_highf'], t_before=1., t_after=40., cc_maxlag=cc_maxlag2)

            dpick2, ccc, uncert, fwm = xcpc2.cross_correlation(plot=plot, fig_dir=fig_dir, plot_name='error',
                                                               min_corr=input_dict['min_corr'])

            logging.debug(f'dpick of second correlation: {dpick2}')

            # if there is a new shift from the second correlation, apply this shift as well (?)
            dpick += dpick2

        err = None
    except Exception as e:
        dpick = np.nan
        ccc = np.nan
        err = str(e)
        fwm = np.nan
        uncert = np.nan

    if err != None:
        err = str(err)

    return {'nwst_id': nwst_id, 'dpick': -1 * dpick, 'ccc': ccc, 'fwm': fwm, 'uncertainty': uncert, 'err': err,
            'channel': channel, }


def get_pick4station(picks, network_code, station_code, method='auto'):
    for pick in picks:
        if pick.waveform_id.network_code == network_code:
            if pick.waveform_id.station_code == station_code:
                if pick.method_id.id.endswith(method):
                    return pick


def get_stations_min_corr(corr_results, min_corr):
    corr_results = {nwst_id: result for nwst_id, result in corr_results.items() if result['ccc'] > min_corr}
    return corr_results


def plot_mastertrace_corr(nwst_id, corr_results, wfdata, stations_dict, picks, trace_master, method, min_corr, title,
                          fig_dir=None):
    nw, st = nwst_id.split('.')
    coords_master = stations_dict[nwst_id]
    pick_master = get_pick4station(picks, nw, st, method)
    stations2plot = get_stations_min_corr(corr_results, min_corr=min_corr)
    # wfdata.select(network=trace_master.stats.network, station=trace_master.stats.station,
    #              channel=trace_master.stats.channel).plot(automerge=False)
    fig, axes = plot_stations(stations_dict, stations2plot, coords_master, corr_results, trace_master, pick_master,
                              window_title=title)
    if fig_dir:
        fig.savefig(os.path.join(fig_dir, 'correlation_all4stack.svg'), dpi=300.)
    fig, axes = plot_stations(stations_dict, stations_dict, coords_master, corr_results, trace_master, pick_master,
                              window_title=title)
    if fig_dir:
        fig.savefig(os.path.join(fig_dir, 'correlation_4stack_threshold.svg'), dpi=300.)

    # stacked_trace = stack_traces(wfdata, trace_master, pick_master, picks, stations_dict, trace_master.stats.channel,
    #                             corr_results, plot_section=True, axes=axes)

    if not fig_dir:
        plt.show()


def make_figure_dirs(fig_dir, trace_master_id):
    fig_dir_traces = ''
    if fig_dir and os.path.isdir(fig_dir):
        fig_dir_traces = os.path.join(fig_dir, 'corr_with_{}'.format(trace_master_id))
        if not os.path.isdir(fig_dir_traces):
            os.mkdir(fig_dir_traces)
    return fig_dir_traces


def plot_section(wfdata, trace_this, pick_this, picks, stations2compare, channels, corr_results, method, dt_pre=20.,
                 dt_post=50, axes=None, max_stations=50.):
    '''Plot a section with all stations used for correlation on ax'''
    ylabels = []
    yticks = []

    trace_this = trace_this.copy()
    trace_this.trim(starttime=pick_this.time - dt_pre, endtime=pick_this.time + dt_post)

    # iterate over all closest stations ("other" stations)
    for index, nwst_id in enumerate(stations2compare.keys()):
        if index >= max_stations:
            break
        nw, st = nwst_id.split('.')

        if not nwst_id in corr_results:
            logging.info('No correlation result for station {}'.format(nwst_id))
            continue

        # get correlation results form dictionary
        ccc = corr_results[nwst_id]['ccc']
        dpick = corr_results[nwst_id]['dpick']

        # station names as y labels
        yticks.append(index)
        ylabels.append('{id}({ccc:1.2}/{dpick:1.2} s)'.format(id=nwst_id, ccc=ccc, dpick=dpick))

        ax_sec = axes['section']

        # plot reference trace (trace_this) with pick
        ax_sec.plot(0.5 * trace_this.data / max(abs(trace_this.data)) + index, 'k--', alpha=0.7, linewidth=0.5)
        ax_sec.vlines((pick_this.time - trace_this.stats.starttime) / trace_this.stats.delta, index - 0.5, index + 0.5,
                      color='k', linestyles='--', alpha=0.7, lw=0.5)

        # get pick of other station
        pick_other = get_pick4station(picks, station_code=st, network_code=nw, method=method)

        # continue if no pick or no pick difference given for other station
        if np.isnan(dpick) or not pick_other or not pick_other.time:
            continue

        # continue if there are no data for station for whatever reason
        for channel in channels:
            stream = wfdata.select(station=st, network=nw, channel=channel)
            if stream:
                break
        if not stream:
            continue

        # plot all stations around correlated time (shift traces onto reference trace)
        correlated_midtime = pick_other.time - dpick
        trace_other = stream[0].copy()
        trace_other.trim(starttime=correlated_midtime - dt_pre, endtime=correlated_midtime + dt_post)

        ax_sec.plot(0.5 * trace_other.data / max(abs(trace_other.data)) + index, color='k', linewidth=0.5)
        ax_sec.vlines((pick_other.time - trace_other.stats.starttime) / trace_other.stats.delta, index - 0.5,
                      index + 0.5, color='r', lw=0.5)

    # Plot desciption
    ax_sec.set_yticks(yticks)
    ax_sec.set_yticklabels(ylabels)
    ax_sec.set_title('Section with corresponding picks.')
    ax_sec.set_xlabel('Samples. Traces are shifted in time.')


def plot_stacked_trace_pick(trace, pick, pylot_parameter):
    trace_filt = trace.copy()
    if pylot_parameter['filter_type'] and pylot_parameter['filter_options']:
        ftype = pylot_parameter['filter_type'][0]
        forder = pylot_parameter['filter_order'][0]
        freqmin, freqmax = pylot_parameter['bpz2']
        trace_filt.detrend('linear')
        trace_filt.taper(0.05, type='cosine')
        trace_filt.filter(type=ftype, corners=forder, freqmin=freqmin, freqmax=freqmax)
    fig = plt.figure(figsize=(16, 9))
    ax_unf = fig.add_subplot(211)
    ax_filt = fig.add_subplot(212)
    time_ax = np.linspace(0, trace.stats.endtime - trace.stats.starttime, trace.stats.npts)
    ax_unf.plot(time_ax, trace.data)
    ax_filt.plot(time_ax, trace_filt.data)
    ax_filt.axvline((pick.time - trace.stats.starttime), c='r', linestyle='-')
    # ax.axvline((pick.time + mean_dpick - trace.stats.starttime), c='r')
    ax_unf.set_title(trace.id)
    ax_filt.set_xlabel('Seconds since {}'.format(trace.stats.starttime))
    return fig


def plot_stations(stations_dict, stations2compare, coords_this, corr_result, trace, pick, window_title=None):
    ''' Plot routine to check proximity algorithm. '''

    title = trace.id

    # create array of longitudes (first column) and latitude (second column) for all stations (row)
    coords_all = np.array([[station['longitude'], station['latitude']] for station in stations_dict.values()])
    # get coords of stations in vicinity of "this" station
    coords_selected = np.array([[stations_dict[station]['longitude'], stations_dict[station]['latitude']] for station in
                                stations2compare.keys()])

    # create arrays with cross correlation coefficient and cc-pick difference as color representation in scatter
    cccs = np.array([corr_result[station]['ccc'] if station in corr_result.keys() else np.nan for station in
                     stations2compare.keys()])
    dpicks = np.array([corr_result[station]['dpick'] if station in corr_result.keys() else np.nan for station in
                       stations2compare.keys()])

    mean_dpick = np.nanmean(dpicks)

    title += ' | mean ccc: {:1.3}'.format(np.nanmean(cccs))
    title += ' | mean dpick: {:2.3} [s]'.format(mean_dpick)

    try:
        abs_max = np.nanmax(np.abs(np.array(dpicks)))
    except ValueError as e:
        logging.error(str(e))
        return

    gs = plt.GridSpec(3, 3, height_ratios=[1, 1, 4], width_ratios=[1, 1, 1])
    fig = plt.figure(figsize=(16, 9))
    fig.canvas.manager.set_window_title(window_title)
    ax_sec = fig.add_subplot(gs[0:, 0])
    ax_trace = fig.add_subplot(gs[0, 1:])
    ax_cc = fig.add_subplot(gs[1, 1])
    ax_dp = fig.add_subplot(gs[1, 2])
    ax_map = fig.add_subplot(gs[2, 1:])

    axes = {'section': ax_sec, 'trace': ax_trace, 'ccc': ax_cc, 'dpick': ax_dp, 'map': ax_map}

    trace = trace.copy()
    trace.trim(starttime=pick.time - 200, endtime=pick.time + 200)

    time_ax = np.linspace(0, trace.stats.endtime - trace.stats.starttime, trace.stats.npts)
    ax_trace.plot(time_ax, trace.data)
    ax_trace.axvline((pick.time - trace.stats.starttime), c='k', linestyle='--')
    ax_trace.axvline((pick.time + mean_dpick - trace.stats.starttime), c='r')
    ax_trace.set_title(title)
    ax_trace.set_xlabel('Seconds since {}'.format(trace.stats.starttime))

    ax_map.scatter(coords_all[:, 0], coords_all[:, 1], c='k', s=3)
    ax_map.set_xlabel('Longitude [$^\circ$]')
    ax_map.set_ylabel('Latitude [$^\circ$]')

    sc = ax_map.scatter(coords_selected[:, 0], coords_selected[:, 1], c=dpicks, s=40 * cccs,
                        cmap=plt.get_cmap('seismic'), edgecolors='grey', vmin=-abs_max, vmax=abs_max, linewidths=0.5)
    # sc.set_facecolor('none')
    plt.colorbar(sc, label='Color codes dpicks [s]. Size codes cross correlation coefficient [0-1]')
    ax_map.scatter(coords_this['longitude'], coords_this['latitude'], marker='*', s=40, c='g')
    ax_cc.hist(cccs, bins=50)
    ax_cc.set_xlabel('Cross Correlation Coefficient')
    ax_dp.hist(dpicks, bins=50)
    ax_dp.set_xlabel('Pick difference distribution [s]')

    return fig, axes


def get_data(eventdir, data_dir, headonly=False):
    ''' Read and return waveformdata read from eventdir/data_dir. '''
    wfdata_path = os.path.join(eventdir, data_dir)
    wfdata = read(os.path.join(wfdata_path, '*'), headonly=headonly)
    return wfdata


def get_closest_stations(coords, stations_dict, n):
    ''' Calculate distances and return n closest stations in stations_dict for station at coords. '''
    distances = {}
    for station_id, st_coords in stations_dict.items():
        dist = gps2dist_azimuth(coords['latitude'], coords['longitude'], st_coords['latitude'], st_coords['longitude'],
                                a=6.371e6, f=0)[0]
        # exclude same coordinate (self or other instrument at same location)
        if dist == 0:
            continue
        distances[dist] = station_id

    closest_distances = sorted(list(distances.keys()))[:n + 1]
    closest_stations = {station: dist for dist, station in distances.items() if dist in closest_distances}
    return closest_stations


def get_random_stations(coords, stations_dict, n):
    ''' Calculate distances and return n randomly selected stations in stations_dict for station at coords. '''
    random_keys = random.sample(list(stations_dict.keys()), n)
    distances = {}
    for station_id, st_coords in stations_dict.items():
        dist = gps2dist_azimuth(coords['latitude'], coords['longitude'], st_coords['latitude'], st_coords['longitude'],
                                a=6.371e6, f=0)[0]
        # exclude same coordinate (self or other instrument at same location)
        if dist == 0:
            continue
        distances[dist] = station_id

    random_stations = {station: dist for dist, station in distances.items() if station in random_keys}
    return random_stations


# Notes: if use_master_trace is set True, traces above a specified ccc threshold will be stacked onto one 'ideal'
# station into a master-trace. This trace will be picked and used again for cross correlation with all other stations
# to find a pick for these stations.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Correlate picks from PyLoT.')
    parser.add_argument('dmt_path', help='path containing dmt_database with PyLoT picks')
    parser.add_argument('pylot_infile', help='path to autoPyLoT inputfile (pylot.in)')
    parser.add_argument('--params', help='path to correlation input parameter file (parameters.yaml)',
                        default='parameters.yaml')
    parser.add_argument('-p', dest='plot', action='store_true', default=False, help='Generate plots')
    parser.add_argument('--show', dest='show_fig', action='store_true', default=False,
                        help='Show plots instead of saving them to directory')
    parser.add_argument('--no_rms_check', dest='no_rms_check', action='store_true', default=False,
                        help='Do not check for RMS outliers.')
    parser.add_argument('--reuse_stacked', dest='use_stacked_trace', action='store_true', default=False,
                        help='Re-use stacked trace of prior correlation. Make sure there is only one stacked trace in'
                             'correlation directory!!!')
    parser.add_argument('--data_dir', default='processed', help='Data subdirectory (e.g. processed or raw)')
    parser.add_argument('--event_blacklist', dest='blacklist', default=None,
                        help='Skip event_ids specified in textfile.')
    parser.add_argument('-pd', dest='plot_detailed', action='store_true', default=False, help='Generate detailed plots')
    parser.add_argument('-u', dest='update', action='store_true', default=False,
                        help='Update picks. Skip events if pickfile exists.')
    parser.add_argument('-pe', dest='pickfile_ext', default='_autopylot',
                        help='Pickfile extension (if taupy is not used), default: _autopylot')
    parser.add_argument('-n', dest='ncores', default=None, help='number of cores for multiprocessing')
    parser.add_argument('-istart', default=0, help='first event index')
    parser.add_argument('-istop', default=1e9, help='last event index')

    args = parser.parse_args()

    # PARAMETER DEFINITIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # plot options
    plot_dict = dict(plot=args.plot, plot_detailed=args.plot_detailed, save_fig=not (args.show_fig))

    # alparray configuration: ZNE ~ 123
    channels = {'P': ['*Z', '*1'], 'S': ['*Q', '*T']}  # '*N', '*E', '*2', '*3',
    debug_levels = {'debug': logging.DEBUG,
                    'info': logging.INFO,
                    'warn': logging.WARNING,
                    'error': logging.ERROR,
                    'critical': logging.CRITICAL}


    if args.ncores is not None:
        _ncores = int(args.ncores)
    else:
        _ncores = None

    with open(args.params) as infile:
        parameters_dict = yaml.safe_load(infile)

    params = {phase: None for phase in parameters_dict['pick_phases']}

    logger = logging.getLogger()
    logger.setLevel(debug_levels.get(parameters_dict['logging']))
    handler = logging.StreamHandler(sys.stdout)
    fhandler = logging.FileHandler('pick_corr.out')
    logger.addHandler(handler)
    logger.addHandler(fhandler)

    for phase in params.keys():
        params_phase = CorrelationParameters(ncores=_ncores)
        params_phase.add_parameters(**plot_dict)
        params_phase.add_parameters(**parameters_dict[phase])
        params[phase] = params_phase

    correlation_main(args.dmt_path, args.pylot_infile, params=params, istart=int(args.istart), istop=int(args.istop),
                     channel_config=channels, update=args.update, event_blacklist=args.blacklist)
