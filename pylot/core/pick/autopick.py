#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Function to run automated picking algorithms using AIC,
HOS and AR prediction. Uses objects CharFuns and Picker and
function conglomerate utils.

:author: MAGS2 EP3 working group / Ludger Kueperkoch
"""
import copy
import traceback

import matplotlib.pyplot as plt
import numpy as np
from obspy import Trace
from obspy.taup import TauPyModel

from pylot.core.pick.charfuns import CharacteristicFunction
from pylot.core.pick.charfuns import HOScf, AICcf, ARZcf, ARHcf, AR3Ccf
from pylot.core.pick.picker import AICPicker, PragPicker
from pylot.core.pick.utils import checksignallength, checkZ4S, earllatepicker, \
    getSNR, fmpicker, checkPonsets, wadaticheck, get_quality_class, PickingFailedException, MissingTraceException
from pylot.core.util.utils import getPatternLine, gen_Pool, \
    get_bool, identifyPhaseID, get_None, correct_iplot


def autopickevent(data, param, iplot=0, fig_dict=None, fig_dict_wadatijack=None, ncores=0, metadata=None, origin=None):
    """
    :param data: ObsPy stream object containing waveform data of all stations in the event
    :type data: ~obspy.core.stream.Stream
    :param param: PylotParameter object containing parameters used for picking
    :type param: pylot.core.io.inputs.PylotParameter
    :param iplot: logical variable for plotting: 0=none, 1=partial, 2=all
    :type iplot: int, Boolean or String
    :param fig_dict: dictionary containing Matplotlib figures used for plotting picking results during tuning
    :type fig_dict: dict
    :param fig_dict_wadatijack: dictionary containing Matplotlib figures used for plotting jackknife-, wadati- and
    mediantest results
    :type fig_dict_wadatijack: dict
    :param ncores: number of cores used for parallel processing. Default (0) uses all available cores
    :type ncores: int
    :param metadata: tuple containing metadata type string and Parser object read from inventory file
    :type metadata: tuple (str, ~obspy.io.xseed.parser.Parser)
    :param origin: list containing origin objects representing origins for all events
    :type origin: list(~obspy.core.event.origin)
    :return: dictionary containing picked stations and pick information
    :rtype: dictionary
    """
    stations = []
    all_onsets = {}
    input_tuples = []
    try:
        iplot = int(iplot)
    except ValueError:
        if iplot is True or iplot == 'True':
            iplot = 2
        else:
            iplot = 0

    # get some parameters for quality control from
    # parameter input file (usually pylot.in).
    wdttolerance = param.get('wdttolerance')
    mdttolerance = param.get('mdttolerance')
    jackfactor = param.get('jackfactor')
    apverbose = param.get('apverbose')
    for n in range(len(data)):
        station = data[n].stats.station
        if station not in stations:
            stations.append(station)
        else:
            continue

    for station in stations:
        topick = data.select(station=station)
        input_tuples.append((topick, param, apverbose, iplot, fig_dict, metadata, origin))

    if iplot > 0:
        print('iPlot Flag active: NO MULTIPROCESSING possible.')
        ncores = 1

    # rename ncores for string representation in case ncores == 0 (use all cores)
    ncores_str = ncores if ncores != 0 else 'all available'
    print('Autopickstation: Distribute autopicking for {} '
          'stations on {} cores.'.format(len(input_tuples), ncores_str))

    if ncores == 1:
        results = serial_picking(input_tuples)
    else:
        results = parallel_picking(input_tuples, ncores)

    for result, station in results:
        if type(result) == dict:
            all_onsets[station] = result
        else:
            if result is None:
                result = 'Picker exited unexpectedly.'
            print('Could not pick a station: {}\nReason: {}'.format(station, result))

    # no Wadati/JK for single station (also valid for tuning mode)
    if len(stations) == 1:
        return all_onsets

    # quality control
    # median check and jackknife on P-onset times
    jk_checked_onsets = checkPonsets(all_onsets, mdttolerance, jackfactor, iplot, fig_dict_wadatijack)
    # check S-P times (Wadati)
    wadationsets = wadaticheck(jk_checked_onsets, wdttolerance, iplot, fig_dict_wadatijack)
    return wadationsets


def serial_picking(input_tuples):
    result = []
    for input_tuple in input_tuples:
        result.append(call_autopickstation(input_tuple))
    return result


def parallel_picking(input_tuples, ncores):
    pool = gen_Pool(ncores)
    result = pool.imap_unordered(call_autopickstation, input_tuples)
    pool.close()
    return result


def call_autopickstation(input_tuple):
    """
    helper function used for multiprocessing
    :param input_tuple: contains all parameters used for autopicking
    :type input_tuple: tuple
    :return: dictionary containing P pick, S pick and station name
    :rtype: dict
    """
    wfstream, pickparam, verbose, iplot, fig_dict, metadata, origin = input_tuple
    if fig_dict:
        print('Running in interactive mode')
    # multiprocessing not possible with interactive plotting
    try:
        return autopickstation(wfstream, pickparam, verbose, fig_dict=fig_dict, iplot=iplot, metadata=metadata,
                               origin=origin)
    except Exception as e:
        tbe = traceback.format_exc()
        return tbe, wfstream[0].stats.station


def get_source_coords(parser, station_id):
    """
    retrieves station coordinates from metadata
    :param parser: Parser object containing metadata read from inventory file
    :type parser: ~obspy.io.xseed.parser.Parser
    :param station_id: station id of which the coordinates should be retrieved, for example 'BW.RJOB..EHZ'. Only
    network and station name is required, channel id (last part) is ignored.
    :type station_id: str
    :return: dictionary containing 'latitude', 'longitude', 'elevation' and 'local_depth' of station
    :rtype: dict
    """
    station_coords = None
    try:
        station_coords = parser.get_coordinates(station_id)
    except Exception as e:
        print('Could not get source coordinates for station {}: {}'.format(station_id, e))
    return station_coords


class PickingResults(dict):
    """
    Used to store picking results.
    PickingResults is a dict like class that adds attribute (dot) access to the dictionaries values.
    """

    def __init__(self):
        """
        Inits default values of picking results. Called to generate a clean
        PickingResults instance with sensible defaults.
        :return:
        :rtype:
        """
        # Magnitude parameters
        self.Mo = None
        self.Mw = None

        # TODO What are those?
        self.w0 = None
        self.fc = None
        self.Ao = None  # Wood-Anderson peak-to-peak amplitude

        # Station information
        self.network = None
        self.channel = None

        # pick information
        self.picker = 'auto'  # type of pick
        self.marked = []

        # pick results
        self.epp = None  # earliest possible pick
        self.mpp = None  # most likely onset
        self.lpp = None  # latest possible pick
        self.fm = 'N'  # first motion polarity, can be set to 'U' (Up) or 'D' (Down)
        self.snr = None  # signal-to-noise ratio of onset
        self.snrdb = None  # signal-to-noise ratio of onset [dB]
        self.spe = None  # symmetrized picking error
        self.weight = 4  # weight of onset

    # to correctly provide dot access to dictionary attributes, all attribute access of the class is forwarded to the
    # dictionary
    def __getattr__(self, item):
        """Override getattr to return an AttributeError instead of a KeyError when the instance doesn't have the
        attribute.
        """
        try:
            attr = dict.__getitem__(self, item)
        except KeyError:
            raise AttributeError('{classname} has no attribute {attrname}'.format(classname=self.__class__.__name__,
                                                                                  attrname=item))
        return attr

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class PickingContainer:
    """
    Keeps intermediary results and values during picking
    """

    def __init__(self):
        # flags for plotting
        self.p_aic_plot_flag = 0
        self.aicSflag = 0
        self.Pflag = 0
        self.Sflag = 0


class AutopickStation(object):

    def __init__(self, wfstream, pickparam, verbose, iplot=0, fig_dict=None, metadata=None, origin=None):
        """
        :param wfstream: stream object containing waveform of all traces
        :type wfstream: ~obspy.core.stream.Stream
        :param pickparam: container of picking parameters from input file, usually pylot.in
        :type pickparam:  pylot.core.io.inputs.PylotParameter
        :param verbose: used to control output to log during picking. True = more information printed
        :type verbose: bool
        :param iplot: logical variable for plotting: 0=none, 1=partial, 2=all
        :type iplot: int, (Boolean or String)
        :param fig_dict: dictionary containing Matplotlib figures used for plotting picking results during tuning
        :type fig_dict: dict
        :param metadata: tuple containing metadata type string and Parser object read from inventory file
        :type metadata: tuple (str, ~obspy.io.xseed.parser.Parser)
        :param origin: list containing origin objects representing origins for all events
        :type origin: list(~obspy.core.event.origin)
        :return: dictionary-like object containing P pick, S pick and station name
        :rtype:
        """
        # save given parameters
        self.wfstream = wfstream
        self.pickparams = copy.deepcopy(pickparam)
        self.verbose = verbose
        self.iplot = correct_iplot(iplot)
        self.fig_dict = get_None(fig_dict)
        self.metadata = metadata
        self.origin = origin

        # initialize picking results
        self.p_results = PickingResults()
        self.s_results = PickingResults()

        # intialize containers that keep intermediary values between picking P- and S phase and plotting
        self.p_data = PickingContainer()
        self.s_data = PickingContainer()

        # extract additional information
        # TODO get channelorder from the pylot preferences
        self.channelorder = {'Z': 3, 'N': 1, 'E': 2}
        self.station_name = wfstream[0].stats.station
        self.network_name = wfstream[0].stats.network
        self.station_id = '{}.{}'.format(self.network_name, self.station_name)

        # save streams and traces
        self.zstream, self.nstream, self.estream = self.get_components_from_waveformstream()
        self.ztrace, self.ntrace, self.etrace = self.get_traces_from_streams()

        # default values used in old autopickstation function
        # #TODO way for user to set those
        self.detrend_type = 'demean'
        self.filter_type = 'bandpass'
        self.zerophase = False
        self.taper_max_percentage = 0.05
        self.taper_type = 'hann'

        # Used during picking to plot results
        self.current_figure = None
        self.current_linecolor = None

    def horizontal_traces_exist(self):
        """
        Return true when at least one horizontal trace exists
        :rtype: bool
        """
        if len(self.nstream) == len(self.estream) == 0:
            return False
        return True

    def vprint(self, s):
        """Only print statement if verbose picking is set to true."""
        if self.verbose:
            print(s)

    def get_components_from_waveformstream(self):
        """
        Splits waveformstream into multiple components zdat, ndat, edat. For traditional orientation (ZNE) these contain
        the vertical, north-south or east-west component. Otherwise they contain components numbered 123 with
        orientation diverging from the traditional orientation.
        :param waveformstream: Stream containing all three components for one station either by ZNE or 123 channel code
        (mixture of both options is handled as well)
        :type waveformstream: obspy.core.stream.Stream
        :return: Tuple containing (z waveform, n waveform, e waveform) selected by the given channels. If no waveform
        could be found for a given channel, an empty obspy Stream will be returned for that channel
        :rtype: (obspy.core.stream.Stream, obspy.core.stream.Stream, obspy.core.stream.Stream)
        """
        waveform_data = {}
        for key in self.channelorder:
            waveform_data[key] = self.wfstream.select(component=key)  # try ZNE first
            if len(waveform_data[key]) == 0:
                waveform_data[key] = self.wfstream.select(
                    component=str(self.channelorder[key]))  # use 123 as second option
        return waveform_data['Z'], waveform_data['N'], waveform_data['E']

    def get_traces_from_streams(self):
        """
        Extract Trace from Stream. If a component has no data, an empty trace will be returned
        :return: Tuple of obspy.Trace instances in order ZNE
        :rtype: (obspy.Trace)
        """
        if len(self.zstream) == 0:
            msg = 'No Z-component found for station {}. STOP'.format(self.wfstream[0].stats.station)
            raise MissingTraceException(msg)
        if not self.horizontal_traces_exist():
            # Both horizontal traces missing, only P pick can be determined
            msg = 'No horizontal traces found for station {}'.format(self.wfstream[0].stats.station)
            self.vprint(msg)
            return self.zstream[0], Trace(), Trace()

        ztrace = self.zstream[0]
        try:
            ntrace = self.nstream[0]
        except IndexError:
            # if N trace is missing, copy E trace
            self.nstream = self.estream
            ntrace = self.nstream[0]
        try:
            # if E trace is missing, copy N trace
            etrace = self.estream[0]
        except IndexError:
            self.estream = self.nstream
            etrace = self.estream[0]
        return ztrace, ntrace, etrace

    def prepare_wfstream(self, wfstream, freqmin=None, freqmax=None):
        """
        Prepare a waveformstream for picking by applying detrending, filtering and tapering. Creates a copy of the
        waveform the leave the original unchanged.
        :param wfstream: waveform stream
        :type wfstream: obspy.core.stream.Stream
        :param freqmin: Lower frequency of bandpass or highpass
        :type freqmin: float
        :param freqmax: Upper frequency of bandpass or lowpass
        :type freqmax: float
        :return: Tuple containing the changed waveform stream and the changed first trace of the stream
        :rtype: (obspy.core.trace.Trace, obspy.core.stream.Stream)
        """
        wfstream_copy = wfstream.copy()
        trace_copy = wfstream[0].copy()
        trace_copy.detrend(type=self.detrend_type)
        trace_copy.filter(self.filter_type, freqmin=freqmin, freqmax=freqmax, zerophase=self.zerophase)
        trace_copy.taper(max_percentage=self.taper_max_percentage, type=self.taper_type)
        wfstream_copy[0].data = trace_copy.data
        return trace_copy, wfstream_copy

    def modify_starttimes_taupy(self):
        """
        Calculate theoretical arrival times for a source at self.origin and a station at self.metadata. Modify
        self.pstart and self.pstop so they are based on a time window around these theoretical arrival times.
        :rtype: None
        """

        def get_seed_id():
            """
            Returns seed id of ztrace
            :return: Seed id with format Network.Station.Location.Channel
            :rtype: str
            """
            stats = self.ztrace.stats
            id = "{network}.{station}.{location}.{channel}"
            id = id.format(network=stats.network, station=stats.station, location=stats.location, channel=stats.channel)
            return id

        def create_arrivals(metadata, origin, taup_model):
            """
            Create List of arrival times for all phases for a given origin and station
            :param metadata: tuple containing metadata type string and Parser object read from inventory file
            :type metadata: tuple (str, ~obspy.io.xseed.parser.Parser)
            :param origin: list containing origin objects representing origins for all events
            :type origin: list(~obspy.core.event.origin)
            :param taup_model: Model name to use. See obspy.taup.tau.TauPyModel for options
            :type taup_model: str
            :return: List of Arrival objects
            :rtype: obspy.taup.tau.Arrivals
            :raises:
                AttributeError when no metadata or source origins is given
            """
            id = get_seed_id()
            station_coords = metadata.get_coordinates(id, self.ztrace.stats.starttime)
            if station_coords is None:
                exit_taupy()
                raise AttributeError('Warning: Could not find station in metadata')
            # TODO: raise when metadata.get_coordinates returns None
            source_origin = origin[0]
            model = TauPyModel(taup_model)
            taup_phases = self.pickparams['taup_phases']
            # for backward compatibility of older input files
            taup_phases = 'ttall' if not taup_phases or taup_phases == 'None' else taup_phases
            phase_list = [item.strip() for item in taup_phases.split(',')]
            arrivals = model.get_travel_times_geo(source_depth_in_km=source_origin.depth,
                                                  source_latitude_in_deg=source_origin.latitude,
                                                  source_longitude_in_deg=source_origin.longitude,
                                                  receiver_latitude_in_deg=station_coords['latitude'],
                                                  receiver_longitude_in_deg=station_coords['longitude'],
                                                  phase_list=phase_list)
            return arrivals

        def first_PS_onsets(arrivals):
            """
            Get first P and S onsets from arrivals list
            :param arrivals: List of Arrival objects
            :type arrivals: obspy.taup.tau.Arrivals
            :return:
            :rtype:
            """
            phases = {'P': [], 'S': []}
            # sort phases in P and S phases
            for arr in arrivals:
                phases[identifyPhaseID(arr.phase.name)].append(arr)
            # get first P and S onsets from arrivals list
            estFirstP = 0
            estFirstS = 0
            if len(phases['P']) > 0:
                arrP, estFirstP = min([(arr, arr.time) for arr in phases['P']], key=lambda t: t[1])
            if len(phases['S']) > 0:
                arrS, estFirstS = min([(arr, arr.time) for arr in phases['S']], key=lambda t: t[1])
            print('autopick: estimated first arrivals for P: {} s, S:{} s after event'
                  ' origin time using TauPy'.format(estFirstP, estFirstS))
            return estFirstP, estFirstS

        def exit_taupy():
            """If taupy failed to calculate theoretical starttimes, picking continues.
            For this a clean exit is required, since the P starttime is no longer relative to the theoretic onset but
            to the vertical trace starttime, eg. it can't be < 0."""
            # TODO here the pickparams is modified, instead of a copy
            if self.pickparams["pstart"] < 0:
                self.pickparams["pstart"] = 0
            if self.pickparams["sstart"] < 0:
                self.pickparams["sstart"] = 0

        if get_bool(self.pickparams["use_taup"]) is False:
            # correct user mistake where a relative cuttime is selected (pstart < 0) but use of taupy is disabled/ has
            # not the required parameters
            exit_taupy()
            return

        print('autopickstation: use_taup flag active.')
        # catch missing metadata or origin information. Onset calculation is stopped, given cuttimes are then used.
        if not self.metadata:
            raise AttributeError('Warning: Could not use TauPy to estimate onsets as there are no metadata given.')
        if not self.origin:
            raise AttributeError('No source origins given!')

        arrivals = create_arrivals(self.metadata, self.origin, self.pickparams["taup_model"])
        estFirstP, estFirstS = first_PS_onsets(arrivals)
        # modifiy pstart and pstop relative to estimated first P arrival (relative to station time axis)
        self.pickparams["pstart"] += (self.origin[0].time + estFirstP) - self.ztrace.stats.starttime
        self.pickparams["pstop"] += (self.origin[0].time + estFirstP) - self.ztrace.stats.starttime
        print('autopick: CF calculation times respectively:'
              ' pstart: {} s, pstop: {} s'.format(self.pickparams["pstart"], self.pickparams["pstop"]))
        # make sure pstart and pstop are inside the starttime/endtime of vertical trace
        self.pickparams["pstart"] = max(self.pickparams["pstart"], 0)
        self.pickparams["pstop"] = min(self.pickparams["pstop"], len(self.ztrace) * self.ztrace.stats.delta)

        if self.horizontal_traces_exist():
            # for the two horizontal components take earliest and latest time to make sure that the s onset is not clipped
            # if start and endtime of horizontal traces differ, the s windowsize will automatically increase
            trace_s_start = min([self.etrace.stats.starttime, self.ntrace.stats.starttime])
            # modifiy sstart and sstop relative to estimated first S arrival (relative to station time axis)
            self.pickparams["sstart"] += (self.origin[0].time + estFirstS) - trace_s_start
            self.pickparams["sstop"] += (self.origin[0].time + estFirstS) - trace_s_start
            print('autopick: CF calculation times respectively:'
                  ' sstart: {} s, sstop: {} s'.format(self.pickparams["sstart"], self.pickparams["sstop"]))
            # make sure pstart and pstop are inside the starttime/endtime of horizontal traces
            self.pickparams["sstart"] = max(self.pickparams["sstart"], 0)
            self.pickparams["sstop"] = min(self.pickparams["sstop"], len(self.ntrace) * self.ntrace.stats.delta,
                                           len(self.etrace) * self.etrace.stats.delta)

    def autopickstation(self):
        """
        Main function of autopickstation, which calculates P and S picks and returns them in a dictionary.
        :return: dict with keys 'P', 'S', and 'station'.
        P's value is a PickingResults instance containing P results.
        S's value is a PickingResults instance containing S results.
        station's value is the station name on which the picks were calculated.
        :rtype: dict
        """

        if get_bool(self.pickparams['use_taup']) is True and self.origin is not None:
            try:
                # modify pstart, pstop, sstart, sstop to be around theoretical onset if taupy should be used,
                # else do nothing
                self.modify_starttimes_taupy()
            except AttributeError as ae:
                print(ae)
            except MissingTraceException as mte:
                print(mte)

        try:
            self.pick_p_phase()
        except MissingTraceException as mte:
            print(mte)
        except PickingFailedException as pfe:
            print(pfe)

        if self.horizontal_traces_exist():
            if (self.p_results.weight is not None and self.p_results.weight < 4) or \
                    get_bool(self.pickparams.get('use_taup')):
                try:
                    self.pick_s_phase()
                except MissingTraceException as mte:
                    print(mte)
                except PickingFailedException as pfe:
                    print(pfe)

        self.plot_pick_results()
        self.finish_picking()
        return [{'P': self.p_results, 'S': self.s_results}, self.ztrace.stats.station]

    def finish_picking(self):

        # calculate "real" onset times, save them in PickingResults
        if self.p_results.lpp is not None and self.p_results.lpp == self.p_results.mpp:
            self.p_results.lpp += self.ztrace.stats.delta
        if self.p_results.epp is not None and self.p_results.epp == self.p_results.mpp:
            self.p_results.epp -= self.ztrace.stats.delta
        if self.p_results.mpp is not None and self.p_results.epp is not None and self.p_results.lpp is not None:
            self.p_results.lpp = self.ztrace.stats.starttime + self.p_results.lpp
            self.p_results.epp = self.ztrace.stats.starttime + self.p_results.epp
            self.p_results.mpp = self.ztrace.stats.starttime + self.p_results.mpp
        else:
            # dummy values (start of seismic trace) in order to derive
            # theoretical onset times for iterative picking
            self.p_results.lpp = self.ztrace.stats.starttime + self.pickparams["timeerrorsP"][3]
            self.p_results.epp = self.ztrace.stats.starttime - self.pickparams["timeerrorsP"][3]
            self.p_results.mpp = self.ztrace.stats.starttime
        self.p_results.channel = self.ztrace.stats.channel
        self.p_results.network = self.ztrace.stats.network
        #
        #   S results
        #
        if not self.horizontal_traces_exist():
            # no horizontal components means there is no S pick to be finished
            return
        if self.etrace:
            hdat = self.etrace
        elif self.ntrace:
            hdat = self.ntrace

        if self.s_results.lpp is not None and self.s_results.lpp == self.s_results.mpp:
            self.s_results.lpp += hdat.stats.delta
        if self.s_results.epp is not None and self.s_results.epp == self.s_results.mpp:
            self.s_results.epp -= hdat.stats.delta
        if self.s_results.mpp is not None and self.s_results.epp is not None and self.s_results.lpp is not None:
            self.s_results.lpp = hdat.stats.starttime + self.s_results.lpp
            self.s_results.epp = hdat.stats.starttime + self.s_results.epp
            self.s_results.mpp = hdat.stats.starttime + self.s_results.mpp
        else:
            # dummy values (start of seismic trace) in order to derive
            # theoretical onset times for iteratve picking
            self.s_results.lpp = hdat.stats.starttime + self.pickparams["timeerrorsS"][3]
            self.s_results.epp = hdat.stats.starttime - self.pickparams["timeerrorsS"][3]
            self.s_results.mpp = hdat.stats.starttime

        self.s_results.channel = self.etrace.stats.channel
        self.s_results.network = self.etrace.stats.network
        self.s_results.fm = None  # override default value 'N'

    def plot_pick_results(self):
        if self.iplot > 0:
            # plot vertical trace
            if self.fig_dict is None:
                fig = plt.figure()
                plt_flag = 1
                linecolor = 'k'
            else:
                fig = self.fig_dict['mainFig']
                linecolor = self.fig_dict['plot_style']['linecolor']['rgba_mpl']
                plt_flag = 0
            fig._tight = True
            ax1 = fig.add_subplot(311)
            tdata = np.linspace(start=0, stop=self.ztrace.stats.endtime - self.ztrace.stats.starttime,
                                num=self.ztrace.stats.npts)
            # plot tapered trace filtered with bpz2 filter settings
            ax1.plot(tdata, self.tr_filt_z_bpz2.data / max(self.tr_filt_z_bpz2.data), color=linecolor, linewidth=0.7,
                     label='Data')
            if self.p_results.weight < 4:
                # plot CF of initial onset (HOScf or ARZcf)
                ax1.plot(self.cf1.getTimeArray(), self.cf1.getCF() / max(self.cf1.getCF()), 'b', label='CF1')
            if self.p_data.p_aic_plot_flag == 1:
                aicpick = self.p_data.aicpick
                refPpick = self.p_data.refPpick
                # plot CF of precise pick (HOScf or ARZcf)
                ax1.plot(self.cf2.getTimeArray(), self.cf2.getCF() / max(self.cf2.getCF()), 'm', label='CF2')
                # plot inital P pick
                ax1.plot([aicpick.getpick(), aicpick.getpick()], [-1, 1], 'r', label='Initial P Onset')
                ax1.plot([aicpick.getpick() - 0.5, aicpick.getpick() + 0.5], [1, 1], 'r')
                ax1.plot([aicpick.getpick() - 0.5, aicpick.getpick() + 0.5], [-1, -1], 'r')
                # plot precise P pick
                ax1.plot([refPpick.getpick(), refPpick.getpick()], [-1.3, 1.3], 'r', linewidth=2, label='Final P Pick')
                ax1.plot([refPpick.getpick() - 0.5, refPpick.getpick() + 0.5], [1.3, 1.3], 'r', linewidth=2)
                ax1.plot([refPpick.getpick() - 0.5, refPpick.getpick() + 0.5], [-1.3, -1.3], 'r', linewidth=2)
                # plot latest possible P pick
                ax1.plot([self.p_results.lpp, self.p_results.lpp], [-1.1, 1.1], 'r--', label='lpp')
                # plot earliest possible P pick
                ax1.plot([self.p_results.epp, self.p_results.epp], [-1.1, 1.1], 'r--', label='epp')
                # add title to plot
                title = '{station}, {channel}, P weight={pweight:d}, SNR={snr:7.2}, SNR[dB]={snrdb:7.2}, Polarity: {polarity}'
                ax1.set_title(title.format(station=self.ztrace.stats.station,
                                           channel=self.ztrace.stats.channel,
                                           pweight=self.p_results.weight,
                                           snr=self.p_results.snr,
                                           snrdb=self.p_results.snrdb,
                                           polarity=self.p_results.fm))
            else:
                title = '{channel}, P weight={pweight}, SNR=None, SNR[dB]=None'
                ax1.set_title(title.format(channel=self.ztrace.stats.channel, pweight=self.p_results.weight))

            ax1.legend(loc=1)
            ax1.set_yticks([])
            ax1.set_ylim([-1.5, 1.5])
            ax1.set_ylabel('Normalized Counts')

            if self.horizontal_traces_exist():# and self.s_data.Sflag == 1:
                # plot E trace
                ax2 = fig.add_subplot(3, 1, 2, sharex=ax1)
                th1data = np.linspace(0, self.etrace.stats.endtime - self.etrace.stats.starttime,
                                      self.etrace.stats.npts)
                # plot filtered and tapered waveform
                ax2.plot(th1data, self.etrace.data / max(self.etrace.data), color=linecolor, linewidth=0.7,
                         label='Data')
                if self.p_results.weight < 4:
                    # plot initial CF (ARHcf or AR3Ccf)
                    ax2.plot(self.arhcf1.getTimeArray(), self.arhcf1.getCF() / max(self.arhcf1.getCF()), 'b',
                             label='CF1')
                    if self.s_data.aicSflag == 1 and self.s_results.weight <= 4:
                        aicarhpick = self.aicarhpick
                        refSpick = self.refSpick
                        # plot second cf, used for determing precise onset (ARHcf or AR3Ccf)
                        ax2.plot(self.arhcf2.getTimeArray(), self.arhcf2.getCF() / max(self.arhcf2.getCF()), 'm',
                                 label='CF2')
                        # plot preliminary onset time, calculated from CF1
                        ax2.plot([aicarhpick.getpick(), aicarhpick.getpick()], [-1, 1], 'g', label='Initial S Onset')
                        ax2.plot([aicarhpick.getpick() - 0.5, aicarhpick.getpick() + 0.5], [1, 1], 'g')
                        ax2.plot([aicarhpick.getpick() - 0.5, aicarhpick.getpick() + 0.5], [-1, -1], 'g')
                        # plot precise onset time, calculated from CF2
                        ax2.plot([refSpick.getpick(), refSpick.getpick()], [-1.3, 1.3], 'g', linewidth=2,
                                 label='Final S Pick')
                        ax2.plot([refSpick.getpick() - 0.5, refSpick.getpick() + 0.5], [1.3, 1.3], 'g', linewidth=2)
                        ax2.plot([refSpick.getpick() - 0.5, refSpick.getpick() + 0.5], [-1.3, -1.3], 'g', linewidth=2)
                        ax2.plot([self.s_results.lpp, self.s_results.lpp], [-1.1, 1.1], 'g--', label='lpp')
                        ax2.plot([self.s_results.epp, self.s_results.epp], [-1.1, 1.1], 'g--', label='epp')
                        title = '{channel}, S weight={sweight}, SNR={snr:7.2}, SNR[dB]={snrdb:7.2}'
                        ax2.set_title(title.format(channel=str(self.etrace.stats.channel),
                                                   sweight=str(self.s_results.weight),
                                                   snr=str(self.s_results.snr),
                                                   snrdb=str(self.s_results.snrdb)))
                    else:
                        title = '{channel}, S weight={sweight}, SNR=None, SNR[dB]=None'
                        ax2.set_title(title.format(channel=str(self.etrace.stats.channel),
                                                   sweight=str(self.s_results.weight)))
                ax2.legend(loc=1)
                ax2.set_yticks([])
                ax2.set_ylim([-1.5, 1.5])
                ax2.set_ylabel('Normalized Counts')

                # plot N trace
                ax3 = fig.add_subplot(3, 1, 3, sharex=ax1)
                th2data = np.linspace(0, self.ntrace.stats.endtime - self.ntrace.stats.starttime,
                                      self.ntrace.stats.npts)
                # plot trace
                ax3.plot(th2data, self.ntrace.data / max(self.ntrace.data), color=linecolor, linewidth=0.7,
                         label='Data')
                if self.p_results.weight < 4:
                    p22, = ax3.plot(self.arhcf1.getTimeArray(), self.arhcf1.getCF() / max(self.arhcf1.getCF()), 'b',
                                    label='CF1')
                    if self.s_data.aicSflag == 1:
                        aicarhpick = self.aicarhpick
                        refSpick = self.refSpick
                        ax3.plot(self.arhcf2.getTimeArray(), self.arhcf2.getCF() / max(self.arhcf2.getCF()), 'm',
                                 label='CF2')
                        ax3.plot([aicarhpick.getpick(), aicarhpick.getpick()], [-1, 1], 'g', label='Initial S Onset')
                        ax3.plot([aicarhpick.getpick() - 0.5, aicarhpick.getpick() + 0.5], [1, 1], 'g')
                        ax3.plot([aicarhpick.getpick() - 0.5, aicarhpick.getpick() + 0.5], [-1, -1], 'g')
                        ax3.plot([refSpick.getpick(), refSpick.getpick()], [-1.3, 1.3], 'g', linewidth=2,
                                 label='Final S Pick')
                        ax3.plot([refSpick.getpick() - 0.5, refSpick.getpick() + 0.5], [1.3, 1.3], 'g', linewidth=2)
                        ax3.plot([refSpick.getpick() - 0.5, refSpick.getpick() + 0.5], [-1.3, -1.3], 'g', linewidth=2)
                        ax3.plot([self.s_results.lpp, self.s_results.lpp], [-1.1, 1.1], 'g--', label='lpp')
                        ax3.plot([self.s_results.epp, self.s_results.epp], [-1.1, 1.1], 'g--', label='epp')
                ax3.legend(loc=1)
                ax3.set_yticks([])
                ax3.set_ylim([-1.5, 1.5])
                ax3.set_xlabel('Time [s] after %s' % self.ntrace.stats.starttime)
                ax3.set_ylabel('Normalized Counts')
                ax3.set_title(self.ntrace.stats.channel)
                if plt_flag == 1:
                    fig.show()
                    try:
                        input()
                    except SyntaxError:
                        pass
                    plt.close(fig)

    def _pick_p_quality_control(self, aicpick, z_copy, tr_filt):
        """
        Quality control of first pick using minseglength and checkZ4S.
        :param aicpick: Instance of AICPicker to run quality control on
        :type aicpick: AICPicker
        :param z_copy: Stream if vertical trace, data replaced with values from from initial CF (HOScf or ARHcf)
        :type z_copy: obspy.core.stream.Stream
        :param tr_filt: Filtered and tapered trace of vertical component
        :type tr_filt: obspy.core.trace.trace
        :return: Flag if P onset passed quality control, 1 if passed, 0 if failed.
        :rtype: int
        """

        self.set_current_figure('slength')
        if aicpick.getpick() is None:
            msg = "Bad initial (AIC) P-pick, skipping this onset!\nAIC-SNR={0}, AIC-Slope={1}counts/s\n " \
                  "(min. AIC-SNR={2}, min. AIC-Slope={3}counts/s)"
            msg = msg.format(aicpick.getSNR(), aicpick.getSlope(), self.pickparams["minAICPSNR"],
                             self.pickparams["minAICPslope"])
            self.vprint(msg)
            return 0
        # Quality check initial pick with minimum signal length
        z_copy[0].data = tr_filt.data  # save filtered, tapered trace in z_copy stream object
        zne = z_copy
        minsiglength = self.pickparams["minsiglength"]
        if len(self.nstream) == 0 or len(self.estream) == 0:
            msg = 'One or more horizontal component(s) missing!\n' \
                  'Signal length only checked on vertical component!\n' \
                  'Decreasing minsiglengh from {0} to {1}' \
                .format(minsiglength, minsiglength / 2)
            self.vprint(msg)
            minsiglength = minsiglength / 2
        else:
            # filter, taper other traces as well since signal length is compared on all traces
            trH1_filt, _ = self.prepare_wfstream(self.estream, freqmin=self.pickparams["bph1"][0],
                                                 freqmax=self.pickparams["bph1"][1])
            trH2_filt, _ = self.prepare_wfstream(self.nstream, freqmin=self.pickparams["bph1"][0],
                                                 freqmax=self.pickparams["bph1"][1])
            zne += trH1_filt
            zne += trH2_filt
            minsiglength = minsiglength
        Pflag = checksignallength(zne, aicpick.getpick(), minsiglength, self.pickparams,
                                  self.iplot, self.current_figure, self.current_linecolor)
        if Pflag == 0:
            self.p_results.marked = 'shortsignallength'
            self.p_results.weight = 9
            return 0

        if self.nstream == self.estream:
            # todo: old implementation skipped this test if one component was missing, why not use one component?
            msg = 'One or more horizontal components missing!\n Skipping control function checkZ4S.'
            self.vprint(msg)
            return 1

        if self.iplot > 1: self.set_current_figure('checkZ4s')
        Pflag = checkZ4S(zne, aicpick.getpick(), self.pickparams, self.iplot,
                         self.current_figure, self.current_linecolor)
        if Pflag == 0:
            self.p_results.marked = 'SinsteadP'
            self.p_results.weight = 9
            return 0
        return 1

    def pick_p_phase(self):
        """
        Pick p phase, store results in self.p_results
        :return: None
        :raises:
            MissingTraceException: If vertical trace is missing.
        """
        try:
            stream_string_repr = str(self.zstream)
        except Exception as e:
            stream_string_repr = None
            print('Could not get string representation of Stream: {}'.format(e))
        msg = '##################################################\nautopickstation:' \
              ' Working on P onset of station {station}\nFiltering vertical ' \
              'trace ...\n{data}'.format(station=self.station_name, data=stream_string_repr)
        self.vprint(msg)

        tr_filt, z_copy = self.prepare_wfstream(self.zstream, self.pickparams["bpz1"][0], self.pickparams["bpz1"][1])
        # save filtered trace in instance for later plotting
        self.tr_filt_z_bpz2 = tr_filt

        Lc = self.pickparams['pstop'] - self.pickparams['pstart']

        Lwf = self.ztrace.stats.endtime - self.ztrace.stats.starttime
        if Lwf < 0:
            raise PickingFailedException('autopickstation: empty trace! Return!')

        Ldiff = Lwf - abs(Lc)
        if Ldiff < self.pickparams['tlta'] or self.pickparams['pstop'] <= self.pickparams['pstart']:
            msg = 'autopickstation: Cutting times are too large for actual waveform!\nUsing entire waveform instead!'
            self.vprint(msg)
            self.pickparams['pstart'] = 0
            self.pickparams['pstop'] = len(self.ztrace.data) * self.ztrace.stats.delta

        cuttimes = self._calculate_cuttimes('P', 1)

        # calculate first CF
        if self.pickparams["algoP"] == 'HOS':
            self.cf1 = HOScf(z_copy, cuttimes, self.pickparams)
        elif self.pickparams["algoP"] == 'ARZ':
            self.cf1 = ARZcf(z_copy, cuttimes, self.pickparams["tdet1z"], self.pickparams["tpred1z"], self.pickparams)
        else:
            self.cf1 = None
        assert isinstance(self.cf1, CharacteristicFunction), 'cf1 is not set correctly: maybe the algorithm name ({})' \
                                                             ' is corrupted'.format(self.pickparams["algoP"])
        # calculate AIC cf from first cf (either HOS or ARZ)
        z_copy[0].data = self.cf1.getCF()
        aiccf = AICcf(z_copy, cuttimes)
        # get preliminary onset time from AIC-CF
        self.set_current_figure('aicFig')
        aicpick = AICPicker(aiccf, self.pickparams["tsnrz"], self.pickparams["pickwinP"], self.iplot,
                            Tsmooth=self.pickparams["aictsmooth"], fig=self.current_figure,
                            linecolor=self.current_linecolor)
        # save aicpick for plotting later
        self.p_data.aicpick = aicpick
        # add pstart and pstop to aic plot
        if self.current_figure:
            # TODO remove plotting from picking, make own plot function
            for ax in self.current_figure.axes:
                ax.vlines(self.pickparams["pstart"], ax.get_ylim()[0], ax.get_ylim()[1], color='c', linestyles='dashed',
                          label='P start')
                ax.vlines(self.pickparams["pstop"], ax.get_ylim()[0], ax.get_ylim()[1], color='c', linestyles='dashed',
                          label='P stop')
                ax.legend(loc=1)

        Pflag = self._pick_p_quality_control(aicpick, z_copy, tr_filt)
        # go on with processing if AIC onset passes quality control
        slope = aicpick.getSlope()
        if not slope: slope = 0
        # todo why did picking fail was saved in the pick dictionary, should this be reimplemented?
        if Pflag != 1:
            raise PickingFailedException('AIC P onset quality control failed')
        if slope <= self.pickparams["minAICPslope"]:
            error_msg = 'AIC P onset slope to small: got {}, min {}'.format(slope, self.pickparams["minAICPslope"])
            raise PickingFailedException(error_msg)
        if aicpick.getSNR() < self.pickparams["minAICPSNR"]:
            error_msg = 'AIC P onset SNR to small: got {}, min {}'.format(aicpick.getSNR(),
                                                                          self.pickparams["minAICPSNR"])
            raise PickingFailedException(error_msg)

        self.p_data.p_aic_plot_flag = 1
        msg = 'AIC P-pick passes quality control: Slope: {0} counts/s, SNR: {1}\nGo on with refined picking ...\n' \
              'autopickstation: re-filtering vertical trace...'.format(aicpick.getSlope(), aicpick.getSNR())
        self.vprint(msg)
        # refilter waveform with larger bandpass
        tr_filt, z_copy = self.prepare_wfstream(self.zstream, freqmin=self.pickparams["bpz2"][0],
                                                freqmax=self.pickparams["bpz2"][1])
        # save filtered trace in instance for later plotting
        self.tr_filt_z_bpz2 = tr_filt
        # determine new times around initial onset
        cuttimes2 = self._calculate_cuttimes('P', 2)
        if self.pickparams["algoP"] == 'HOS':
            self.cf2 = HOScf(z_copy, cuttimes2, self.pickparams)
        elif self.pickparams["algoP"] == 'ARZ':
            self.cf2 = ARZcf(z_copy, cuttimes2, self.pickparams["tdet2z"], self.pickparams["tpred2z"], self.pickparams)
        else:
            self.cf2 = None
        assert isinstance(self.cf2, CharacteristicFunction), 'cf2 is not set correctly: maybe the algorithm name () is ' \
                                                             'corrupted'.format(self.pickparams["algoP"])
        self.set_current_figure('refPpick')
        # get refined onset time from CF2
        refPpick = PragPicker(self.cf2, self.pickparams["tsnrz"], self.pickparams["pickwinP"], self.iplot,
                              self.pickparams["ausP"],
                              self.pickparams["tsmoothP"], aicpick.getpick(), self.current_figure,
                              self.current_linecolor)
        # save PragPicker result for plotting
        self.p_data.refPpick = refPpick
        self.p_results.mpp = refPpick.getpick()
        if self.p_results.mpp is None:
            msg = 'Bad initial (AIC) P-pick, skipping this onset!\n AIC-SNR={}, AIC-Slope={}counts/s\n' \
                  '(min. AIC-SNR={}, min. AIC-Slope={}counts/s)'
            msg.format(aicpick.getSNR(), aicpick.getSlope(), self.pickparams["minAICPSNR"],
                       self.pickparams["minAICPslope"])
            self.vprint(msg)
            self.s_data.Sflag = 0
            raise PickingFailedException(msg)
        # quality assessment, get earliest/latest pick and symmetrized uncertainty
        # todo quality assessment in own function
        self.set_current_figure('el_Ppick')
        elpicker_results = earllatepicker(z_copy, self.pickparams["nfacP"], self.pickparams["tsnrz"],
                                          self.p_results.mpp,
                                          self.iplot, fig=self.current_figure, linecolor=self.current_linecolor)
        self.p_results.epp, self.p_results.lpp, self.p_results.spe = elpicker_results
        snr_results = getSNR(z_copy, self.pickparams["tsnrz"], self.p_results.mpp)
        self.p_results.snr, self.p_results.snrdb, _ = snr_results

        # weight P-onset using symmetric error
        self.p_results.weight = get_quality_class(self.p_results.spe, self.pickparams["timeerrorsP"])
        if self.p_results.weight <= self.pickparams["minfmweight"] and self.p_results.snr >= self.pickparams[
            "minFMSNR"]:
            # if SNR is high enough, try to determine first motion of onset
            self.set_current_figure('fm_picker')
            self.p_results.fm = fmpicker(self.zstream.copy(), z_copy, self.pickparams["fmpickwin"], self.p_results.mpp,
                                         self.iplot, self.current_figure, self.current_linecolor)
        msg = "autopickstation: P-weight: {}, SNR: {}, SNR[dB]: {}, Polarity: {}"
        msg = msg.format(self.p_results.weight, self.p_results.snr, self.p_results.snrdb, self.p_results.fm)
        print(msg)
        msg = 'autopickstation: Refined P-Pick: {} s | P-Error: {} s'
        msg = msg.format(self.p_results.mpp, self.p_results.spe)
        print(msg)
        self.s_data.Sflag = 1

    def _calculate_cuttimes(self, type, iteration):
        """
        Calculate cuttimes for a trace
        :param type: 'P' or 'S', denoting the pick for which cuttime should be calculated
        :type type: str
        :param iteration: Calculate cut times for initial pick or for the smaller window of the precise pick around
        the initial pick
        :type iteration: int
        :return: tuple of (starttime, endtime) in seconds
        :rtype: (int, int)
        """
        # extract parameters
        Precalcwin = self.pickparams["Precalcwin"]
        Srecalcwin = self.pickparams["Srecalcwin"]

        if type.upper() == 'P':
            if iteration == 1:
                return [self.pickparams["pstart"], self.pickparams["pstop"]]
            if iteration == 2:
                starttime2 = round(max(self.p_data.aicpick.getpick() - Precalcwin, 0))
                endtime2 = round(
                    min(len(self.ztrace.data) * self.ztrace.stats.delta, self.p_data.aicpick.getpick() + Precalcwin))
                return [starttime2, endtime2]
        elif type.upper() == 'S':
            if iteration == 1:
                # Calculate start times for preliminary S onset
                start = round(max(self.p_results.mpp + self.pickparams["sstart"], 0))  # limit start time to >0 seconds
                stop = round(min([
                    self.p_results.mpp + self.pickparams["sstop"],
                    self.etrace.stats.endtime - self.etrace.stats.starttime,
                    self.ntrace.stats.endtime - self.ntrace.stats.starttime
                ]))
                cuttimesh = (start, stop)
                if cuttimesh[1] <= cuttimesh[0]:
                    raise PickingFailedException('Cut window for horizontal phases too small! Will not pick S onsets.')
                return cuttimesh
            if iteration == 2:
                # recalculate cf from refiltered trace in vicinity of initial onset
                start = round(self.aicarhpick.getpick() - Srecalcwin)
                stop = round(self.aicarhpick.getpick() + Srecalcwin)
                return (start, stop)
        else:
            raise ValueError('Wrong type given, can only be P or S')

    def _calculate_autoregressive_cf_s_pick(self, cuttimesh):
        algoS = self.pickparams["algoS"]
        filter_freq_min, filter_freq_max = self.pickparams["bph1"]
        # prepare traces for picking by filtering, taper
        if algoS == 'ARH':
            self.hdat = self.nstream.copy() + self.estream.copy()
            trH1_filt, _ = self.prepare_wfstream(self.estream, filter_freq_min, filter_freq_max)
            trH2_filt, _ = self.prepare_wfstream(self.nstream, filter_freq_min, filter_freq_max)
            h_copy = self.hdat.copy()
            h_copy[0].data = trH1_filt.data
            h_copy[1].data = trH2_filt.data
        if algoS == 'AR3':
            self.hdat = self.zstream.copy() + self.estream.copy() + self.nstream.copy()
            trH1_filt, _ = self.prepare_wfstream(self.zstream, filter_freq_min, filter_freq_max)
            trH2_filt, _ = self.prepare_wfstream(self.estream, filter_freq_min, filter_freq_max)
            trH3_filt, _ = self.prepare_wfstream(self.nstream, filter_freq_min, filter_freq_max)
            h_copy = self.hdat.copy()
            h_copy[0].data = trH1_filt.data
            h_copy[1].data = trH2_filt.data
            h_copy[2].data = trH3_filt.data

        self.trH1_filt = trH1_filt
        self.h_copy = h_copy

        # calculate initial CF based on autoregression
        if algoS == 'ARH':
            arhcf1 = ARHcf(h_copy, cuttimesh, self.pickparams["tdet1h"], self.pickparams["tpred1h"], self.pickparams)
        elif algoS == 'AR3':
            arhcf1 = AR3Ccf(h_copy, cuttimesh, self.pickparams["tdet1h"], self.pickparams["tpred1h"], self.pickparams)
        return arhcf1

    def _calculate_aic_cf_s_pick(self, cuttimesh):
        stream = self.estream.copy()
        stream[0].data = self.arhcf1.getCF()
        haiccf = AICcf(stream, cuttimesh)
        return haiccf

    def _pick_s_quality_control(self):
        """
        Check quality of pick. Function will raise a PickingFailedException if the S pick does not fullfill all quality
        criteria. Else nothing happens and picking can continue.
        """
        # go on with processing if AIC onset passes quality control
        slope = self.aicarhpick.getSlope()
        minSlope = self.pickparams["minAICSslope"]
        minSNR = self.pickparams["minAICSSNR"]

        if not slope:
            slope = 0

        if slope < minSlope:
            error_msg = error_msg = 'AIC S onset slope to small: got {}, min {}'.format(slope, minSlope)
            raise PickingFailedException(error_msg)
        if self.aicarhpick.getSNR() < minSNR:
            error_msg = 'AIC S onset SNR to small: got {}, min {}'.format(self.aicarhpick.getSNR(), minSNR)
            raise PickingFailedException(error_msg)
        if self.aicarhpick.getpick() is None:
            error_msg = 'Invalid AIC S pick!'
            raise PickingFailedException(error_msg)
        self.s_data.aicSflag = 1
        msg = 'AIC S-pick passes quality control: Slope: {0} counts/s, ' \
              'SNR: {1}\nGo on with refined picking ...\n' \
              'autopickstation: re-filtering horizontal traces ' \
              '...'.format(self.aicarhpick.getSlope(), self.aicarhpick.getSNR())
        self.vprint(msg)

    def _pick_s_calculate_ar_cf_2(self):
        algoS = self.pickparams["algoS"]
        filter_freq_min, filter_freq_max = self.pickparams["bph2"]
        cuttimesh2 = self._calculate_cuttimes('S', 2)
        # refilter waveform with larger bandpass
        trH1_filt, _ = self.prepare_wfstream(self.estream, filter_freq_min, filter_freq_max)
        trH2_filt, _ = self.prepare_wfstream(self.nstream, filter_freq_min, filter_freq_max)
        if algoS == 'ARH':
            h_copy = self.hdat.copy()
            h_copy[0].data = trH1_filt.data
            h_copy[1].data = trH2_filt.data
        elif algoS == 'AR3':
            trH3_filt, _ = self.prepare_wfstream(self.zstream, filter_freq_min, filter_freq_max)
            h_copy = self.hdat.copy()
            h_copy[0].data = trH3_filt.data
            h_copy[1].data = trH1_filt.data
            h_copy[2].data = trH2_filt.data

        # save filtered traces for plotting
        self.estream_bph2 = trH1_filt
        self.nstream_bph2 = trH2_filt

        # calculate second cf
        if algoS == 'ARH':
            arhcf2 = ARHcf(h_copy, cuttimesh2, self.pickparams["tdet2h"], self.pickparams["tpred2h"], self.pickparams)
        elif algoS == 'AR3':
            arhcf2 = AR3Ccf(h_copy, cuttimesh2, self.pickparams["tdet2h"], self.pickparams["tpred2h"], self.pickparams)
        # save cf for later plotting
        self.arhcf2 = arhcf2

        self.h_copy = h_copy
        return arhcf2

    def _pick_s_quality_assessment(self, h_copy):
        """
        quality assessment: get earliest/latest possible pick and symmetrized uncertainty
        """
        h_copy[0].data = self.estream_bph2.data
        if self.iplot:
            self.set_current_figure('el_S1pick')
        epickS1, lpickS1, Serror1 = earllatepicker(h_copy, self.pickparams["nfacS"], self.pickparams["tsnrh"],
                                                   self.s_results.mpp, self.iplot, fig=self.current_figure,
                                                   linecolor=self.current_linecolor)

        h_copy[0].data = self.nstream_bph2.data
        if self.iplot:
            self.set_current_figure('el_S2pick')
        else:
            # why is it set to empty here? DA
            linecolor = ''
        epickS2, lpickS2, Serror2 = earllatepicker(h_copy, self.pickparams["nfacS"], self.pickparams["tsnrh"],
                                                   self.s_results.mpp, self.iplot, fig=self.current_figure,
                                                   linecolor=self.current_linecolor)

        if epickS1 is not None and epickS2 is not None:
            if self.pickparams["algoS"] == 'ARH':
                # get earliest pick of both earliest possible picks
                epick = [epickS1, epickS2]
                lpick = [lpickS1, lpickS2]
                pickerr = [Serror1, Serror2]
                ipick = np.argmin(epick)
            if self.pickparams["algoS"] == 'AR3':
                epickS3, lpickS3, Serror3 = earllatepicker(h_copy, self.pickparams["nfacS"], self.pickparams["tsnrh"],
                                                           self.s_results.mpp, self.iplot)
                # get earliest of all three picks
                epick = [epickS1, epickS2, epickS3]
                lpick = [lpickS1, lpickS2, lpickS3]
                pickerr = [Serror1, Serror2, Serror3]

                if epickS3 is not None:
                    ipick = np.argmin(epick)
                else:
                    ipick = np.argmin([epickS1, epickS2])
            self.s_results.epp = epick[ipick]
            self.s_results.lpp = lpick[ipick]
            self.s_results.spe = pickerr[ipick]

            msg = 'autopickstation: Refined S-Pick: {} s | S-Error: {} s'.format(self.s_results.mpp,
                                                                                 self.s_results.spe)
            print(msg)

            # get SNR
            self.s_results.snr, self.s_results.snrdb, _ = getSNR(h_copy, self.pickparams["tsnrh"], self.s_results.mpp)

            self.s_results.weight = get_quality_class(self.s_results.spe, self.pickparams["timeerrorsS"])

            print('autopickstation: S-weight: {0}, SNR: {1}, '
                  'SNR[dB]: {2}\n'
                  '##################################################'
                  ''.format(self.s_results.weight, self.s_results.snr, self.s_results.snrdb))

    def pick_s_phase(self):
        if get_bool(self.pickparams.get('use_taup')) is True:
            cuttimesh = (self.pickparams.get('sstart'), self.pickparams.get('sstop'))
        else:
            # determine time window for calculating CF after P onset
            cuttimesh = self._calculate_cuttimes(type='S', iteration=1)

        # calculate autoregressive CF
        self.arhcf1 = self._calculate_autoregressive_cf_s_pick(cuttimesh)

        # calculate AIC cf
        haiccf = self._calculate_aic_cf_s_pick(cuttimesh)

        # get preliminary onset time from AIC cf
        self.set_current_figure('aicARHfig')
        aicarhpick = AICPicker(haiccf, self.pickparams["tsnrh"], self.pickparams["pickwinS"], self.iplot,
                               Tsmooth=self.pickparams["aictsmoothS"], fig=self.current_figure,
                               linecolor=self.current_linecolor)
        # save pick for later plotting
        self.aicarhpick = aicarhpick

        # check quality of pick
        self._pick_s_quality_control()

        arhcf2 = self._pick_s_calculate_ar_cf_2()

        # get refined onset time from CF2
        self.set_current_figure('refSpick')
        refSpick = PragPicker(arhcf2, self.pickparams["tsnrh"], self.pickparams["pickwinS"], self.iplot,
                              self.pickparams["ausS"],
                              self.pickparams["tsmoothS"], aicarhpick.getpick(), self.current_figure,
                              self.current_linecolor)
        # save refSpick for later plotitng
        self.refSpick = refSpick
        self.s_results.mpp = refSpick.getpick()

        if self.s_results.mpp is not None:
            self._pick_s_quality_assessment(self.h_copy)

    def set_current_figure(self, figkey):
        """
        Extracts a figure by name from dictionary and set it as the currently active figure.
        All functions that create plots during picking will use the currently active figure to plot them.
        :param figkey:
        :type figkey:
        :return:
        :rtype:
        """
        if self.fig_dict is None:
            return None, None
        self.current_figure = self.fig_dict.get(figkey, None)
        plot_style = self.fig_dict.get('plot_style', 'k')
        self.current_linecolor = plot_style['linecolor']['rgba_mpl']


def autopickstation(wfstream, pickparam, verbose=False, iplot=0, fig_dict=None, metadata=None, origin=None):
    """
    Main function to calculate picks for the station.
    :return:
    :rtype: dict
    """
    try:
        station = AutopickStation(wfstream, pickparam, verbose, iplot, fig_dict, metadata, origin)
        return station.autopickstation()
    except MissingTraceException as e:
        # Either vertical or both horizontal traces are missing
        print(e)
        try:
            station_name = wfstream[0].stats.station
        except IndexError:
            station_name = 'None'
        return None, station_name


def iteratepicker(wf, NLLocfile, picks, badpicks, pickparameter, fig_dict=None):
    """
    Repicking of bad onsets. Uses theoretical onset times from NLLoc-location file.
    :param wf: waveform, obspy stream object
    :type wf: ~obspy.core.stream.Stream
    :param NLLocfile: path/name of NLLoc-location file
    :type NLLocfile: str
    :param picks: dictionary of available onset times
    :type picks: dict
    :param badpicks: picks to be repicked
    :type badpicks:
    :param pickparameter: picking parameters from autoPyLoT-input file
    :type pickparameter: pylot.core.io.inputs.PylotParameter
    :param fig_dict: dictionary containing Matplotlib figures used for plotting results
    :type fig_dict: dict
    :return: dictionary containing iterative picks
    :rtype: dict
    """

    msg = '##################################################\n' \
          'autoPyLoT: Found {0} bad onsets at station(s) {1}, ' \
          'starting re-picking them ...'.format(len(badpicks), badpicks)
    print(msg)

    newpicks = {}
    for i in range(0, len(badpicks)):
        if len(badpicks[i][0]) > 4:
            Ppattern = '%s  ?    ?    ? P' % badpicks[i][0]
        elif len(badpicks[i][0]) == 4:
            Ppattern = '%s   ?    ?    ? P' % badpicks[i][0]
        elif len(badpicks[i][0]) < 4:
            Ppattern = '%s    ?    ?    ? P' % badpicks[i][0]
        nllocline = getPatternLine(NLLocfile, Ppattern)
        res = nllocline.split(None)[16]
        # get theoretical P-onset time from residuum
        badpicks[i][1] = picks[badpicks[i][0]]['P']['mpp'] - float(res)

        # get corresponding waveform stream
        msg = '##################################################\n' \
              'iteratepicker: Re-picking station {0}'.format(badpicks[i][0])
        print(msg)
        wf2pick = wf.select(station=badpicks[i][0])

        # modify some picking parameters
        pstart_old = pickparameter.get('pstart')
        pstop_old = pickparameter.get('pstop')
        sstop_old = pickparameter.get('sstop')
        pickwinP_old = pickparameter.get('pickwinP')
        Precalcwin_old = pickparameter.get('Precalcwin')
        noisefactor_old = pickparameter.get('noisefactor')
        zfac_old = pickparameter.get('zfac')
        twindows = pickparameter.get('tsnrz')
        tsafety = twindows[1]
        pstart = max([0, badpicks[i][1] - wf2pick[0].stats.starttime - pickparameter.get('tlta')])
        if abs(float(res)) <= tsafety / 2 or pstart == 0:
            print("iteratepicker: Small residuum, leave parameters unchanged for this phase!")
        else:
            pickparameter.setParam(pstart=pstart)
            pickparameter.setParam(pstop=pickparameter.get('pstart') + (pickparameter.get('Precalcwin')))
            pickparameter.setParam(sstop=pickparameter.get('sstop') / 2)
            pickparameter.setParam(pickwinP=pickparameter.get('pickwinP') / 2)
            pickparameter.setParam(Precalcwin=pickparameter.get('Precalcwin') / 2)
            pickparameter.setParam(noisefactor=1.0)
            pickparameter.setParam(zfac=1.0)

        print(
            "iteratepicker: The following picking parameters have been modified for iterative picking:")
        print(
            "pstart: %fs => %fs" % (pstart_old, pickparameter.get('pstart')))
        print(
            "pstop: %fs => %fs" % (pstop_old, pickparameter.get('pstop')))
        print(
            "sstop: %fs => %fs" % (sstop_old, pickparameter.get('sstop')))
        print("pickwinP: %fs => %fs" % (
            pickwinP_old, pickparameter.get('pickwinP')))
        print("Precalcwin: %fs => %fs" % (
            Precalcwin_old, pickparameter.get('Precalcwin')))
        print("noisefactor: %f => %f" % (
            noisefactor_old, pickparameter.get('noisefactor')))
        print("zfac: %f => %f" % (zfac_old, pickparameter.get('zfac')))

        # repick station
        newpicks, _ = autopickstation(wf2pick, pickparameter, fig_dict=fig_dict)

        # replace old dictionary with new one
        picks[badpicks[i][0]] = newpicks

        # reset temporary change of picking parameters
        print("iteratepicker: Resetting picking parameters ...")
        pickparameter.setParam(pstart=pstart_old)
        pickparameter.setParam(pstop=pstop_old)
        pickparameter.setParam(sstop=sstop_old)
        pickparameter.setParam(pickwinP=pickwinP_old)
        pickparameter.setParam(Precalcwin=Precalcwin_old)
        pickparameter.setParam(noisefactor=noisefactor_old)
        pickparameter.setParam(zfac=zfac_old)

    return picks
