#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Function to run automated picking algorithms using AIC,
HOS and AR prediction. Uses objects CharFuns and Picker and
function conglomerate utils.

:author: MAGS2 EP3 working group / Ludger Kueperkoch
"""

import matplotlib.pyplot as plt
import numpy as np
from pylot.core.pick.charfuns import CharacteristicFunction
from pylot.core.pick.charfuns import HOScf, AICcf, ARZcf, ARHcf, AR3Ccf
from pylot.core.pick.picker import AICPicker, PragPicker
from pylot.core.pick.utils import checksignallength, checkZ4S, earllatepicker, \
    getSNR, fmpicker, checkPonsets, wadaticheck, get_pickparams, get_quality_class
from pylot.core.util.utils import getPatternLine, gen_Pool,\
    real_Bool, identifyPhaseID, real_None, correct_iplot

from obspy.taup import TauPyModel


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

        if iplot is None or iplot == 'None' or iplot == 0:
            input_tuples.append((topick, param, apverbose, metadata, origin))
        if iplot > 0:
            all_onsets[station] = autopickstation(topick, param, verbose=apverbose,
                                                  iplot=iplot, fig_dict=fig_dict,
                                                  metadata=metadata, origin=origin)

    if iplot > 0:
        print('iPlot Flag active: NO MULTIPROCESSING possible.')
        return all_onsets

    # rename str for ncores in case ncores == 0 (use all cores)
    ncores_str = ncores if ncores != 0 else 'all available'

    print('Autopickstation: Distribute autopicking for {} '
          'stations on {} cores.'.format(len(input_tuples), ncores_str))

    pool = gen_Pool(ncores)
    result = pool.map(call_autopickstation, input_tuples)
    pool.close()

    for pick in result:
        if pick:
            station = pick['station']
            pick.pop('station')
            all_onsets[station] = pick

    # quality control
    # median check and jackknife on P-onset times
    jk_checked_onsets = checkPonsets(all_onsets, mdttolerance, jackfactor, 1, fig_dict_wadatijack)
    # check S-P times (Wadati)
    wadationsets = wadaticheck(jk_checked_onsets, wdttolerance, 1, fig_dict_wadatijack)
    return wadationsets


def call_autopickstation(input_tuple):
    """
    helper function used for multiprocessing
    :param input_tuple: contains all parameters used for autopicking
    :type input_tuple: tuple
    :return: dictionary containing P pick, S pick and station name
    :rtype: dict
    """
    wfstream, pickparam, verbose, metadata, origin = input_tuple
    # multiprocessing not possible with interactive plotting
    return autopickstation(wfstream, pickparam, verbose, iplot=0, metadata=metadata, origin=origin)


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

class PickingParameters(object):
    """
    Stores parameters used for picking a single station.
    """

    def __init__(self, *args, **kwargs):
        """
        Add dictionaries given as positional arguments and the keyword argument dictionary to the instance
        as attributes. Positional arguments with types differing from dict are ignored.
        """
        # add entries from dictionaries given as positional arguments
        for arg in args:
            if type(arg) == dict:
                self.add_params_from_dict(arg)
        # add values given as keyword arguments
        self.add_params_from_dict(kwargs)

    def add_params_from_dict(self, d):
        """
        Add all key-value pairs from dictionary d to the class namespace as attributes.
        :param d:
        :type d: dict
        :rtype: None
        """
        for key, value in d.items():
            setattr(self, key, value)

class PickingResults(dict):

    def __init__(self):
        # initialize output
        self.Pweight = 4  # weight for P onset
        self.Sweight = 4  # weight for S onset
        self.FM = 'N'  # first motion (polarity)
        self.SNRP = None  # signal-to-noise ratio of P onset
        self.SNRPdB = None  # signal-to-noise ratio of P onset [dB]
        self.SNRS = None  # signal-to-noise ratio of S onset
        self.SNRSdB = None  # signal-to-noise ratio of S onset [dB]
        self.mpickP = None  # most likely P onset
        self.lpickP = None  # latest possible P onset
        self.epickP = None  # earliest possible P onset
        self.mpickS = None  # most likely S onset
        self.lpickS = None  # latest possible S onset
        self.epickS = None  # earliest possible S onset
        self.Perror = None  # symmetrized picking error P onset
        self.Serror = None  # symmetrized picking error S onset

        #self.aicSflag = 0
        #self.aicPflag = 0
        #self.Pflag = 0
        self.Pmarker = []
        self.Ao = None  # Wood-Anderson peak-to-peak amplitude
        self.picker = 'auto'  # type of picks

        # these get added during the construction of the p pick results dictionary
        self.w0 = None
        self.fc = None
        self.Mw = None

    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]

class MissingTraceException(ValueError):
    """
    Used to indicate missing traces in a obspy.core.stream.Stream object
    """
    pass

class PickingFailedException(Exception):
    """
    Raised when picking fails due to missing values etc.
    """
    pass

class AutopickStation(object):

    def __init__(self, wfstream, pickparam, verbose, iplot, fig_dict, metadata, origin):
        # save given parameters
        self.wfstream = wfstream
        self.pickparam = pickparam
        self.verbose = verbose
        self.iplot = correct_iplot(iplot)
        self.fig_dict = real_None(fig_dict)
        self.metadata = metadata
        self.origin = origin

        # extract additional information
        pickparams = self.extract_pickparams(pickparam)
        self.p_params, self.s_params, self.first_motion_params, self.signal_length_params = pickparams
        self.results = PickingResults()
        self.channelorder = {'Z': 3, 'N': 1, 'E': 2} # TODO get this from the pylot preferences
        self.station_name = wfstream[0].stats.station
        self.network_name = wfstream[0].stats.network
        self.station_id = '{}.{}'.format(self.network_name, self.station_name)

        # save streams and traces
        # TODO: error handling of missing traces
        self.zstream, self.nstream, self.estream = self.get_components_from_waveformstream()
        self.ztrace = self.zstream[0]
        try:
            self.ntrace = self.nstream[0]
        except IndexError:
            self.nstream = self.estream
            self.ntrace = self.nstream[0]
        try:
            self.etrace = self.estream[0]
        except IndexError:
            self.estream = self.nstream
            self.etrace = self.estream[0]

        # default values used in old autopickstation function #TODO way for user to set those
        self.detrend_type = 'demean'
        self.filter_type = 'bandpass'
        self.zerophase = False
        self.taper_max_percentage = 0.05
        self.taper_type = 'hann'

        # initialize picking results
        self.p_results = PickingResults()
        self.s_results = PickingResults()

    def vprint(self, s):
        """Only print statement if verbose picking is set to true."""
        if self.verbose:
            print(s)

    def extract_pickparams(self, pickparam):
        """
        Get parameter names out of pickparam dictionary into PickingParameters objects and return them.
        :return: PickingParameters objects containing 1. p pick parameters, 2. s pick parameters, 3. first motion determinatiion
        parameters, 4. signal length parameters
        :rtype: (PickingParameters, PickingParameters, PickingParameters,  PickingParameters)
        """
        # Define names of all parameters in different groups
        p_parameter_names = 'algoP pstart pstop use_taup taup_model tlta tsnrz hosorder bpz1 bpz2 pickwinP aictsmooth tsmoothP ausP nfacP tpred1z tdet1z Parorder addnoise Precalcwin minAICPslope minAICPSNR timeerrorsP'.split(
            ' ')
        s_parameter_names = 'algoS sstart sstop bph1 bph2 tsnrh pickwinS tpred1h tdet1h tpred2h tdet2h Sarorder aictsmoothS tsmoothS ausS minAICSslope minAICSSNR Srecalcwin nfacS timeerrorsS zfac'.split(
            ' ')
        first_motion_names = 'minFMSNR fmpickwin minfmweight'.split(' ')
        signal_length_names = 'minsiglength minpercent noisefactor'.split(' ')
        # Get list of values from pickparam by name
        p_parameter_values = map(pickparam.get, p_parameter_names)
        s_parameter_values = map(pickparam.get, s_parameter_names)
        fm_parameter_values = map(pickparam.get, first_motion_names)
        sl_parameter_values = map(pickparam.get, signal_length_names)
        # construct dicts from names and values
        p_params = dict(zip(p_parameter_names, p_parameter_values))
        s_params = dict(zip(s_parameter_names, s_parameter_values))
        first_motion_params = dict(zip(first_motion_names, fm_parameter_values))
        signal_length_params = dict(zip(signal_length_names, sl_parameter_values))

        p_params['use_taup'] = real_Bool(p_params['use_taup'])

        return PickingParameters(p_params), PickingParameters(s_params), PickingParameters(first_motion_params), PickingParameters(signal_length_params)

    def get_components_from_waveformstream(self):
        """
        Splits waveformstream into multiple components zdat, ndat, edat. For traditional orientation (ZNE) these contain
        the vertical, north-south or east-west component. Otherwise they contain components numbered 123 with
        orientation diverging from the traditional orientation.
        :param waveformstream: Stream containing all three components for one station either by ZNE or 123 channel code
        (mixture of both options is handled as well)
        :type waveformstream: obspy.core.stream.Stream
        :return: Tuple containing (z waveform, n waveform, e waveform) selected by the given channels
        :rtype: (obspy.core.stream.Stream, obspy.core.stream.Stream, obspy.core.stream.Stream)
        """

        waveform_data = {}
        for key in self.channelorder:
            waveform_data[key] = self.wfstream.select(component=key) # try ZNE first
            if len(waveform_data[key]) == 0:
                waveform_data[key] = self.wfstream.select(component=str(self.channelorder[key])) # use 123 as second option
        return waveform_data['Z'], waveform_data['N'], waveform_data['E']

    def prepare_wfstream(self, wfstream, freqmin=None, freqmax=None):
        """
        Prepare a waveformstream for picking by applying detrending, filtering and tapering. Creates a copy of the
        waveform the leave the original unchanged.
        :param wfstream:
        :type wfstream:
        :param freqmin:
        :type freqmin:
        :param freqmax:
        :type freqmax:
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

        def create_arrivals(metadata, origin, station_id, taup_model):
            """
            Create List of arrival times for all phases for a given origin and station
            :param metadata: tuple containing metadata type string and Parser object read from inventory file
            :type metadata: tuple (str, ~obspy.io.xseed.parser.Parser)
            :param origin: list containing origin objects representing origins for all events
            :type origin: list(~obspy.core.event.origin)
            :param station_id: Station id with format NETWORKNAME.STATIONNAME
            :type station_id: str
            :param taup_model: Model name to use. See obspy.taup.tau.TauPyModel for options
            :type taup_model: str
            :return: List of Arrival objects
            :rtype: obspy.taup.tau.Arrivals
            :raises:
                AttributeError when no metadata or source origins is given
            """
            parser = metadata[1]
            station_coords = get_source_coords(parser, station_id)
            source_origin = origin[0]
            model = TauPyModel(taup_model)
            arrivals = model.get_travel_times_geo(source_depth_in_km=source_origin.depth,
                                                  source_latitude_in_deg=source_origin.latitude,
                                                  source_longitude_in_deg=source_origin.longitude,
                                                  receiver_latitude_in_deg=station_coords['latitude'],
                                                  receiver_longitude_in_deg=station_coords['longitude'])
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
            arrP, estFirstP = min([(arr, arr.time) for arr in phases['P']], key=lambda t: t[1])
            arrS, estFirstS = min([(arr, arr.time) for arr in phases['S']], key=lambda t: t[1])
            print('autopick: estimated first arrivals for P: {} s, S:{} s after event'
                  ' origin time using TauPy'.format(estFirstP, estFirstS))
            return estFirstP, estFirstS

        if self.p_params.use_taup is False and self.p_params.pstart < 0:
            # correct user mistake where a relative cuttime is selected (pstart < 0) but use of taupy is disabled/ has
            # not the required parameters
            self.p_params.pstart = 0
            return

        print('autopickstation: use_taup flag active.')
        # catch missing metadata or origin information. Onset calculation is stopped, given cuttimes are then used.
        if not self.metadata[1]:
            raise AttributeError('Warning: Could not use TauPy to estimate onsets as there are no metadata given.')
        if not self.origin:
            raise AttributeError('No source origins given!')

        arrivals = create_arrivals(self.metadata, self.origin, self.station_id, self.p_params.taup_model)
        estFirstP, estFirstS = first_PS_onsets(arrivals)
        # modifiy pstart and pstop relative to estimated first P arrival (relative to station time axis)
        self.p_params.pstart += (self.origin[0].time + estFirstP) - self.ztrace.stats.starttime
        self.p_params.pstop += (self.origin[0].time + estFirstP) - self.ztrace.stats.starttime
        print('autopick: CF calculation times respectively:'
              ' pstart: {} s, pstop: {} s'.format(self.p_params.pstart, self.p_params.pstop))
        # make sure pstart and pstop are inside the starttime/endtime of vertical trace
        self.p_params.pstart = max(self.p_params.pstart, 0)
        self.p_params.pstop = min(self.p_params.pstop, len(self.ztrace) * self.ztrace.stats.delta)

    def autopickstation(self):
        try:
            self.pick_p_phase()
        #TODO handle exceptions correctly (goal is to be compatible with old code first)
        # requires an overlook of what should be returned in case picking fails at various stages
        except MissingTraceException as mte:
            print(mte)
        except PickingFailedException as pfe:
            print(pfe)

        if self.estream is not None and self.nstream is not None and len(self.estream) > 0 and len(self.nstream) > 0 and self.p_results.Pweight is not None and self.p_results.Pweight < 4:
            try:
                self.pick_s_phase()
            # TODO: when an exception occurs, return picking results so far. This requires that the pick methods save their results in the instance member of PickingResults
            except MissingTraceException as mte:
                print(mte)
            except PickingFailedException as pfe:
                print(pfe)

        self.finish_picking()
        return {'P': self.p_results, 'S':self.s_results, 'station':self.ztrace.stats.station} #TODO method to format picking results as a dict correctly

    def finish_picking(self):

        # calculate "real" onset times, save them in PickingResults
        if self.p_results.lpickP is not None and self.p_results.lpickP == self.p_results.mpickP:
            self.p_results.lpickP += self.ztrace.stats.delta
        if self.p_results.epickP is not None and self.p_results.epickP == self.p_results.mpickP:
            self.p_results.epickP -= self.ztrace.stats.delta
        if self.p_results.mpickP is not None and self.p_results.epickP is not None and self.p_results.lpickP is not None:
            lpickP = self.ztrace.stats.starttime + self.p_results.lpickP
            epickP = self.ztrace.stats.starttime + self.p_results.epickP
            mpickP = self.ztrace.stats.starttime + self.p_results.mpickP
        else:
            # dummy values (start of seismic trace) in order to derive
            # theoretical onset times for iterative picking
            lpickP = self.ztrace.stats.starttime + self.p_params.timeerrorsP[3]
            epickP = self.ztrace.stats.starttime - self.p_params.timeerrorsP[3]
            mpickP = self.ztrace.stats.starttime

        # add picking results to PickingResults instance
        self.p_results.channel = self.ztrace.stats.channel
        self.p_results.network = self.ztrace.stats.network
        self.p_results.lpp = lpickP
        self.p_results.epp = epickP
        self.p_results.mpp = mpickP
        self.p_results.snrdb = self.p_results.SNRPdB
        self.p_results.snr = self.p_results.SNRP
        self.p_results.marked = self.p_results.Pmarker # TODO add equal implementation as in old function
        self.p_results.weight = self.p_results.Pweight
        self.p_results.fm = self.p_results.FM
        self.p_results.spe = self.p_results.Perror
        self.p_results.Mo = None  # TODO what is this parameter for?

        #
        #   S results
        #
        if self.etrace:
            hdat = self.etrace
        elif self.ntrace:
            hdat = self.ntrace

        # TODO picks are calculated with stats from E trace if it exists, otherwise stats from N trace is used.
        # What if the stats are different?
        if self.s_results.lpickS is not None and self.s_results.lpickS == self.s_results.mpickS:
            self.s_results.lpickS += hdat.stats.delta
        if self.s_results.epickS is not None and self.s_results.epickS == self.s_results.mpickS:
            self.s_results.epickS -= hdat.stats.delta
        if self.s_results.mpickS is not None and self.s_results.epickS is not None and self.s_results.lpickS is not None:
            lpickS = hdat.stats.starttime + self.s_results.lpickS
            epickS = hdat.stats.starttime + self.s_results.epickS
            mpickS = hdat.stats.starttime + self.s_results.mpickS
        else:
            # dummy values (start of seismic trace) in order to derive
            # theoretical onset times for iteratve picking
            lpickS = hdat.stats.starttime + self.s_params.timeerrorsS[3]
            epickS = hdat.stats.starttime - self.s_params.timeerrorsS[3]
            mpickS = hdat.stats.starttime

        self.s_results.channel = self.etrace.stats.channel
        self.s_results.network = self.etrace.stats.network
        self.s_results.lpp = lpickS
        self.s_results.epp = epickS
        self.s_results.mpp = mpickS
        self.s_results.spe = self.s_results.Serror
        self.s_results.snr = self.s_results.SNRS
        self.s_results.snrdb = self.s_results.SNRSdB
        self.s_results.weight = self.s_results.Sweight
        self.s_results.fm = None
        self.s_results.picker='auto'
        self.s_results.Ao = None

    def pick_p_qc1(self, aicpick, z_copy, tr_filt):
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

        fig, linecolor = get_fig_from_figdict(self.fig_dict, 'slength')
        if aicpick.getpick() is None:
            msg = "Bad initial (AIC) P-pick, skipping this onset!\nAIC-SNR={0}, AIC-Slope={1}counts/s\n " \
                  "(min. AIC-SNR={2}, min. AIC-Slope={3}counts/s)"
            msg = msg.format(aicpick.getSNR(), aicpick.getSlope(), self.p_params.minAICPSNR, self.p_params.minAICPslope)
            self.vprint(msg)
            return 0
        # Quality check initial pick with minimum signal length
        z_copy[0].data = tr_filt.data  # save filtered, tapered trace in z_copy stream object
        zne = z_copy
        if len(self.nstream) == 0 or len(self.estream) == 0:
            msg = 'One or more horizontal component(s) missing!\n' \
                  'Signal length only checked on vertical component!\n' \
                  'Decreasing minsiglengh from {0} to {1}'\
                  .format(self.signal_length_params.minsiglength, self.signal_length_params.minsiglength / 2)
            self.vprint(msg)
            minsiglength = self.signal_length_params.minsiglength / 2
        else:
            # filter, taper other traces as well since signal length is compared on all traces
            trH1_filt, _ = self.prepare_wfstream(self.estream, freqmin=self.s_params.bph1[0], freqmax=self.s_params.bph1[1])
            trH2_filt, _ = self.prepare_wfstream(self.nstream, freqmin=self.s_params.bph1[0], freqmax=self.s_params.bph1[1])
            zne += trH1_filt
            zne += trH2_filt
            minsiglength = self.signal_length_params.minsiglength
        Pflag = checksignallength(zne, aicpick.getpick(), self.p_params.tsnrz, minsiglength,
                                  self.signal_length_params.noisefactor, self.signal_length_params.minpercent,
                                  self.iplot, fig, linecolor)
        if Pflag == 0:
            self.p_results.Pmarker = 'shortsignallength'
            self.p_results.Pweight = 9
            return 0

        if self.nstream == self.estream:
            # todo: old implementation skipped this test if one component was misisng, why not use one component?
            msg = 'One or more horizontal components missing!\n Skipping control function checkZ4S.'
            self.vprint(msg)
            return 1

        if self.iplot > 1: fig, linecolor = get_fig_from_figdict(self.fig_dict, 'checkZ4s')
        Pflag = checkZ4S(zne, aicpick.getpick(), self.s_params.zfac, self.p_params.tsnrz[2], self.iplot, fig, linecolor)
        if Pflag == 0:
            self.p_results.Pmarker = 'SinsteadP'
            self.p_results.Pweight = 9
            return 0
        return 1

    def pick_p_phase(self):
        """
        Pick p phase, return results
        :return: P pick results
        :rtype: PickingResults
        :raises:
            MissingTraceException: If vertical trace is missing.
        """
        if not self.zstream or self.zstream is None:
            raise MissingTraceException('No z-component found for station {}'.format(self.station_name))

        msg = '##################################################\nautopickstation:' \
              ' Working on P onset of station {station}\nFiltering vertical ' \
              'trace ...\n{data}'.format(station=self.station_name, data=str(self.zstream))
        self.vprint(msg)

        tr_filt, z_copy = self.prepare_wfstream(self.zstream, self.p_params.bpz1[0], self.p_params.bpz1[1])
        # save filtered trace in instance for later plotting
        self.tr_filt_z = tr_filt
        try:
            # modify pstart, pstop to be around theoretical onset if taupy should be used, else does nothing
            self.modify_starttimes_taupy()
        except AttributeError as ae:
            print(ae)
        except MissingTraceException as mte:
            print(mte)
        cuttimes = [self.p_params.pstart, self.p_params.pstop]

        # calculate first CF
        if self.p_params.algoP == 'HOS':
            cf1 = HOScf(z_copy, cuttimes, self.p_params.tlta, self.p_params.hosorder)
        elif self.p_params.algoP == 'ARZ':
            cf1 = ARZcf(z_copy, cuttimes, self.p_params.tpred1z, self.p_params.Parorder, self.p_params.tdet1z,
                        self.p_params.addnoise)
        else:
            cf1 = None
        # save cf1 for plotting
        self.cf1 = cf1
        assert isinstance(cf1, CharacteristicFunction), 'cf2 is not set ' \
                                                        'correctly: maybe the algorithm name ({algoP}) is ' \
                                                        'corrupted'.format(algoP=self.p_params.algoP)
        # AICcf needs stream object -> build it
        tr_aic = tr_filt.copy()
        tr_aic.data = cf1.getCF()
        z_copy[0].data = tr_aic.data
        aiccf = AICcf(z_copy, cuttimes)
        # get preliminary onset time from AIC-CF
        fig, linecolor = get_fig_from_figdict(self.fig_dict, 'aicFig')
        aicpick = AICPicker(aiccf, self.p_params.tsnrz, self.p_params.pickwinP, self.iplot,
                            Tsmooth=self.p_params.aictsmooth, fig=fig, linecolor=linecolor)
        # add pstart and pstop to aic plot
        if fig:
            for ax in fig.axes:
                ax.vlines(self.p_params.pstart, ax.get_ylim()[0], ax.get_ylim()[1], color='c', linestyles='dashed', label='P start')
                ax.vlines(self.p_params.pstop, ax.get_ylim()[0], ax.get_ylim()[1], color='c', linestyles='dashed', label='P stop')
                ax.legend(loc=1)

        Pflag = self.pick_p_qc1(aicpick, z_copy, tr_filt)
        # go on with processing if AIC onset passes quality control
        slope = aicpick.getSlope()
        if not slope: slope = 0
        if slope >= self.p_params.minAICPslope and aicpick.getSNR() >= self.p_params.minAICPSNR and Pflag == 1:
            aicPflag = 1
            msg = 'AIC P-pick passes quality control: Slope: {0} counts/s, ' \
                  'SNR: {1}\nGo on with refined picking ...\n' \
                  'autopickstation: re-filtering vertical trace ' \
                  '...'.format(aicpick.getSlope(), aicpick.getSNR())
            self.vprint(msg)
            # refilter waveform with larger bandpass
            tr_filt, z_copy = self.prepare_wfstream(self.zstream, freqmin=self.p_params.bpz2[0], freqmax=self.p_params.bpz2[1])
            cuttimes2 = [round(max([aicpick.getpick() - self.p_params.Precalcwin, 0])),
                         round(min([len(self.ztrace.data) * self.ztrace.stats.delta,
                                    aicpick.getpick() + self.p_params.Precalcwin]))]
            if self.p_params.algoP == 'HOS':
                cf2 = HOScf(z_copy, cuttimes2, self.p_params.tlta, self.p_params.hosorder)
            elif self.p_params.algoP == 'ARZ':
                cf2 = ARZcf(z_copy, cuttimes2, self.p_params.tpred1z, self.p_params.Parorder, self.p_params.tdet1z, self.p_params.addnoise)
            else:
                cf2 = None
            # save cf2 for plotting
            self.cf2 = cf2
            # get refined onset time from CF2
            assert isinstance(cf2, CharacteristicFunction), 'cf2 is not set ' \
                                                            'correctly: maybe the algorithm name ({algoP}) is ' \
                                                            'corrupted'.format(algoP=self.p_params.algoP)
            fig, linecolor = get_fig_from_figdict(self.fig_dict, 'refPpick')
            refPpick = PragPicker(cf2, self.p_params.tsnrz, self.p_params.pickwinP, self.iplot, self.p_params.ausP,
                                  self.p_params.tsmoothP, aicpick.getpick(), fig, linecolor)
            self.p_results.mpickP = refPpick.getpick()
            if self.p_results.mpickP is not None:
                # quality assessment, get earliest/latest pick and symmetrized uncertainty
                fig, linecolor = get_fig_from_figdict(self.fig_dict, 'el_Ppick')
                self.p_results.epickP, self.p_results.lpickP, self.p_results.Perror = earllatepicker(z_copy, self.p_params.nfacP, self.p_params.tsnrz, self.p_results.mpickP,
                                                        self.iplot, fig=fig, linecolor=linecolor)
                self.p_results.SNRP, self.p_results.SNRPdB, self.p_results.Pnoiselevel = getSNR(z_copy, self.p_params.tsnrz, self.p_results.mpickP)

                # weight P-onset using symmetric error
                self.p_results.Pweight = get_quality_class(self.p_results.Perror, self.p_params.timeerrorsP)
                if self.p_results.Pweight <= self.first_motion_params.minfmweight and self.p_results.SNRP >= self.first_motion_params.minFMSNR:
                    fig, linecolor = get_fig_from_figdict(self.fig_dict, 'fm_picker')
                    self.p_results.FM = fmpicker(self.zstream, z_copy, self.first_motion_params.fmpickwin, self.p_results.mpickP, self.iplot,
                                  fig, linecolor)
                else:
                    self.p_results.FM = 'N'
                msg = "autopickstation: P-weight: {0}, " \
                      "SNR: {1}, SNR[dB]: {2}, Polarity: {3}".format(self.p_results.Pweight, self.p_results.SNRP, self.p_results.SNRPdB, self.p_results.FM)
                print(msg)
                msg = 'autopickstation: Refined P-Pick: {} s | P-Error: {} s'.format(self.p_results.mpickP, self.p_results.Perror)
                print(msg)
                Sflag = 1

                self.p_results.aicpick = aicpick
            else:
                msg = 'Bad initial (AIC) P-pick, skipping this onset!\n' \
                      'AIC-SNR={0}, AIC-Slope={1}counts/s\n' \
                      '(min. AIC-SNR={2}, ' \
                      'min. AIC-Slope={3}counts/s)'.format(aicpick.getSNR(), aicpick.getSlope(),
                                                           self.p_params.minAICPSNR, self.p_params.minAICPslope)
                self.vprint(msg)
                Sflag = 0
        else:
            #todo add why did picking fail, which should be saved in the pick dictionary
            raise PickingFailedException("Why did it fail")

    def pick_s_phase(self):

        def check_existence_ne_traces():
            """
            Checks if both N and E stream contain traces and copy over the trace to the missing stream if only
            one trace is missing
            :rtype: None
            :raises:
                MissingTraceException when both traces are missing
            """
            # check existence of vertical traces #Maybe do this in __init__?
            nstream_exists, estream_exists = bool(len(self.nstream)), bool(len(self.estream))
            if not all((nstream_exists, estream_exists)):
                msg = 'Go on picking S onset ...\n' \
                      '##################################################\n' \
                      'Only one horizontal component available!\n' \
                      'ARH prediction requires at least 2 components!\n' \
                      'Copying existing horizontal component ...'
                self.vprint(msg)
                if estream_exists == 0:
                    self.estream = self.nstream
                elif nstream_exists == 0:
                    self.nstream = self.estream
            else:
                raise MissingTraceException("N and E stream missing, can't pick S phase.")


        #check_existence_ne_traces()
        # determine time window for calculating CF after P onset
        start = round(max(self.p_results.mpickP + self.s_params.sstart, 0))
        stop = round(min([
            self.p_results.mpickP + self.s_params.sstop,
            self.etrace.stats.endtime - self.etrace.stats.starttime,
            self.ntrace.stats.endtime - self.ntrace.stats.starttime
        ]))
        cuttimesh = (start, stop)

        if cuttimesh[1] <= cuttimesh[0]:
            pickSonset = False
            raise PickingFailedException('Cut window for horizontal phases too small! Will not pick S onsets.')

        # prepare traces for picking by filtering, taper
        if self.s_params.algoS == 'ARH':
            hdat = self.nstream.copy() + self.estream.copy()
            trH1_filt, _ = self.prepare_wfstream(self.estream, freqmin=self.s_params.bph1[0], freqmax=self.s_params.bph1[1])
            trH2_filt, _ = self.prepare_wfstream(self.nstream, freqmin=self.s_params.bph1[0], freqmax=self.s_params.bph1[1])
            h_copy = hdat.copy()
            h_copy[0].data = trH1_filt.data
            h_copy[1].data = trH2_filt.data
        if self.s_params.algoS == 'AR3':
            hdat = self.zstream.copy() + self.estream.copy() + self.nstream.copy()
            trH1_filt, _ = self.prepare_wfstream(self.zstream, freqmin=self.s_params.bph1[0], freqmax=self.s_params.bph1[1])
            trH2_filt, _ = self.prepare_wfstream(self.estream, freqmin=self.s_params.bph1[0], freqmax=self.s_params.bph1[1])
            trH3_filt, _ = self.prepare_wfstream(self.nstream, freqmin=self.s_params.bph1[0], freqmax=self.s_params.bph1[1])
            h_copy = hdat.copy()
            h_copy[0].data = trH1_filt.data
            h_copy[1].data = trH2_filt.data
            h_copy[2].data = trH3_filt.data

        # calculate inital CF based on autoregression
        if self.s_params.algoS == 'ARH':
            arhcf1 = ARHcf(h_copy, cuttimesh, self.s_params.tpred1h, self.s_params.Sarorder, self.s_params.tdet1h, self.p_params.addnoise)
        elif self.s_params.algoS == 'AR3':
            arhcf1 = AR3Ccf(h_copy, cuttimesh, self.s_params.tpred1h, self.s_params.Sarorder, self.s_params.tdet1h, self.p_params.addnoise)

        tr_arhaic = trH1_filt.copy()
        tr_arhaic.data = arhcf1.getCF()
        h_copy[0].data = tr_arhaic.data

        # calculate AIC cf
        haiccf = AICcf(h_copy, cuttimesh)

        # get preliminary onset time from AIC cf
        fig, linecolor = get_fig_from_figdict(self.fig_dict, 'aicARHfig')
        aicarhpick = AICPicker(haiccf, self.s_params.tsnrh, self.s_params.pickwinS, self.iplot, Tsmooth=self.s_params.aictsmoothS, fig=fig, linecolor=linecolor)

        # go on with processing if AIC onset passes quality control
        slope = aicarhpick.getSlope()

        if not slope:
            slope = 0

        if slope >= self.s_params.minAICSslope and aicarhpick.getSNR() >= self.s_params.minAICSSNR and aicarhpick.getpick() is not None:
            aicSflag = 1
            msg = 'AIC S-pick passes quality control: Slope: {0} counts/s, ' \
                  'SNR: {1}\nGo on with refined picking ...\n' \
                  'autopickstation: re-filtering horizontal traces ' \
                  '...'.format(aicarhpick.getSlope(), aicarhpick.getSNR())
            self.vprint(msg)

            # recalculate cf from refiltered trace in vicinity of initial onset
            start = round(aicarhpick.getpick() - self.s_params.Srecalcwin)
            stop = round(aicarhpick.getpick() + self.s_params.Srecalcwin)
            cuttimesh2 = (start, stop)

            # refilter waveform with larger bandpass
            if self.s_params.algoS == 'ARH':
                trH1_filt, _ = self.prepare_wfstream(self.estream, freqmin=self.s_params.bph2[0], freqmax=self.s_params.bph2[1])
                trH2_filt, _ = self.prepare_wfstream(self.nstream, freqmin=self.s_params.bph2[0], freqmax=self.s_params.bph2[1])
                h_copy = hdat.copy()
                h_copy[0].data = trH1_filt.data
                h_copy[1].data = trH2_filt.data
            elif self.s_params.algoS == 'AR3':
                trH1_filt, _ = self.prepare_wfstream(self.zstream, freqmin=self.s_params.bph2[0], freqmax=self.s_params.bph2[1])
                trH2_filt, _ = self.prepare_wfstream(self.estream, freqmin=self.s_params.bph2[0], freqmax=self.s_params.bph2[1])
                trH3_filt, _ = self.prepare_wfstream(self.nstream, freqmin=self.s_params.bph2[0], freqmax=self.s_params.bph2[1])
                h_copy = hdat.copy()
                h_copy[0].data = trH1_filt.data
                h_copy[1].data = trH2_filt.data
                h_copy[2].data = trH3_filt.data

            # calculate second cd
            if self.s_params.algoS == 'ARH':
                arhcf2 = ARHcf(h_copy, cuttimesh2, self.s_params.tpred2h, self.s_params.Sarorder, self.s_params.tdet2h, self.p_params.addnoise)
            elif self.s_params.algoS == 'AR3':
                arhcf2 = AR3Ccf(h_copy, cuttimesh2, self.s_params.tpred2h, self.s_params.Sarorder, self.s_params.tdet2h, self.p_params.addnoise)

            # get refined onset time from CF2
            fig, linecolor = get_fig_from_figdict(self.fig_dict, 'refSpick')
            refSpick = PragPicker(arhcf2, self.s_params.tsnrh, self.s_params.pickwinS, self.iplot, self.s_params.ausS, self.s_params.tsmoothS, aicarhpick.getpick(), fig, linecolor)
            self.s_results.mpickS = refSpick.getpick()

            if self.s_results.mpickS is not None:
                # quality assessment
                # get earliest/latest possible pick and symmetrized uncertainty
                h_copy[0].data = trH1_filt.data
                if self.iplot:
                    fig, linecolor = get_fig_from_figdict(self.fig_dict, 'el_Slpick')
                epickS1, lpickS1, Serror1 = earllatepicker(h_copy, self.s_params.nfacS, self.s_params.tsnrh, self.s_results.mpickS, self.iplot, fig=fig, linecolor=linecolor)

                h_copy[0].data = trH2_filt.data
                if self.iplot:
                    fig, linecolor = get_fig_from_figdict(self.fig_dict, 'el_S2pick')
                else:
                    # why is it set to empty here? DA
                    linecolor = ''
                epickS2, lpickS2, Serror2 = earllatepicker(h_copy, self.s_params.nfacS, self.s_params.tsnrh, self.s_results.mpickS, self.iplot, fig=fig, linecolor=linecolor)

                if epickS1 is not None and epickS2 is not None:
                    if self.s_params.algoS == 'ARH':
                        # get earliest pick of both earliest possible picks
                        epick = [epickS1, epickS2]
                        lpick = [lpickS1, lpickS2]
                        pickerr = [Serror1, Serror2]
                        ipick = np.argmin(epick)
                    if self.s_params.algoS == 'AR3':
                        epickS3, lpickS3, Serror3 = earllatepicker(h_copy, self.s_params.nfacS, self.s_params.tsnrh, self.s_results.mpickS, self.iplot)
                        # get earliest of all three picks
                        epick = [epickS1, epickS2, epickS3]
                        lpick = [lpickS1, lpickS2, lpickS3]
                        pickerr = [Serror1, Serror2, Serror3]

                        if epickS3 is not None:
                            ipick = np.argmin(epick)
                        else:
                            ipick = np.argmin([epickS1, epickS2])
                    self.s_results.epickS = epick[ipick]
                    self.s_results.lpickS = lpick[ipick]
                    self.s_results.Serror = pickerr[ipick]

                    msg = 'autopickstation: Refined S-Pick: {} s | S-Error: {} s'.format(self.s_results.mpickS, self.s_results.Serror)
                    print(msg)

                    # get SNR
                    self.s_results.SNRS, self.s_results.SNRSdB, self.s_results.Snoiselevel = getSNR(h_copy, self.s_params.tsnrh, self.s_results.mpickS)

                    self.s_results.Sweight = get_quality_class(self.s_results.Serror, self.s_params.timeerrorsS)

                    print('autopickstation: S-weight: {0}, SNR: {1}, '
                          'SNR[dB]: {2}\n'
                          '##################################################'
                          ''.format(self.s_results.Sweight, self.s_results.SNRS, self.s_results.SNRSdB))

        else:
            print('autopickstation: No horizontal component data available or '
                  'bad P onset, skipping S picking!')

        ### plotting ###
        # TODO finish converting the plotting code
        if self.iplot > 0:
            # plot vertical trace
            if self.fig_dict is None:
                fig = plt.figure()
                plt_flag = 1
                linecolor = 'k'
            else:
                fig, linecolor = get_fig_from_figdict(self.fig_dict, 'mainFig')
            fig._tight = True
            ax1 = fig.add_subplot(311)
            # create time axis
            tdata = np.linspace(start=0., stop=self.ztrace.stats.npts*self.tr_filt_z.stats.delta, num=self.tr_filt_z.stats.npts)
            # plot filtered waveform of z trace
            ax1.plot(tdata, self.tr_filt_z.data / max(self.tr_filt_z.data), color=linecolor, linewidth=0.7, label='Data')
            if self.p_results.Pweight < 4:
                # plot first HOS/ARZ cf of z trace
                ax1.plot(self.cf1.getTimeArray(), self.cf1.getCF / max(self.cf1.getCF()), color='b', label='CF1')
                if self.p_results.aicPflag == 1:
                    ax1.plot(self.cf2.getTimeArray(), self.cf2.getCF() / max(self.cf2.getCF()), color='m', label='CF2')
                    ax1.plot([self.p_results.aicpick.getpick(), self.p_results.aicpick.getpick()], [-1, 1], 'r', label='Initial P Onset')


def get_fig_from_figdict(figdict, figkey):
    """
    Helper method to extract a figure by name from dictionary
    :param figdict:
    :type figdict: dict
    :param figkey:
    :type figkey: str
    :return:
    :rtype:
    """
    if figdict is None:
        return None, None
    fig = figdict.get(figkey, None)
    linecolor = figdict.get('plot_style', 'k')
    linecolor = linecolor['linecolor']['rgba_mpl']
    return fig, linecolor


def autopickstation(wfstream, pickparam, verbose=False, iplot=0, fig_dict=None, metadata=None, origin=None):
    station = AutopickStation(wfstream, pickparam, verbose, iplot, fig_dict, metadata, origin)
    return station.autopickstation()


def nautopickstation(wfstream, pickparam, verbose=False,
                    iplot=0, fig_dict=None, metadata=None, origin=None):
    """
    picks a single station
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
    :return: dictionary containing P pick, S pick and station name
    :rtype: dict
    """

    # declaring pickparam variables (only for convenience)
    # read your pylot.in for details!
    plt_flag = 0

    # get picking parameter dictionaries
    p_params, s_params, first_motion_params, signal_length_params = get_pickparams(pickparam)

    # initialize output
    Pweight = 4  # weight for P onset
    Sweight = 4  # weight for S onset
    FM = 'N'  # first motion (polarity)
    SNRP = None  # signal-to-noise ratio of P onset
    SNRPdB = None  # signal-to-noise ratio of P onset [dB]
    SNRS = None  # signal-to-noise ratio of S onset
    SNRSdB = None  # signal-to-noise ratio of S onset [dB]
    mpickP = None  # most likely P onset
    lpickP = None  # latest possible P onset
    epickP = None  # earliest possible P onset
    mpickS = None  # most likely S onset
    lpickS = None  # latest possible S onset
    epickS = None  # earliest possible S onset
    Perror = None  # symmetrized picking error P onset
    Serror = None  # symmetrized picking error S onset

    aicSflag = 0
    aicPflag = 0
    Pflag = 0
    Sflag = 0
    Pmarker = []
    Ao = None  # Wood-Anderson peak-to-peak amplitude
    picker = 'auto'  # type of picks

    def get_components_from_waveformstream(waveformstream):
        """
        Splits waveformstream into multiple components zdat, ndat, edat. For traditional orientation (ZNE) these contain
        the vertical, north-south or east-west component. Otherwise they contain components numbered 123 with
        orientation diverging from the traditional orientation.
        :param waveformstream: Stream containing all three components for one station either by ZNE or 123 channel code
        (mixture of both options is handled as well)
        :type waveformstream: obspy.core.stream.Stream
        :return: Tuple containing (z waveform, n waveform, e waveform) selected by the given channels
        :rtype: (obspy.core.stream.Stream, obspy.core.stream.Stream, obspy.core.stream.Stream)
        """

        #TODO: get this order from the pylot preferences
        channelorder_default = {'Z': 3, 'N': 1, 'E': 2}
        waveform_data = {}
        for key in channelorder_default:
            waveform_data[key] = waveformstream.select(component=key) # try ZNE first
            if len(waveform_data[key]) == 0:
                waveform_data[key] = waveformstream.select(component=str(channelorder_default[key])) # use 123 as second option
        return waveform_data['Z'], waveform_data['N'], waveform_data['E']


    def prepare_wfstream_component(wfstream, detrend_type='demean', filter_type='bandpass', freqmin=None, freqmax=None, zerophase=False, taper_max_percentage=0.05, taper_type='hann'):
        """
        Prepare a waveformstream for picking by applying detrending, filtering and tapering. Creates a copy of the
        waveform the leave the original unchanged.
        :param wfstream:
        :type wfstream:
        :param detrend_type:
        :type detrend_type:
        :param filter_type:
        :type filter_type:
        :param freqmin:
        :type freqmin:
        :param freqmax:
        :type freqmax:
        :param zerophase:
        :type zerophase:
        :param taper_max_percentage:
        :type taper_max_percentage:
        :param taper_type:
        :type taper_type:
        :return: Tuple containing the changed waveform stream and the first trace of the stream
        :rtype: (obspy.core.stream.Stream, obspy.core.trace.Trace)
        """
        wfstream_copy = wfstream.copy()
        trace_copy = wfstream[0].copy()
        trace_copy.detrend(type=detrend_type)
        trace_copy.filter(filter_type, freqmin=freqmin, freqmax=freqmax, zerophase=zerophase)
        trace_copy.taper(max_percentage=taper_max_percentage, type=taper_type)
        wfstream_copy[0].data = trace_copy.data
        return wfstream_copy, trace_copy

    # split components
    zdat, ndat, edat = get_components_from_waveformstream(wfstream)

    if not zdat:
        print('No z-component found for station {}. STOP'.format(wfstream[0].stats.station))
        return

    if p_params['algoP'] == 'HOS' or p_params['algoP'] == 'ARZ' and zdat is not None:
        msg = '##################################################\nautopickstation:' \
              ' Working on P onset of station {station}\nFiltering vertical ' \
              'trace ...\n{data}'.format(station=wfstream[0].stats.station, data=str(zdat))
        if verbose: print(msg)
        z_copy, tr_filt = prepare_wfstream_component(zdat, freqmin=p_params['bpz1'][0], freqmax=p_params['bpz1'][1])
        ##############################################################
        # check length of waveform and compare with cut times

        # for global seismology: use tau-p method for estimating travel times (needs source and station coords.)
        # if not given: sets Lc to infinity to use full stream
        if p_params['use_taup'] is True:
            Lc = np.inf
            print('autopickstation: use_taup flag active.')
            if not metadata[1]:
                print('Warning: Could not use TauPy to estimate onsets as there are no metadata given.')
            else:
                station_id = wfstream[0].get_id()
                parser = metadata[1]
                station_coords = get_source_coords(parser, station_id)
                if station_coords and origin:
                    source_origin = origin[0]
                    model = TauPyModel(p_params['taup_model'])
                    arrivals = model.get_travel_times_geo(
                        source_origin.depth,
                        source_origin.latitude,
                        source_origin.longitude,
                        station_coords['latitude'],
                        station_coords['longitude']
                    )
                    phases = {'P': [],
                              'S': []}
                    for arr in arrivals:
                        phases[identifyPhaseID(arr.phase.name)].append(arr)

                    # get first P and S onsets from arrivals list
                    arrP, estFirstP = min([(arr, arr.time) for arr in phases['P']], key=lambda t: t[1])
                    arrS, estFirstS = min([(arr, arr.time) for arr in phases['S']], key=lambda t: t[1])
                    print('autopick: estimated first arrivals for P: {} s, S:{} s after event'
                          ' origin time using TauPy'.format(estFirstP, estFirstS))

                    # modifiy pstart and pstop relative to estimated first P arrival (relative to station time axis)
                    p_params['pstart'] += (source_origin.time + estFirstP) - zdat[0].stats.starttime
                    p_params['pstop']+= (source_origin.time + estFirstP) - zdat[0].stats.starttime
                    print('autopick: CF calculation times respectively:'
                          ' pstart: {} s, pstop: {} s'.format(p_params['pstart'], p_params['pstop']))
                elif not origin:
                    print('No source origins given!')

        # make sure pstart and pstop are inside zdat[0]
        pstart = max(p_params['pstart'], 0)
        pstop = min(p_params['pstop'], len(zdat[0])*zdat[0].stats.delta)

        if p_params['use_taup'] is False or origin:
            Lc = p_params['pstop'] - p_params['pstart']

        Lwf = zdat[0].stats.endtime - zdat[0].stats.starttime
        if not Lwf > 0:
            print('autopickstation: empty trace! Return!')
            return

        Ldiff = Lwf - abs(Lc)
        if Ldiff < 0 or pstop <= pstart:
            msg = 'autopickstation: Cutting times are too large for actual ' \
                  'waveform!\nUsing entire waveform instead!'
            if verbose: print(msg)
            pstart = 0
            pstop = len(zdat[0].data) * zdat[0].stats.delta
        cuttimes = [pstart, pstop]
        cf1 = None
        if p_params['algoP'] == 'HOS':
            # calculate HOS-CF using subclass HOScf of class
            # CharacteristicFunction
            cf1 = HOScf(z_copy, cuttimes, p_params['tlta'], p_params['hosorder'])  # instance of HOScf
        elif p_params['algoP'] == 'ARZ':
            # calculate ARZ-CF using subclass ARZcf of class
            # CharcteristicFunction
            cf1 = ARZcf(z_copy, cuttimes, p_params['tpred1z'], p_params['Parorder'], p_params['tdet1z'],
                        p_params['addnoise'])  # instance of ARZcf
        ##############################################################
        # calculate AIC-HOS-CF using subclass AICcf of class
        # CharacteristicFunction
        # class needs stream object => build it
        assert isinstance(cf1, CharacteristicFunction), 'cf2 is not set ' \
                                                        'correctly: maybe the algorithm name ({algoP}) is ' \
                                                        'corrupted'.format(algoP=p_params['algoP'])
        tr_aic = tr_filt.copy()
        tr_aic.data = cf1.getCF()
        z_copy[0].data = tr_aic.data
        aiccf = AICcf(z_copy, cuttimes)  # instance of AICcf
        ##############################################################
        # get preliminary onset time from AIC-HOS-CF using subclass AICPicker
        # of class AutoPicking
        key = 'aicFig'
        if fig_dict:
            fig = fig_dict[key]
            linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
        else:
            fig = None
            linecolor = 'k'
        aicpick = AICPicker(aiccf, p_params['tsnrz'], p_params['pickwinP'], iplot, Tsmooth=p_params['aictsmooth'],
                            fig=fig, linecolor=linecolor)
        # add pstart and pstop to aic plot
        if fig:
            for ax in fig.axes:
                ax.vlines(pstart, ax.get_ylim()[0], ax.get_ylim()[1], color='c', linestyles='dashed', label='P start')
                ax.vlines(pstop, ax.get_ylim()[0], ax.get_ylim()[1], color='c', linestyles='dashed', label='P stop')
                ax.legend(loc=1)
        ##############################################################
        if aicpick.getpick() is not None:
            # check signal length to detect spuriously picked noise peaks
            # use all available components to avoid skipping correct picks
            # on vertical traces with weak P coda
            z_copy[0].data = tr_filt.data
            zne = z_copy
            if len(ndat) == 0 or len(edat) == 0:
                msg = 'One or more horizontal component(s) missing!\n' \
                      'Signal length only checked on vertical component!\n' \
                      'Decreasing minsiglengh from {0} to {1}' \
                      .format(signal_length_params['minsiglength'], signal_length_params['minsiglength'] / 2)
                if verbose: print(msg)
                key = 'slength'
                if fig_dict:
                    fig = fig_dict[key]
                    linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                else:
                    fig = None
                    linecolor = 'k'
                Pflag = checksignallength(zne, aicpick.getpick(), p_params['tsnrz'],
                                          signal_length_params['minsiglength'] / 2,
                                          signal_length_params['noisefactor'], signal_length_params['minpercent'], iplot,
                                          fig, linecolor)
            else:
                trH1_filt, _ = prepare_wfstream_component(edat, freqmin=s_params['bph1'][0], freqmax=s_params['bph1'][1])
                trH2_filt, _ = prepare_wfstream_component(ndat, freqmin=s_params['bph1'][0], freqmax=s_params['bph1'][1])
                zne += trH1_filt
                zne += trH2_filt
                if fig_dict:
                    fig = fig_dict['slength']
                    linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                else:
                    fig = None
                    linecolor = 'k'
                Pflag = checksignallength(zne, aicpick.getpick(), p_params['tsnrz'],
                                          signal_length_params['minsiglength'],
                                          signal_length_params['noisefactor'], signal_length_params['minpercent'], iplot,
                                          fig, linecolor)

            if Pflag == 1:
                # check for spuriously picked S onset
                # both horizontal traces needed
                if len(ndat) == 0 or len(edat) == 0:
                    msg = 'One or more horizontal components missing!\n' \
                          'Skipping control function checkZ4S.'
                    if verbose: print(msg)
                else:
                    if iplot > 1:
                        if fig_dict:
                            fig = fig_dict['checkZ4s']
                            linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                        else:
                            fig = None
                            linecolor = 'k'
                    Pflag = checkZ4S(zne, aicpick.getpick(), s_params['zfac'],
                                     p_params['tsnrz'][2], iplot, fig, linecolor)
                    if Pflag == 0:
                        Pmarker = 'SinsteadP'
                        Pweight = 9
            else:
                Pmarker = 'shortsignallength'
                Pweight = 9
        ##############################################################
        # go on with processing if AIC onset passes quality control
        slope = aicpick.getSlope()
        if not slope:
            slope = 0
        if slope >= p_params['minAICPslope'] and aicpick.getSNR() >= p_params['minAICPSNR'] and Pflag == 1:
            aicPflag = 1
            msg = 'AIC P-pick passes quality control: Slope: {0} counts/s, ' \
                  'SNR: {1}\nGo on with refined picking ...\n' \
                  'autopickstation: re-filtering vertical trace ' \
                  '...'.format(aicpick.getSlope(), aicpick.getSNR())
            if verbose: print(msg)
            # re-filter waveform with larger bandpass
            z_copy, tr_filt = prepare_wfstream_component(zdat, freqmin=p_params['bpz2'][0], freqmax=p_params['bpz2'][1])
            #############################################################
            # re-calculate CF from re-filtered trace in vicinity of initial
            # onset
            cuttimes2 = [round(max([aicpick.getpick() - p_params['Precalcwin'], 0])),
                         round(min([len(zdat[0].data) * zdat[0].stats.delta,
                                    aicpick.getpick() + p_params['Precalcwin']]))]
            cf2 = None
            if p_params['algoP'] == 'HOS':
                # calculate HOS-CF using subclass HOScf of class
                # CharacteristicFunction
                cf2 = HOScf(z_copy, cuttimes2, p_params['tlta'],
                            p_params['hosorder'])  # instance of HOScf
            elif p_params['algoP'] == 'ARZ':
                # calculate ARZ-CF using subclass ARZcf of class
                # CharcteristicFunction
                cf2 = ARZcf(z_copy, cuttimes2, p_params['tpred1z'], p_params['Parorder'], p_params['tdet1z'],
                            p_params['addnoise'])  # instance of ARZcf
            ##############################################################
            # get refined onset time from CF2 using class Picker
            assert isinstance(cf2, CharacteristicFunction), 'cf2 is not set ' \
                                                            'correctly: maybe the algorithm name ({algoP}) is ' \
                                                            'corrupted'.format(algoP=p_params['algoP'])
            if fig_dict:
                fig = fig_dict['refPpick']
                linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
            else:
                fig = None
                linecolor = 'k'
            refPpick = PragPicker(cf2, p_params['tsnrz'], p_params['pickwinP'], iplot, p_params['ausP'],
                                  p_params['tsmoothP'], aicpick.getpick(), fig, linecolor)
            mpickP = refPpick.getpick()
            #############################################################
            if mpickP is not None:
                # quality assessment
                # get earliest/latest possible pick and symmetrized uncertainty
                if iplot:
                    if fig_dict:
                        fig = fig_dict['el_Ppick']
                        linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                    else:
                        fig = None
                        linecolor = 'k'
                epickP, lpickP, Perror = earllatepicker(z_copy, p_params['nfacP'], p_params['tsnrz'],
                                                        mpickP, iplot, fig=fig,
                                                        linecolor=linecolor)

                # get SNR
                SNRP, SNRPdB, Pnoiselevel = getSNR(z_copy, p_params['tsnrz'], mpickP)

                # weight P-onset using symmetric error
                Pweight = get_quality_class(Perror, p_params['timeerrorsP'])

                ##############################################################
                # get first motion of P onset
                # certain quality required
                if Pweight <= first_motion_params['minfmweight'] and SNRP >= first_motion_params['minFMSNR']:
                    if iplot:
                        if fig_dict:
                            fig = fig_dict['fm_picker']
                            linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                        else:
                            fig = None
                        FM = fmpicker(zdat, z_copy, first_motion_params['fmpickwin'], mpickP, iplot, fig, linecolor)
                    else:
                        FM = fmpicker(zdat, z_copy, first_motion_params['fmpickwin'], mpickP, iplot)
                else:
                    FM = 'N'

                msg = "autopickstation: P-weight: {0}, " \
                      "SNR: {1}, SNR[dB]: {2}, Polarity: {3}".format(Pweight, SNRP, SNRPdB, FM)
                print(msg)
                msg = 'autopickstation: Refined P-Pick: {} s | P-Error: {} s'.format(mpickP, Perror)
                print(msg)
                Sflag = 1

        else:
            msg = 'Bad initial (AIC) P-pick, skipping this onset!\n' \
                  'AIC-SNR={0}, AIC-Slope={1}counts/s\n' \
                  '(min. AIC-SNR={2}, ' \
                  'min. AIC-Slope={3}counts/s)'.format(aicpick.getSNR(),
                                                       aicpick.getSlope(),
                                                       p_params['minAICPSNR'],
                                                       p_params['minAICPslope'])
            if verbose: print(msg)
            Sflag = 0

    else:
        print('autopickstation: No vertical component data available!, '
              'Skipping station!')

    if ((len(edat) > 0 and len(ndat) == 0) or (len(ndat) > 0 and len(edat) == 0)) and Pweight < 4:
        msg = 'Go on picking S onset ...\n' \
              '##################################################\n' \
              'Only one horizontal component available!\n' \
              'ARH prediction requires at least 2 components!\n' \
              'Copying existing horizontal component ...'
        if verbose: print(msg)

        # check which component is missing
        if len(edat) == 0:
            edat = ndat
        else:
            ndat = edat

    pickSonset = (edat is not None and ndat is not None and len(edat) > 0 and len(
                  ndat) > 0 and Pweight < 4)

    if pickSonset:
        # determine time window for calculating CF after P onset
        cuttimesh = [
            round(max([mpickP + s_params['sstart'], 0])),  # MP MP relative time axis
            round(min([
                mpickP + s_params['sstop'],
                edat[0].stats.endtime-edat[0].stats.starttime,
                ndat[0].stats.endtime-ndat[0].stats.starttime
            ]))
        ]

        if not cuttimesh[1] >= cuttimesh[0]:
            print('Cut window for horizontal phases too small! Will not pick S onsets.')
            pickSonset = False

    if pickSonset:
        msg = 'Go on picking S onset ...\n' \
              '##################################################\n' \
              'Working on S onset of station {0}\nFiltering horizontal ' \
              'traces ...'.format(edat[0].stats.station)
        if verbose: print(msg)

        if s_params['algoS'] == 'ARH':
            # re-create stream object including both horizontal components
            hdat = edat.copy()
            hdat += ndat
            h_copy = hdat.copy()
            # filter and taper data
            trH1_filt = hdat[0].copy()
            trH2_filt = hdat[1].copy()
            trH1_filt.detrend(type='demean')
            trH2_filt.detrend(type='demean')
            trH1_filt.filter('bandpass', freqmin=s_params['bph1'][0], freqmax=s_params['bph1'][1],
                             zerophase=False)
            trH2_filt.filter('bandpass', freqmin=s_params['bph1'][0], freqmax=s_params['bph1'][1],
                             zerophase=False)
            trH1_filt.taper(max_percentage=0.05, type='hann')
            trH2_filt.taper(max_percentage=0.05, type='hann')
            h_copy[0].data = trH1_filt.data
            h_copy[1].data = trH2_filt.data
        elif s_params['algoS'] == 'AR3':
            # re-create stream object including all components
            hdat = zdat.copy()
            hdat += edat
            hdat += ndat
            h_copy = hdat.copy()
            # filter and taper data
            trH1_filt = hdat[0].copy()
            trH2_filt = hdat[1].copy()
            trH3_filt = hdat[2].copy()
            trH1_filt.detrend(type='demean')
            trH2_filt.detrend(type='demean')
            trH3_filt.detrend(type='demean')
            trH1_filt.filter('bandpass', freqmin=s_params['bph1'][0], freqmax=s_params['bph1'][1],
                             zerophase=False)
            trH2_filt.filter('bandpass', freqmin=s_params['bph1'][0], freqmax=s_params['bph1'][1],
                             zerophase=False)
            trH3_filt.filter('bandpass', freqmin=s_params['bph1'][0], freqmax=s_params['bph1'][1],
                             zerophase=False)
            trH1_filt.taper(max_percentage=0.05, type='hann')
            trH2_filt.taper(max_percentage=0.05, type='hann')
            trH3_filt.taper(max_percentage=0.05, type='hann')
            h_copy[0].data = trH1_filt.data
            h_copy[1].data = trH2_filt.data
            h_copy[2].data = trH3_filt.data
        ##############################################################
        if s_params['algoS'] == 'ARH':
            # calculate ARH-CF using subclass ARHcf of class
            # CharcteristicFunction
            arhcf1 = ARHcf(h_copy, cuttimesh, s_params['tpred1h'], s_params['Sarorder'], s_params['tdet1h'],
                           p_params['addnoise'])  # instance of ARHcf
        elif s_params['algoS'] == 'AR3':
            # calculate ARH-CF using subclass AR3cf of class
            # CharcteristicFunction
            arhcf1 = AR3Ccf(h_copy, cuttimesh, s_params['tpred1h'], s_params['Sarorder'], s_params['tdet1h'],
                            p_params['addnoise'])  # instance of ARHcf
        ##############################################################
        # calculate AIC-ARH-CF using subclass AICcf of class
        # CharacteristicFunction
        # class needs stream object => build it
        tr_arhaic = trH1_filt.copy()
        tr_arhaic.data = arhcf1.getCF()
        h_copy[0].data = tr_arhaic.data
        # calculate ARH-AIC-CF
        haiccf = AICcf(h_copy, cuttimesh)  # instance of AICcf
        ##############################################################
        # get prelimenary onset time from AIC-HOS-CF using subclass AICPicker
        # of class AutoPicking
        if fig_dict:
            fig = fig_dict['aicARHfig']
            linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
        else:
            fig = None
            linecolor = 'k'
        aicarhpick = AICPicker(haiccf, s_params['tsnrh'], s_params['pickwinS'], iplot, None,
                               s_params['aictsmoothS'], fig=fig, linecolor=linecolor)
        ###############################################################
        # go on with processing if AIC onset passes quality control
        slope = aicarhpick.getSlope()
        if not slope:
            slope = 0
        if (slope >= s_params['minAICSslope'] and
            aicarhpick.getSNR() >= s_params['minAICSSNR'] and aicarhpick.getpick() is not None):
            aicSflag = 1
            msg = 'AIC S-pick passes quality control: Slope: {0} counts/s, ' \
                  'SNR: {1}\nGo on with refined picking ...\n' \
                  'autopickstation: re-filtering horizontal traces ' \
                  '...'.format(aicarhpick.getSlope(), aicarhpick.getSNR())
            if verbose: print(msg)
            # re-calculate CF from re-filtered trace in vicinity of initial
            # onset
            cuttimesh2 = [round(aicarhpick.getpick() - s_params['Srecalcwin']),
                          round(aicarhpick.getpick() + s_params['Srecalcwin'])]
            # re-filter waveform with larger bandpass
            h_copy = hdat.copy()
            # filter and taper data
            if s_params['algoS']== 'ARH':
                trH1_filt = hdat[0].copy()
                trH2_filt = hdat[1].copy()
                trH1_filt.detrend(type='demean')
                trH2_filt.detrend(type='demean')
                trH1_filt.filter('bandpass', freqmin=s_params['bph2'][0], freqmax=s_params['bph2'][1],
                                 zerophase=False)
                trH2_filt.filter('bandpass', freqmin=s_params['bph2'][0], freqmax=s_params['bph2'][1],
                                 zerophase=False)
                trH1_filt.taper(max_percentage=0.05, type='hann')
                trH2_filt.taper(max_percentage=0.05, type='hann')
                h_copy[0].data = trH1_filt.data
                h_copy[1].data = trH2_filt.data
                #############################################################
                arhcf2 = ARHcf(h_copy, cuttimesh2, s_params['tpred2h'], s_params['Sarorder'], s_params['tdet2h'],
                               p_params['addnoise'])  # instance of ARHcf
            elif s_params['algoS'] == 'AR3':
                trH1_filt = hdat[0].copy()
                trH2_filt = hdat[1].copy()
                trH3_filt = hdat[2].copy()
                trH1_filt.detrend(type='demean')
                trH2_filt.detrend(type='demean')
                trH3_filt.detrend(type='demean')
                trH1_filt.filter('bandpass', freqmin=s_params['bph2'][0], freqmax=s_params['bph2'][1],
                                 zerophase=False)
                trH2_filt.filter('bandpass', freqmin=s_params['bph2'][0], freqmax=s_params['bph2'][1],
                                 zerophase=False)
                trH3_filt.filter('bandpass', freqmin=s_params['bph2'][0], freqmax=s_params['bph2'][1],
                                 zerophase=False)
                trH1_filt.taper(max_percentage=0.05, type='hann')
                trH2_filt.taper(max_percentage=0.05, type='hann')
                trH3_filt.taper(max_percentage=0.05, type='hann')
                h_copy[0].data = trH1_filt.data
                h_copy[1].data = trH2_filt.data
                h_copy[2].data = trH3_filt.data
                #############################################################
                arhcf2 = AR3Ccf(h_copy, cuttimesh2, s_params['tpred2h'], s_params['Sarorder'], s_params['tdet2h'],
                                p_params['addnoise'])  # instance of ARHcf

            # get refined onset time from CF2 using class Picker
            if fig_dict:
                fig = fig_dict['refSpick']
                linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
            else:
                fig = None
                linecolor = 'k'
            refSpick = PragPicker(arhcf2, s_params['tsnrh'], s_params['pickwinS'], iplot, s_params['ausS'],
                                  s_params['tsmoothS'], aicarhpick.getpick(), fig, linecolor)
            mpickS = refSpick.getpick()
            #############################################################
            if mpickS is not None:
                # quality assessment
                # get earliest/latest possible pick and symmetrized uncertainty
                h_copy[0].data = trH1_filt.data
                if iplot:
                    if fig_dict:
                        fig = fig_dict['el_S1pick']
                        linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                    else:
                        fig = None
                        linecolor = 'k'
                    epickS1, lpickS1, Serror1 = earllatepicker(h_copy, s_params['nfacS'], s_params['tsnrh'], mpickS,
                                                               iplot, fig=fig, linecolor=linecolor)
                else:
                    epickS1, lpickS1, Serror1 = earllatepicker(h_copy, s_params['nfacS'], s_params['tsnrh'], mpickS, iplot)

                h_copy[0].data = trH2_filt.data
                if iplot:
                    if fig_dict:
                        fig = fig_dict['el_S2pick']
                        linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                    else:
                        fig = None
                        linecolor = ''
                    epickS2, lpickS2, Serror2 = earllatepicker(h_copy, s_params['nfacS'], s_params['tsnrh'],  mpickS,
                                                               iplot, fig=fig, linecolor=linecolor)
                else:
                    epickS2, lpickS2, Serror2 = earllatepicker(h_copy, s_params['nfacS'], s_params['tsnrh'], mpickS, iplot)
                if epickS1 is not None and epickS2 is not None:
                    if s_params['algoS'] == 'ARH':
                        # get earliest pick of both earliest possible picks
                        epick = [epickS1, epickS2]
                        lpick = [lpickS1, lpickS2]
                        pickerr = [Serror1, Serror2]
                        if epickS1 is None and epickS2 is not None:
                            ipick = 1
                        elif epickS1 is not None and epickS2 is None:
                            ipick = 0
                        elif epickS1 is not None and epickS2 is not None:
                            ipick = np.argmin([epickS1, epickS2])
                    elif s_params['algoS'] == 'AR3':
                        [epickS3, lpickS3, Serror3] = earllatepicker(h_copy, s_params['nfacS'], s_params['tsnrh'],
                                                                     mpickS, iplot)
                        # get earliest pick of all three picks
                        epick = [epickS1, epickS2, epickS3]
                        lpick = [lpickS1, lpickS2, lpickS3]
                        pickerr = [Serror1, Serror2, Serror3]
                        if epickS1 is None and epickS2 is not None \
                                and epickS3 is not None:
                            ipick = np.argmin([epickS2, epickS3])
                        elif epickS1 is not None and epickS2 is None \
                                and epickS3 is not None:
                            ipick = np.argmin([epickS2, epickS3])
                        elif epickS1 is not None and epickS2 is not None \
                                and epickS3 is None:
                            ipick = np.argmin([epickS1, epickS2])
                        elif epickS1 is not None and epickS2 is not None \
                                and epickS3 is not None:
                            ipick = np.argmin([epickS1, epickS2, epickS3])

                    epickS = epick[ipick]
                    lpickS = lpick[ipick]
                    Serror = pickerr[ipick]

                    msg = 'autopickstation: Refined S-Pick: {} s | S-Error: {} s'.format(mpickS, Serror)
                    print(msg)

                    # get SNR
                    [SNRS, SNRSdB, Snoiselevel] = getSNR(h_copy, s_params['tsnrh'], mpickS)

                    # weight S-onset using symmetric error
                    if Serror <= s_params['timeerrorsS'][0]:
                        Sweight = 0
                    elif s_params['timeerrorsS'][0] < Serror <= s_params['timeerrorsS'][1]:
                        Sweight = 1
                    elif Perror > s_params['timeerrorsS'][1] and Serror <= s_params['timeerrorsS'][2]:
                        Sweight = 2
                    elif s_params['timeerrorsS'][2] < Serror <= s_params['timeerrorsS'][3]:
                        Sweight = 3
                    elif Serror > s_params['timeerrorsS'][3]:
                        Sweight = 4

                    print('autopickstation: S-weight: {0}, SNR: {1}, '
                          'SNR[dB]: {2}\n'
                          '##################################################'
                          ''.format(Sweight, SNRS, SNRSdB))
                ################################################################
                # get Wood-Anderson peak-to-peak amplitude
                # initialize Data object
                # re-create stream object including both horizontal components
                hdat = edat.copy()
                hdat += ndat
        else:
            msg = 'Bad initial (AIC) S-pick, skipping this onset!\n' \
                  'AIC-SNR={0}, AIC-Slope={1}counts/s\n' \
                  '(min. AIC-SNR={2}, ' \
                  'min. AIC-Slope={3}counts/s)\n' \
                  '##################################################' \
                  ''.format(aicarhpick.getSNR(), aicarhpick.getSlope(), s_params['minAICSSNR'], s_params['minAICSslope'])
            if verbose: print(msg)

            ############################################################
            # get Wood-Anderson peak-to-peak amplitude
            # initialize Data object
            # re-create stream object including both horizontal components
            hdat = edat.copy()
            hdat += ndat
 
    else:
        print('autopickstation: No horizontal component data available or '
              'bad P onset, skipping S picking!')

    ##############################################################
    try:
        iplot = int(iplot)
    except ValueError:
        if iplot is True or iplot == 'True':
            iplot = 2
        else:
            iplot = 0

    if iplot > 0:
        # plot vertical trace
        if fig_dict is None or fig_dict == 'None':
            fig = plt.figure()
            plt_flag = 1
            linecolor = 'k'
        else:
            fig = fig_dict['mainFig']
            linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
        fig._tight = True
        ax1 = fig.add_subplot(311)
        tdata = np.arange(0, zdat[0].stats.npts / tr_filt.stats.sampling_rate,
                          tr_filt.stats.delta)
        # check equal length of arrays, sometimes they are different!?
        wfldiff = len(tr_filt.data) - len(tdata)
        if wfldiff < 0:
            tdata = tdata[0:len(tdata) - abs(wfldiff)]
        ax1.plot(tdata, tr_filt.data / max(tr_filt.data), color=linecolor, linewidth=0.7, label='Data')
        if Pweight < 4:
            ax1.plot(cf1.getTimeArray(), cf1.getCF() / max(cf1.getCF()),
                     'b', label='CF1')
            if aicPflag == 1:
                ax1.plot(cf2.getTimeArray(),
                         cf2.getCF() / max(cf2.getCF()), 'm', label='CF2')
                ax1.plot([aicpick.getpick(), aicpick.getpick()], [-1, 1],
                         'r', label='Initial P Onset')
                ax1.plot([aicpick.getpick() - 0.5, aicpick.getpick() + 0.5],
                         [1, 1], 'r')
                ax1.plot([aicpick.getpick() - 0.5, aicpick.getpick() + 0.5],
                         [-1, -1], 'r')
                ax1.plot([refPpick.getpick(), refPpick.getpick()],
                         [-1.3, 1.3], 'r', linewidth=2, label='Final P Pick')
                ax1.plot([refPpick.getpick() - 0.5, refPpick.getpick() + 0.5],
                         [1.3, 1.3], 'r', linewidth=2)
                ax1.plot([refPpick.getpick() - 0.5, refPpick.getpick() + 0.5],
                         [-1.3, -1.3], 'r', linewidth=2)
                ax1.plot([lpickP, lpickP], [-1.1, 1.1], 'r--', label='lpp')
                ax1.plot([epickP, epickP], [-1.1, 1.1], 'r--', label='epp')
                ax1.set_title('%s, %s, P Weight=%d, SNR=%7.2f, SNR[dB]=%7.2f '
                              'Polarity: %s' % (tr_filt.stats.station,
                                                tr_filt.stats.channel,
                                                Pweight,
                                                SNRP,
                                                SNRPdB,
                                                FM))
            else:
                ax1.set_title('%s, P Weight=%d, SNR=None, '
                              'SNRdB=None' % (tr_filt.stats.channel, Pweight))
        else:
            ax1.set_title('%s, %s, P Weight=%d' % (tr_filt.stats.station,
                                                   tr_filt.stats.channel,
                                                   Pweight))
        ax1.legend(loc=1)
        ax1.set_yticks([])
        ax1.set_ylim([-1.5, 1.5])
        ax1.set_ylabel('Normalized Counts')
        # fig.suptitle(tr_filt.stats.starttime)
        try:
            len(edat[0])
        except:
            edat = ndat
        try:
            len(ndat[0])
        except:
            ndat = edat
        if len(edat[0]) > 1 and len(ndat[0]) > 1 and Sflag == 1:
            # plot horizontal traces
            ax2 = fig.add_subplot(3, 1, 2, sharex=ax1)
            th1data = np.arange(0,
                                trH1_filt.stats.npts /
                                trH1_filt.stats.sampling_rate,
                                trH1_filt.stats.delta)
            # check equal length of arrays, sometimes they are different!?
            wfldiff = len(trH1_filt.data) - len(th1data)
            if wfldiff < 0:
                th1data = th1data[0:len(th1data) - abs(wfldiff)]
            ax2.plot(th1data, trH1_filt.data / max(trH1_filt.data), color=linecolor, linewidth=0.7, label='Data')
            if Pweight < 4:
                ax2.plot(arhcf1.getTimeArray(),
                         arhcf1.getCF() / max(arhcf1.getCF()), 'b', label='CF1')
                if aicSflag == 1 and Sweight < 4:
                    ax2.plot(arhcf2.getTimeArray(),
                             arhcf2.getCF() / max(arhcf2.getCF()), 'm', label='CF2')
                    ax2.plot(
                        [aicarhpick.getpick(), aicarhpick.getpick()],
                        [-1, 1], 'g', label='Initial S Onset')
                    ax2.plot(
                        [aicarhpick.getpick() - 0.5,
                         aicarhpick.getpick() + 0.5],
                        [1, 1], 'g')
                    ax2.plot(
                        [aicarhpick.getpick() - 0.5,
                         aicarhpick.getpick() + 0.5],
                        [-1, -1], 'g')
                    ax2.plot([refSpick.getpick(), refSpick.getpick()],
                             [-1.3, 1.3], 'g', linewidth=2, label='Final S Pick')
                    ax2.plot(
                        [refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                        [1.3, 1.3], 'g', linewidth=2)
                    ax2.plot(
                        [refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                        [-1.3, -1.3], 'g', linewidth=2)
                    ax2.plot([lpickS, lpickS], [-1.1, 1.1], 'g--', label='lpp')
                    ax2.plot([epickS, epickS], [-1.1, 1.1], 'g--', label='epp')
                    ax2.set_title('%s, S Weight=%d, SNR=%7.2f, SNR[dB]=%7.2f' % (
                        trH1_filt.stats.channel,
                        Sweight, SNRS, SNRSdB))
                else:
                    ax2.set_title('%s, S Weight=%d, SNR=None, SNRdB=None' % (
                        trH1_filt.stats.channel, Sweight))
            ax2.legend(loc=1)
            ax2.set_yticks([])
            ax2.set_ylim([-1.5, 1.5])
            ax2.set_ylabel('Normalized Counts')
            # fig.suptitle(trH1_filt.stats.starttime)

            ax3 = fig.add_subplot(3, 1, 3, sharex=ax1)
            th2data = np.arange(0,
                                trH2_filt.stats.npts /
                                trH2_filt.stats.sampling_rate,
                                trH2_filt.stats.delta)
            # check equal length of arrays, sometimes they are different!?
            wfldiff = len(trH2_filt.data) - len(th2data)
            if wfldiff < 0:
                th2data = th2data[0:len(th2data) - abs(wfldiff)]
            ax3.plot(th2data, trH2_filt.data / max(trH2_filt.data), color=linecolor, linewidth=0.7, label='Data')
            if Pweight < 4:
                p22, = ax3.plot(arhcf1.getTimeArray(),
                                arhcf1.getCF() / max(arhcf1.getCF()), 'b', label='CF1')
                if aicSflag == 1:
                    ax3.plot(arhcf2.getTimeArray(),
                             arhcf2.getCF() / max(arhcf2.getCF()), 'm', label='CF2')
                    ax3.plot(
                        [aicarhpick.getpick(), aicarhpick.getpick()],
                        [-1, 1], 'g', label='Initial S Onset')
                    ax3.plot(
                        [aicarhpick.getpick() - 0.5,
                         aicarhpick.getpick() + 0.5],
                        [1, 1], 'g')
                    ax3.plot(
                        [aicarhpick.getpick() - 0.5,
                         aicarhpick.getpick() + 0.5],
                        [-1, -1], 'g')
                    ax3.plot([refSpick.getpick(), refSpick.getpick()],
                             [-1.3, 1.3], 'g', linewidth=2, label='Final S Pick')
                    ax3.plot(
                        [refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                        [1.3, 1.3], 'g', linewidth=2)
                    ax3.plot(
                        [refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                        [-1.3, -1.3], 'g', linewidth=2)
                    ax3.plot([lpickS, lpickS], [-1.1, 1.1], 'g--', label='lpp')
                    ax3.plot([epickS, epickS], [-1.1, 1.1], 'g--', label='epp')
            ax3.legend(loc=1)
            ax3.set_yticks([])
            ax3.set_ylim([-1.5, 1.5])
            ax3.set_xlabel('Time [s] after %s' % tr_filt.stats.starttime)
            ax3.set_ylabel('Normalized Counts')
            ax3.set_title(trH2_filt.stats.channel)
            if plt_flag == 1:
                fig.show()
                try:
                    input()
                except SyntaxError:
                    pass
                plt.close(fig)
    ##########################################################################
    # calculate "real" onset times
    if lpickP is not None and lpickP == mpickP:
        lpickP += zdat[0].stats.delta
    if epickP is not None and epickP == mpickP:
        epickP -= zdat[0].stats.delta
    if mpickP is not None and epickP is not None and lpickP is not None:
        lpickP = zdat[0].stats.starttime + lpickP
        epickP = zdat[0].stats.starttime + epickP
        mpickP = zdat[0].stats.starttime + mpickP
    else:
        # dummy values (start of seismic trace) in order to derive
        # theoretical onset times for iteratve picking
        lpickP = zdat[0].stats.starttime + p_params['timeerrorsP'][3]
        epickP = zdat[0].stats.starttime - p_params['timeerrorsP'][3]
        mpickP = zdat[0].stats.starttime

    if edat:
        hdat = edat[0]
    elif ndat:
        hdat = ndat[0]
    else:
        return

    if lpickS is not None and lpickS == mpickS:
        lpickS += hdat.stats.delta
    if epickS is not None and epickS == mpickS:
        epickS -= hdat.stats.delta
    if mpickS is not None and epickS is not None and lpickS is not None:
        lpickS = hdat.stats.starttime + lpickS
        epickS = hdat.stats.starttime + epickS
        mpickS = hdat.stats.starttime + mpickS
    else:
        # dummy values (start of seismic trace) in order to derive
        # theoretical onset times for iteratve picking
        lpickS = hdat.stats.starttime + s_params['timeerrorsS'][3]
        epickS = hdat.stats.starttime - s_params['timeerrorsS'][3]
        mpickS = hdat.stats.starttime

    # create dictionary
    # for P phase
    ccode = zdat[0].stats.channel
    ncode = zdat[0].stats.network
    ppick = dict(channel=ccode, network=ncode, lpp=lpickP, epp=epickP, mpp=mpickP, spe=Perror, snr=SNRP,
                 snrdb=SNRPdB, weight=Pweight, fm=FM, w0=None, fc=None, Mo=None,
                 Mw=None, picker=picker, marked=Pmarker)
    # add S phase
    ccode = hdat.stats.channel
    ncode = hdat.stats.network
    spick = dict(channel=ccode, network=ncode, lpp=lpickS, epp=epickS, mpp=mpickS, spe=Serror, snr=SNRS,
                 snrdb=SNRSdB, weight=Sweight, fm=None, picker=picker, Ao=Ao)
    # merge picks into returning dictionary
    picks = dict(P=ppick, S=spick, station=zdat[0].stats.station)
    return picks


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
        newpicks = autopickstation(wf2pick, pickparameter, fig_dict=fig_dict)

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
