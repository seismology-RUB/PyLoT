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
from pylot.core.io.data import Data
from pylot.core.io.inputs import PylotParameter
from pylot.core.pick.charfuns import CharacteristicFunction
from pylot.core.pick.charfuns import HOScf, AICcf, ARZcf, ARHcf, AR3Ccf
from pylot.core.pick.picker import AICPicker, PragPicker
from pylot.core.pick.utils import checksignallength, checkZ4S, earllatepicker, \
    getSNR, fmpicker, checkPonsets, wadaticheck
from pylot.core.util.utils import getPatternLine, gen_Pool,\
    real_Bool, identifyPhaseID

from obspy.taup import TauPyModel


def autopickevent(data, param, iplot=0, fig_dict=None, fig_dict_wadatijack=None, ncores=0, metadata=None, origin=None):
    stations = []
    all_onsets = {}
    input_tuples = []
    try:
        iplot = int(iplot)
    except:
        if iplot == True or iplot == 'True':
           iplot = 2
        else:
           iplot = 0


    # get some parameters for quality control from
    # parameter input file (usually autoPyLoT.in).
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

        if iplot == None or iplot == 'None' or iplot == 0:
            input_tuples.append((topick, param, apverbose, metadata, origin))
        if iplot > 0:
            all_onsets[station] = autopickstation(topick, param, verbose=apverbose,
                                                  iplot=iplot, fig_dict=fig_dict,
                                                  metadata=metadata, origin=origin)

    if iplot > 0:
        print('iPlot Flag active: NO MULTIPROCESSING possible.')
        return all_onsets

    print('Autopickstation: Distribute autopicking for {} '
          'stations on {} cores.'.format(len(input_tuples), ncores))

    pool = gen_Pool(ncores)
    result = pool.map(call_autopickstation, input_tuples)
    pool.close()

    for pick in result:
        if pick:
            station = pick['station']
            pick.pop('station')
            all_onsets[station] = pick

    #return all_onsets

    # quality control
    # median check and jackknife on P-onset times
    jk_checked_onsets = checkPonsets(all_onsets, mdttolerance, jackfactor, 1, fig_dict_wadatijack)
    #return jk_checked_onsets
    # check S-P times (Wadati)
    wadationsets = wadaticheck(jk_checked_onsets, wdttolerance, 1, fig_dict_wadatijack)
    return wadationsets


def call_autopickstation(input_tuple):
    wfstream, pickparam, verbose, metadata, origin = input_tuple
    # multiprocessing not possible with interactive plotting
    return autopickstation(wfstream, pickparam, verbose, iplot=0, metadata=metadata, origin=origin)


def get_source_coords(parser, station_id):
    station_coords = None
    try:
        station_coords = parser.get_coordinates(station_id)
    except Exception as e:
        print('Could not get source coordinates for station {}: {}'.format(station_id, e))
    return station_coords


def autopickstation(wfstream, pickparam, verbose=False,
                    iplot=0, fig_dict=None, metadata=None, origin=None):
    """
    :param wfstream: `~obspy.core.stream.Stream`  containing waveform
    :type wfstream: obspy.core.stream.Stream

    :param pickparam: container of picking parameters from input file,
           usually autoPyLoT.in
    :type pickparam: PylotParameter
    :param verbose:
    :type verbose: bool

    """

    # declaring pickparam variables (only for convenience)
    # read your autoPyLoT.in for details!
    plt_flag = 0

    # special parameters for P picking
    algoP = pickparam.get('algoP')
    pstart = pickparam.get('pstart')
    pstop = pickparam.get('pstop')
    thosmw = pickparam.get('tlta')
    tsnrz = pickparam.get('tsnrz')
    hosorder = pickparam.get('hosorder')
    bpz1 = pickparam.get('bpz1')
    bpz2 = pickparam.get('bpz2')
    pickwinP = pickparam.get('pickwinP')
    aictsmoothP = pickparam.get('aictsmooth')
    tsmoothP = pickparam.get('tsmoothP')
    ausP = pickparam.get('ausP')
    nfacP = pickparam.get('nfacP')
    tpred1z = pickparam.get('tpred1z')
    tdet1z = pickparam.get('tdet1z')
    Parorder = pickparam.get('Parorder')
    addnoise = pickparam.get('addnoise')
    Precalcwin = pickparam.get('Precalcwin')
    minAICPslope = pickparam.get('minAICPslope')
    minAICPSNR = pickparam.get('minAICPSNR')
    timeerrorsP = pickparam.get('timeerrorsP')
    # special parameters for S picking
    algoS = pickparam.get('algoS')
    sstart = pickparam.get('sstart')
    sstop = pickparam.get('sstop')
    use_taup = real_Bool(pickparam.get('use_taup'))
    taup_model = pickparam.get('taup_model')
    bph1 = pickparam.get('bph1')
    bph2 = pickparam.get('bph2')
    tsnrh = pickparam.get('tsnrh')
    pickwinS = pickparam.get('pickwinS')
    tpred1h = pickparam.get('tpred1h')
    tdet1h = pickparam.get('tdet1h')
    tpred2h = pickparam.get('tpred2h')
    tdet2h = pickparam.get('tdet2h')
    Sarorder = pickparam.get('Sarorder')
    aictsmoothS = pickparam.get('aictsmoothS')
    tsmoothS = pickparam.get('tsmoothS')
    ausS = pickparam.get('ausS')
    minAICSslope = pickparam.get('minAICSslope')
    minAICSSNR = pickparam.get('minAICSSNR')
    Srecalcwin = pickparam.get('Srecalcwin')
    nfacS = pickparam.get('nfacS')
    timeerrorsS = pickparam.get('timeerrorsS')
    # parameters for first-motion determination
    minFMSNR = pickparam.get('minFMSNR')
    fmpickwin = pickparam.get('fmpickwin')
    minfmweight = pickparam.get('minfmweight')
    # parameters for checking signal length
    minsiglength = pickparam.get('minsiglength')
    minpercent = pickparam.get('minpercent')
    nfacsl = pickparam.get('noisefactor')
    # parameter to check for spuriously picked S onset
    zfac = pickparam.get('zfac')
    # path to inventory-, dataless- or resp-files

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

    # split components
    zdat = wfstream.select(component="Z")
    if len(zdat) == 0:  # check for other components
        zdat = wfstream.select(component="3")
    edat = wfstream.select(component="E")
    if len(edat) == 0:  # check for other components
        edat = wfstream.select(component="2")
    ndat = wfstream.select(component="N")
    if len(ndat) == 0:  # check for other components
        ndat = wfstream.select(component="1")

    if not zdat:
        print('No z-component found for station {}. STOP'.format(wfstream[0].stats.station))
        return

    if algoP == 'HOS' or algoP == 'ARZ' and zdat is not None:
        msg = '##################################################\nautopickstation:' \
              ' Working on P onset of station {station}\nFiltering vertical ' \
              'trace ...\n{data}'.format(station=wfstream[0].stats.station,
                                         data=str(zdat))
        if verbose: print(msg)
        z_copy = zdat.copy()
        # filter and taper data
        tr_filt = zdat[0].copy()
        tr_filt.filter('bandpass', freqmin=bpz1[0], freqmax=bpz1[1],
                       zerophase=False)
        tr_filt.taper(max_percentage=0.05, type='hann')
        z_copy[0].data = tr_filt.data
        ##############################################################
        # check length of waveform and compare with cut times

        # for global seismology: use tau-p method for estimating travel times (needs source and station coords.)
        # if not given: sets Lc to infinity to use full stream
        if use_taup == True:
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
                    model = TauPyModel(taup_model)
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
                    arrP, estFirstP = min([(arr, arr.time) for arr in phases['P']], key = lambda t: t[1])
                    arrS, estFirstS = min([(arr, arr.time) for arr in phases['S']], key = lambda t: t[1])
                    print('autopick: estimated first arrivals for P: {} s, S:{} s after event'
                          ' origin time using TauPy'.format(estFirstP, estFirstS))

                    # modifiy pstart and pstop relative to estimated first P arrival (relative to station time axis)
                    pstart += (source_origin.time + estFirstP) - zdat[0].stats.starttime
                    pstop += (source_origin.time + estFirstP) - zdat[0].stats.starttime
                    print('autopick: CF calculation times respectively:'
                          ' pstart: {} s, pstop: {} s'.format(pstart, pstop))
                elif not origin:
                    print('No source origins given!')

        # make sure pstart and pstop are inside zdat[0]
        pstart = max(pstart, 0)
        pstop = min(pstop, len(zdat[0])*zdat[0].stats.delta)

        if not use_taup == True or origin:
            Lc = pstop - pstart

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
        if algoP == 'HOS':
            # calculate HOS-CF using subclass HOScf of class
            # CharacteristicFunction
            cf1 = HOScf(z_copy, cuttimes, thosmw, hosorder)  # instance of HOScf
        elif algoP == 'ARZ':
            # calculate ARZ-CF using subclass ARZcf of class
            # CharcteristicFunction
            cf1 = ARZcf(z_copy, cuttimes, tpred1z, Parorder, tdet1z,
                        addnoise)  # instance of ARZcf
        ##############################################################
        # calculate AIC-HOS-CF using subclass AICcf of class
        # CharacteristicFunction
        # class needs stream object => build it
        assert isinstance(cf1, CharacteristicFunction), 'cf2 is not set ' \
                                                        'correctly: maybe the algorithm name ({algoP}) is ' \
                                                        'corrupted'.format(
            algoP=algoP)
        tr_aic = tr_filt.copy()
        tr_aic.data = cf1.getCF()
        z_copy[0].data = tr_aic.data
        aiccf = AICcf(z_copy, cuttimes)  # instance of AICcf
        ##############################################################
        # get prelimenary onset time from AIC-HOS-CF using subclass AICPicker
        # of class AutoPicking
        key = 'aicFig'
        if fig_dict:
            fig = fig_dict[key]
            linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
        else:
            fig = None
            linecolor = 'k'
        aicpick = AICPicker(aiccf, tsnrz, pickwinP, iplot, None, aictsmoothP, fig=fig, linecolor=linecolor)
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
                msg = 'One or more horizontal component(s) missing!\nSignal ' \
                      'length only checked on vertical component!\n' \
                      'Decreasing minsiglengh from {0} to ' \
                      '{1}'.format(minsiglength, minsiglength / 2)
                if verbose: print(msg)
                key = 'slength'
                if fig_dict:
                    fig = fig_dict[key]
                    linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                else:
                    fig = None
                    linecolor = 'k'
                Pflag = checksignallength(zne, aicpick.getpick(), tsnrz,
                                          minsiglength / 2,
                                          nfacsl, minpercent, iplot,
                                          fig, linecolor)
            else:
                # filter and taper horizontal traces
                trH1_filt = edat.copy()
                trH2_filt = ndat.copy()
                trH1_filt.filter('bandpass', freqmin=bph1[0],
                                 freqmax=bph1[1],
                                 zerophase=False)
                trH2_filt.filter('bandpass', freqmin=bph1[0],
                                 freqmax=bph1[1],
                                 zerophase=False)
                trH1_filt.taper(max_percentage=0.05, type='hann')
                trH2_filt.taper(max_percentage=0.05, type='hann')
                zne += trH1_filt
                zne += trH2_filt
                if fig_dict:
                    fig = fig_dict['slength']
                    linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                else:
                    fig = None
                    linecolor = 'k'
                Pflag = checksignallength(zne, aicpick.getpick(), tsnrz,
                                          minsiglength,
                                          nfacsl, minpercent, iplot,
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
                    Pflag = checkZ4S(zne, aicpick.getpick(), zfac,
                                     tsnrz[2], iplot, fig, linecolor)
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
        if (slope >= minAICPslope and
            aicpick.getSNR() >= minAICPSNR and Pflag == 1):
            aicPflag = 1
            msg = 'AIC P-pick passes quality control: Slope: {0} counts/s, ' \
                  'SNR: {1}\nGo on with refined picking ...\n' \
                  'autopickstation: re-filtering vertical trace ' \
                  '...'.format(aicpick.getSlope(), aicpick.getSNR())
            if verbose: print(msg)
            # re-filter waveform with larger bandpass
            z_copy = zdat.copy()
            tr_filt = zdat[0].copy()
            tr_filt.filter('bandpass', freqmin=bpz2[0], freqmax=bpz2[1],
                           zerophase=False)
            tr_filt.taper(max_percentage=0.05, type='hann')
            z_copy[0].data = tr_filt.data
            #############################################################
            # re-calculate CF from re-filtered trace in vicinity of initial
            # onset
            cuttimes2 = [round(max([aicpick.getpick() - Precalcwin, 0])),
                         round(min([len(zdat[0].data) * zdat[0].stats.delta,
                                    aicpick.getpick() + Precalcwin]))]
            cf2 = None
            if algoP == 'HOS':
                # calculate HOS-CF using subclass HOScf of class
                # CharacteristicFunction
                cf2 = HOScf(z_copy, cuttimes2, thosmw,
                            hosorder)  # instance of HOScf
            elif algoP == 'ARZ':
                # calculate ARZ-CF using subclass ARZcf of class
                # CharcteristicFunction
                cf2 = ARZcf(z_copy, cuttimes2, tpred1z, Parorder, tdet1z,
                            addnoise)  # instance of ARZcf
            ##############################################################
            # get refined onset time from CF2 using class Picker
            assert isinstance(cf2, CharacteristicFunction), 'cf2 is not set ' \
                                                            'correctly: maybe the algorithm name ({algoP}) is ' \
                                                            'corrupted'.format(
                algoP=algoP)
            if fig_dict:
                fig = fig_dict['refPpick']
                linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
            else:
                fig = None
                linecolor = 'k'
            refPpick = PragPicker(cf2, tsnrz, pickwinP, iplot, ausP, tsmoothP,
                                  aicpick.getpick(), fig, linecolor)
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
                    epickP, lpickP, Perror = earllatepicker(z_copy, nfacP, tsnrz,
                                                            mpickP, iplot, fig=fig,
                                                            linecolor=linecolor)
                else:
                    epickP, lpickP, Perror = earllatepicker(z_copy, nfacP, tsnrz,
                                                            mpickP, iplot)

                # get SNR
                [SNRP, SNRPdB, Pnoiselevel] = getSNR(z_copy, tsnrz, mpickP)

                # weight P-onset using symmetric error
                if Perror == None:
                    Pweight = 4
                else:
                    if Perror <= timeerrorsP[0]:
                        Pweight = 0
                    elif timeerrorsP[0] < Perror <= timeerrorsP[1]:
                        Pweight = 1
                    elif timeerrorsP[1] < Perror <= timeerrorsP[2]:
                        Pweight = 2
                    elif timeerrorsP[2] < Perror <= timeerrorsP[3]:
                        Pweight = 3
                    elif Perror > timeerrorsP[3]:
                        Pweight = 4

                ##############################################################
                # get first motion of P onset
                # certain quality required
                if Pweight <= minfmweight and SNRP >= minFMSNR:
                    if iplot:
                        if fig_dict:
                            fig = fig_dict['fm_picker']
                            linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                        else:
                            fig = None
                        FM = fmpicker(zdat, z_copy, fmpickwin, mpickP, iplot, fig, linecolor)
                    else:
                        FM = fmpicker(zdat, z_copy, fmpickwin, mpickP, iplot)
                else:
                    FM = 'N'

                msg = "autopickstation: P-weight: {0}, " \
                      "SNR: {1}, SNR[dB]: {2}, Polarity: {3}".format(Pweight,
                                                                     SNRP,
                                                                     SNRPdB,
                                                                     FM)
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
                                                       minAICPSNR,
                                                       minAICPslope)
            if verbose: print(msg)
            Sflag = 0

    else:
        print('autopickstation: No vertical component data available!, '
              'Skipping station!')

    if ((len(edat) > 0 and len(ndat) == 0) or (
        len(ndat) > 0 and len(edat) == 0)) and Pweight < 4:
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
            round(max([mpickP + sstart, 0])), # MP MP relative time axis
            round(min([
                mpickP + sstop,
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


        if algoS == 'ARH':
            # re-create stream object including both horizontal components
            hdat = edat.copy()
            hdat += ndat
            h_copy = hdat.copy()
            # filter and taper data
            trH1_filt = hdat[0].copy()
            trH2_filt = hdat[1].copy()
            trH1_filt.filter('bandpass', freqmin=bph1[0], freqmax=bph1[1],
                             zerophase=False)
            trH2_filt.filter('bandpass', freqmin=bph1[0], freqmax=bph1[1],
                             zerophase=False)
            trH1_filt.taper(max_percentage=0.05, type='hann')
            trH2_filt.taper(max_percentage=0.05, type='hann')
            h_copy[0].data = trH1_filt.data
            h_copy[1].data = trH2_filt.data
        elif algoS == 'AR3':
            # re-create stream object including all components
            hdat = zdat.copy()
            hdat += edat
            hdat += ndat
            h_copy = hdat.copy()
            # filter and taper data
            trH1_filt = hdat[0].copy()
            trH2_filt = hdat[1].copy()
            trH3_filt = hdat[2].copy()
            trH1_filt.filter('bandpass', freqmin=bph1[0], freqmax=bph1[1],
                             zerophase=False)
            trH2_filt.filter('bandpass', freqmin=bph1[0], freqmax=bph1[1],
                             zerophase=False)
            trH3_filt.filter('bandpass', freqmin=bph1[0], freqmax=bph1[1],
                             zerophase=False)
            trH1_filt.taper(max_percentage=0.05, type='hann')
            trH2_filt.taper(max_percentage=0.05, type='hann')
            trH3_filt.taper(max_percentage=0.05, type='hann')
            h_copy[0].data = trH1_filt.data
            h_copy[1].data = trH2_filt.data
            h_copy[2].data = trH3_filt.data
        ##############################################################
        if algoS == 'ARH':
            # calculate ARH-CF using subclass ARHcf of class
            # CharcteristicFunction
            arhcf1 = ARHcf(h_copy, cuttimesh, tpred1h, Sarorder, tdet1h,
                           addnoise)  # instance of ARHcf
        elif algoS == 'AR3':
            # calculate ARH-CF using subclass AR3cf of class
            # CharcteristicFunction
            arhcf1 = AR3Ccf(h_copy, cuttimesh, tpred1h, Sarorder, tdet1h,
                            addnoise)  # instance of ARHcf
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
        aicarhpick = AICPicker(haiccf, tsnrh, pickwinS, iplot, None,
                               aictsmoothS, fig=fig, linecolor=linecolor)
        ###############################################################
        # go on with processing if AIC onset passes quality control
        slope = aicarhpick.getSlope()
        if not slope:
            slope = 0
        if (slope >= minAICSslope and
            aicarhpick.getSNR() >= minAICSSNR and
            aicarhpick.getpick() is not None):
            aicSflag = 1
            msg = 'AIC S-pick passes quality control: Slope: {0} counts/s, ' \
                  'SNR: {1}\nGo on with refined picking ...\n' \
                  'autopickstation: re-filtering horizontal traces ' \
                  '...'.format(aicarhpick.getSlope(), aicarhpick.getSNR())
            if verbose: print(msg)
            # re-calculate CF from re-filtered trace in vicinity of initial
            # onset
            cuttimesh2 = [round(aicarhpick.getpick() - Srecalcwin),
                          round(aicarhpick.getpick() + Srecalcwin)]
            # re-filter waveform with larger bandpass
            h_copy = hdat.copy()
            # filter and taper data
            if algoS == 'ARH':
                trH1_filt = hdat[0].copy()
                trH2_filt = hdat[1].copy()
                trH1_filt.filter('bandpass', freqmin=bph2[0], freqmax=bph2[1],
                                 zerophase=False)
                trH2_filt.filter('bandpass', freqmin=bph2[0], freqmax=bph2[1],
                                 zerophase=False)
                trH1_filt.taper(max_percentage=0.05, type='hann')
                trH2_filt.taper(max_percentage=0.05, type='hann')
                h_copy[0].data = trH1_filt.data
                h_copy[1].data = trH2_filt.data
                #############################################################
                arhcf2 = ARHcf(h_copy, cuttimesh2, tpred2h, Sarorder, tdet2h,
                               addnoise)  # instance of ARHcf
            elif algoS == 'AR3':
                trH1_filt = hdat[0].copy()
                trH2_filt = hdat[1].copy()
                trH3_filt = hdat[2].copy()
                trH1_filt.filter('bandpass', freqmin=bph2[0], freqmax=bph2[1],
                                 zerophase=False)
                trH2_filt.filter('bandpass', freqmin=bph2[0], freqmax=bph2[1],
                                 zerophase=False)
                trH3_filt.filter('bandpass', freqmin=bph2[0], freqmax=bph2[1],
                                 zerophase=False)
                trH1_filt.taper(max_percentage=0.05, type='hann')
                trH2_filt.taper(max_percentage=0.05, type='hann')
                trH3_filt.taper(max_percentage=0.05, type='hann')
                h_copy[0].data = trH1_filt.data
                h_copy[1].data = trH2_filt.data
                h_copy[2].data = trH3_filt.data
                #############################################################
                arhcf2 = AR3Ccf(h_copy, cuttimesh2, tpred2h, Sarorder, tdet2h,
                                addnoise)  # instance of ARHcf

            # get refined onset time from CF2 using class Picker
            if fig_dict:
                fig = fig_dict['refSpick']
                linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
            else:
                fig = None
                linecolor = 'k'
            refSpick = PragPicker(arhcf2, tsnrh, pickwinS, iplot, ausS,
                                  tsmoothS, aicarhpick.getpick(), fig, linecolor)
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
                    epickS1, lpickS1, Serror1 = earllatepicker(h_copy, nfacS,
                                                               tsnrh,
                                                               mpickS, iplot,
                                                               fig=fig,
                                                               linecolor=linecolor)
                else:
                    epickS1, lpickS1, Serror1 = earllatepicker(h_copy, nfacS,
                                                               tsnrh,
                                                               mpickS, iplot)

                h_copy[0].data = trH2_filt.data
                if iplot:
                    if fig_dict:
                        fig = fig_dict['el_S2pick']
                        linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
                    else:
                        fig = None
                        linecolor = ''
                    epickS2, lpickS2, Serror2 = earllatepicker(h_copy, nfacS,
                                                               tsnrh,
                                                               mpickS, iplot,
                                                               fig=fig,
                                                               linecolor=linecolor)
                else:
                    epickS2, lpickS2, Serror2 = earllatepicker(h_copy, nfacS,
                                                               tsnrh,
                                                               mpickS, iplot)
                if epickS1 is not None and epickS2 is not None:
                    if algoS == 'ARH':
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
                    elif algoS == 'AR3':
                        [epickS3, lpickS3, Serror3] = earllatepicker(h_copy,
                                                                     nfacS,
                                                                     tsnrh,
                                                                     mpickS,
                                                                     iplot)
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
                    [SNRS, SNRSdB, Snoiselevel] = getSNR(h_copy, tsnrh, mpickS)

                    # weight S-onset using symmetric error
                    if Serror <= timeerrorsS[0]:
                        Sweight = 0
                    elif timeerrorsS[0] < Serror <= timeerrorsS[1]:
                        Sweight = 1
                    elif Perror > timeerrorsS[1] and Serror <= timeerrorsS[2]:
                        Sweight = 2
                    elif timeerrorsS[2] < Serror <= timeerrorsS[3]:
                        Sweight = 3
                    elif Serror > timeerrorsS[3]:
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
                  ''.format(aicarhpick.getSNR(),
                            aicarhpick.getSlope(),
                            minAICSSNR,
                            minAICSslope)
            if verbose: print(msg)

            ############################################################
            # get Wood-Anderson peak-to-peak amplitude
            # initialize Data object
            # re-create stream object including both horizontal components
            hdat = edat.copy()
            hdat += ndat
 
    else:
        print('autopickstation: No horizontal component data available or ' \
              'bad P onset, skipping S picking!')

    ##############################################################
    try:
       iplot = int(iplot)
    except:
       if iplot == True or iplot == 'True':
           iplot = 2
       else:
           iplot = 0

    if iplot > 0:
        # plot vertical trace
        if fig_dict == None or fig_dict == 'None':
            fig = plt.figure()
            plt_flag = 1
            linecolor = 'k'
        else:
            fig = fig_dict['mainFig']
            linecolor = fig_dict['plot_style']['linecolor']['rgba_mpl']
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
                try: input()
                except SyntaxError: pass
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
        lpickP = zdat[0].stats.starttime + timeerrorsP[3]
        epickP = zdat[0].stats.starttime - timeerrorsP[3]
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
        lpickS = hdat.stats.starttime + timeerrorsS[3]
        epickS = hdat.stats.starttime - timeerrorsS[3]
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
    '''
    Repicking of bad onsets. Uses theoretical onset times from NLLoc-location file.

    :param wf: waveform, obspy stream object

    :param NLLocfile: path/name of NLLoc-location file

    :param picks: dictionary of available onset times

    :param badpicks: picks to be repicked

    :param pickparameter: picking parameters from autoPyLoT-input file
    '''

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
            pickparameter.setParam(pstop=pickparameter.get('pstart') + \
                                         (pickparameter.get('Precalcwin')))
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
