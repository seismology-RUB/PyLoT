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
from pylot.core.io.inputs import AutoPickParameter
from pylot.core.pick.picker import AICPicker, PragPicker
from pylot.core.pick.charfuns import CharacteristicFunction
from pylot.core.pick.charfuns import HOScf, AICcf, ARZcf, ARHcf, AR3Ccf
from pylot.core.pick.utils import checksignallength, checkZ4S, earllatepicker, \
    getSNR, fmpicker, checkPonsets, wadaticheck
from pylot.core.util.utils import getPatternLine
from pylot.core.io.data import Data


def autopickevent(data, param):
    stations = []
    all_onsets = {}

    # get some parameters for quality control from
    # parameter input file (usually autoPyLoT.in).
    wdttolerance = param.get('wdttolerance')
    mdttolerance = param.get('mdttolerance')
    iplot = param.get('iplot')
    apverbose = param.get('apverbose')
    for n in range(len(data)):
        station = data[n].stats.station
        if station not in stations:
            stations.append(station)
        else:
            continue

    for station in stations:
        topick = data.select(station=station)
        all_onsets[station] = autopickstation(topick, param, verbose=apverbose)

    # quality control
    # median check and jackknife on P-onset times
    jk_checked_onsets = checkPonsets(all_onsets, mdttolerance, iplot)
    # check S-P times (Wadati)
    return wadaticheck(jk_checked_onsets, wdttolerance, iplot)


def autopickstation(wfstream, pickparam, verbose=False):
    """
    :param wfstream: `~obspy.core.stream.Stream`  containing waveform
    :type wfstream: obspy.core.stream.Stream

    :param pickparam: container of picking parameters from input file,
           usually autoPyLoT.in
    :type pickparam: AutoPickParameter
    :param verbose:
    :type verbose: bool

    """

    # declaring pickparam variables (only for convenience)
    # read your autoPyLoT.in for details!

    # special parameters for P picking
    algoP = pickparam.get('algoP')
    iplot = pickparam.get('iplot')
    pstart = pickparam.get('pstart')
    pstop = pickparam.get('pstop')
    thosmw = pickparam.get('tlta')
    tsnrz = pickparam.get('tsnrz')
    hosorder = pickparam.get('hosorder')
    bpz1 = pickparam.get('bpz1')
    bpz2 = pickparam.get('bpz2')
    pickwinP = pickparam.get('pickwinP')
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
    picker = 'autoPyLoT'  # name of the picking programm

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

    if algoP == 'HOS' or algoP == 'ARZ' and zdat is not None:
        msg = '##########################################\nautopickstation:' \
              ' Working on P onset of station {station}\nFiltering vertical ' \
              'trace ...\n{data}'.format(station=zdat[0].stats.station,
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
        Lc = pstop - pstart
        Lwf = zdat[0].stats.endtime - zdat[0].stats.starttime
        Ldiff = Lwf - Lc
        if Ldiff < 0:
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
        aicpick = AICPicker(aiccf, tsnrz, pickwinP, iplot, None, tsmoothP)
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
                Pflag = checksignallength(zne, aicpick.getpick(), tsnrz,
                                          minsiglength / 2,
                                          nfacsl, minpercent, iplot)
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
                Pflag = checksignallength(zne, aicpick.getpick(), tsnrz,
                                          minsiglength,
                                          nfacsl, minpercent, iplot)

            if Pflag == 1:
                # check for spuriously picked S onset
                # both horizontal traces needed
                if len(ndat) == 0 or len(edat) == 0:
                    msg = 'One or more horizontal components missing!\n' \
                          'Skipping control function checkZ4S.'
                    if verbose: print(msg)
                else:
                    Pflag = checkZ4S(zne, aicpick.getpick(), zfac,
                                     tsnrz[3], iplot)
                    if Pflag == 0:
                        Pmarker = 'SinsteadP'
                        Pweight = 9
            else:
                Pmarker = 'shortsignallength'
                Pweight = 9
        ##############################################################
        # go on with processing if AIC onset passes quality control
        if (aicpick.getSlope() >= minAICPslope and
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
            refPpick = PragPicker(cf2, tsnrz, pickwinP, iplot, ausP, tsmoothP,
                                  aicpick.getpick())
            mpickP = refPpick.getpick()
            #############################################################
            if mpickP is not None:
                # quality assessment
                # get earliest/latest possible pick and symmetrized uncertainty
                [epickP, lpickP, Perror] = earllatepicker(z_copy, nfacP, tsnrz,
                                                          mpickP, iplot)

                # get SNR
                [SNRP, SNRPdB, Pnoiselevel] = getSNR(z_copy, tsnrz, mpickP)

                # weight P-onset using symmetric error
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
                    FM = fmpicker(zdat, z_copy, fmpickwin, mpickP, iplot)
                else:
                    FM = 'N'

                msg = "autopickstation: P-weight: {0}, " \
                      "SNR: {1}, SNR[dB]: {2}, Polarity: {3}".format(Pweight,
                                                                     SNRP,
                                                                     SNRPdB,
                                                                     FM)
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

    if edat is not None and ndat is not None and len(edat) > 0 and len(
            ndat) > 0 and Pweight < 4:
        msg = 'Go on picking S onset ...\n' \
              '##################################################\n' \
              'Working on S onset of station {0}\nFiltering horizontal ' \
              'traces ...'.format(edat[0].stats.station)
        if verbose: print(msg)
        # determine time window for calculating CF after P onset
        cuttimesh = [round(max([mpickP + sstart, 0])),
                     round(min([mpickP + sstop, Lwf]))]

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
        aicarhpick = AICPicker(haiccf, tsnrh, pickwinS, iplot, None,
                               aictsmoothS)
        ###############################################################
        # go on with processing if AIC onset passes quality control
        if (aicarhpick.getSlope() >= minAICSslope and
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
            refSpick = PragPicker(arhcf2, tsnrh, pickwinS, iplot, ausS,
                                  tsmoothS, aicarhpick.getpick())
            mpickS = refSpick.getpick()
            #############################################################
            if mpickS is not None:
                # quality assessment
                # get earliest/latest possible pick and symmetrized uncertainty
                h_copy[0].data = trH1_filt.data
                [epickS1, lpickS1, Serror1] = earllatepicker(h_copy, nfacS,
                                                             tsnrh,
                                                             mpickS, iplot)

                h_copy[0].data = trH2_filt.data
                [epickS2, lpickS2, Serror2] = earllatepicker(h_copy, nfacS,
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
                          '################################################'
                          ''.format(Sweight, SNRS, SNRSdB))
                ################################################################
                # get Wood-Anderson peak-to-peak amplitude
                # initialize Data object
                data = Data()
                # re-create stream object including both horizontal components
                hdat = edat.copy()
                hdat += ndat
        else:
            msg = 'Bad initial (AIC) S-pick, skipping this onset!\n' \
                  'AIC-SNR={0}, AIC-Slope={1}counts/s\n' \
                  '(min. AIC-SNR={2}, ' \
                  'min. AIC-Slope={3}counts/s)\n' \
                  '################################################' \
                  ''.format(aicarhpick.getSNR(),
                            aicarhpick.getSlope(),
                            minAICSSNR,
                            minAICSslope)
            if verbose: print(msg)

            ############################################################
            # get Wood-Anderson peak-to-peak amplitude
            # initialize Data object
            data = Data()
            # re-create stream object including both horizontal components
            hdat = edat.copy()
            hdat += ndat
    else:
        print('autopickstation: No horizontal component data available or ' \
              'bad P onset, skipping S picking!')

    ##############################################################
    if iplot > 0:
        # plot vertical trace
        plt.figure()
        plt.subplot(3, 1, 1)
        tdata = np.arange(0, zdat[0].stats.npts / tr_filt.stats.sampling_rate,
                          tr_filt.stats.delta)
        # check equal length of arrays, sometimes they are different!?
        wfldiff = len(tr_filt.data) - len(tdata)
        if wfldiff < 0:
            tdata = tdata[0:len(tdata) - abs(wfldiff)]
        p1, = plt.plot(tdata, tr_filt.data / max(tr_filt.data), 'k')
        if Pweight < 4:
            p2, = plt.plot(cf1.getTimeArray(), cf1.getCF() / max(cf1.getCF()),
                           'b')
            if aicPflag == 1:
                p3, = plt.plot(cf2.getTimeArray(),
                               cf2.getCF() / max(cf2.getCF()), 'm')
                p4, = plt.plot([aicpick.getpick(), aicpick.getpick()], [-1, 1],
                               'r')
                plt.plot([aicpick.getpick() - 0.5, aicpick.getpick() + 0.5],
                         [1, 1], 'r')
                plt.plot([aicpick.getpick() - 0.5, aicpick.getpick() + 0.5],
                         [-1, -1], 'r')
                p5, = plt.plot([refPpick.getpick(), refPpick.getpick()],
                               [-1.3, 1.3], 'r', linewidth=2)
                plt.plot([refPpick.getpick() - 0.5, refPpick.getpick() + 0.5],
                         [1.3, 1.3], 'r', linewidth=2)
                plt.plot([refPpick.getpick() - 0.5, refPpick.getpick() + 0.5],
                         [-1.3, -1.3], 'r', linewidth=2)
                plt.plot([lpickP, lpickP], [-1.1, 1.1], 'r--')
                plt.plot([epickP, epickP], [-1.1, 1.1], 'r--')
                plt.legend([p1, p2, p3, p4, p5],
                           ['Data', 'CF1', 'CF2', 'Initial P Onset',
                            'Final P Pick'])
                plt.title('%s, %s, P Weight=%d, SNR=%7.2f, SNR[dB]=%7.2f '
                          'Polarity: %s' % (tr_filt.stats.station,
                                            tr_filt.stats.channel,
                                            Pweight,
                                            SNRP,
                                            SNRPdB,
                                            FM))
            else:
                plt.legend([p1, p2], ['Data', 'CF1'])
                plt.title('%s, P Weight=%d, SNR=None, '
                          'SNRdB=None' % (tr_filt.stats.channel, Pweight))
        else:
            plt.title('%s, %s, P Weight=%d' % (tr_filt.stats.station,
                                               tr_filt.stats.channel,
                                               Pweight))

        plt.yticks([])
        plt.ylim([-1.5, 1.5])
        plt.ylabel('Normalized Counts')
        plt.suptitle(tr_filt.stats.starttime)

        if len(edat[0]) > 1 and len(ndat[0]) > 1 and Sflag == 1:
            # plot horizontal traces
            plt.subplot(3, 1, 2)
            th1data = np.arange(0,
                                trH1_filt.stats.npts /
                                trH1_filt.stats.sampling_rate,
                                trH1_filt.stats.delta)
            # check equal length of arrays, sometimes they are different!?
            wfldiff = len(trH1_filt.data) - len(th1data)
            if wfldiff < 0:
                th1data = th1data[0:len(th1data) - abs(wfldiff)]
            p21, = plt.plot(th1data, trH1_filt.data / max(trH1_filt.data), 'k')
            if Pweight < 4:
                p22, = plt.plot(arhcf1.getTimeArray(),
                                arhcf1.getCF() / max(arhcf1.getCF()), 'b')
                if aicSflag == 1:
                    p23, = plt.plot(arhcf2.getTimeArray(),
                                    arhcf2.getCF() / max(arhcf2.getCF()), 'm')
                    p24, = plt.plot(
                        [aicarhpick.getpick(), aicarhpick.getpick()],
                        [-1, 1], 'g')
                    plt.plot(
                        [aicarhpick.getpick() - 0.5,
                         aicarhpick.getpick() + 0.5],
                        [1, 1], 'g')
                    plt.plot(
                        [aicarhpick.getpick() - 0.5,
                         aicarhpick.getpick() + 0.5],
                        [-1, -1], 'g')
                    p25, = plt.plot([refSpick.getpick(), refSpick.getpick()],
                                    [-1.3, 1.3], 'g', linewidth=2)
                    plt.plot(
                        [refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                        [1.3, 1.3], 'g', linewidth=2)
                    plt.plot(
                        [refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                        [-1.3, -1.3], 'g', linewidth=2)
                    plt.plot([lpickS, lpickS], [-1.1, 1.1], 'g--')
                    plt.plot([epickS, epickS], [-1.1, 1.1], 'g--')
                    plt.legend([p21, p22, p23, p24, p25],
                               ['Data', 'CF1', 'CF2', 'Initial S Onset',
                                'Final S Pick'])
                    plt.title('%s, S Weight=%d, SNR=%7.2f, SNR[dB]=%7.2f' % (
                        trH1_filt.stats.channel,
                        Sweight, SNRS, SNRSdB))
                else:
                    plt.legend([p21, p22], ['Data', 'CF1'])
                    plt.title('%s, S Weight=%d, SNR=None, SNRdB=None' % (
                        trH1_filt.stats.channel, Sweight))
            plt.yticks([])
            plt.ylim([-1.5, 1.5])
            plt.ylabel('Normalized Counts')
            plt.suptitle(trH1_filt.stats.starttime)

            plt.subplot(3, 1, 3)
            th2data = np.arange(0,
                                trH2_filt.stats.npts /
                                trH2_filt.stats.sampling_rate,
                                trH2_filt.stats.delta)
            # check equal length of arrays, sometimes they are different!?
            wfldiff = len(trH2_filt.data) - len(th2data)
            if wfldiff < 0:
                th2data = th2data[0:len(th2data) - abs(wfldiff)]
            plt.plot(th2data, trH2_filt.data / max(trH2_filt.data), 'k')
            if Pweight < 4:
                p22, = plt.plot(arhcf1.getTimeArray(),
                                arhcf1.getCF() / max(arhcf1.getCF()), 'b')
                if aicSflag == 1:
                    p23, = plt.plot(arhcf2.getTimeArray(),
                                    arhcf2.getCF() / max(arhcf2.getCF()), 'm')
                    p24, = plt.plot(
                        [aicarhpick.getpick(), aicarhpick.getpick()],
                        [-1, 1], 'g')
                    plt.plot(
                        [aicarhpick.getpick() - 0.5,
                         aicarhpick.getpick() + 0.5],
                        [1, 1], 'g')
                    plt.plot(
                        [aicarhpick.getpick() - 0.5,
                         aicarhpick.getpick() + 0.5],
                        [-1, -1], 'g')
                    p25, = plt.plot([refSpick.getpick(), refSpick.getpick()],
                                    [-1.3, 1.3], 'g', linewidth=2)
                    plt.plot(
                        [refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                        [1.3, 1.3], 'g', linewidth=2)
                    plt.plot(
                        [refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                        [-1.3, -1.3], 'g', linewidth=2)
                    plt.plot([lpickS, lpickS], [-1.1, 1.1], 'g--')
                    plt.plot([epickS, epickS], [-1.1, 1.1], 'g--')
                    plt.legend([p21, p22, p23, p24, p25],
                               ['Data', 'CF1', 'CF2', 'Initial S Onset',
                                'Final S Pick'])
                else:
                    plt.legend([p21, p22], ['Data', 'CF1'])
            plt.yticks([])
            plt.ylim([-1.5, 1.5])
            plt.xlabel('Time [s] after %s' % tr_filt.stats.starttime)
            plt.ylabel('Normalized Counts')
            plt.title(trH2_filt.stats.channel)
            plt.show()
            raw_input()
            plt.close()
    ##########################################################################
    # calculate "real" onset times
    if lpickP is not None and lpickP == mpickP:
        lpickP += timeerrorsP[0]
    if epickP is not None and epickP == mpickP:
        epickP -= timeerrorsP[0]
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

    if lpickS is not None and lpickS == mpickS:
        lpickS += timeerrorsS[0]
    if epickS is not None and epickS == mpickS:
        epickS -= timeerrorsS[0]
    if mpickS is not None and epickS is not None and lpickS is not None:
        lpickS = edat[0].stats.starttime + lpickS
        epickS = edat[0].stats.starttime + epickS
        mpickS = edat[0].stats.starttime + mpickS
    else:
        # dummy values (start of seismic trace) in order to derive
        # theoretical onset times for iteratve picking
        lpickS = edat[0].stats.starttime + timeerrorsS[3]
        epickS = edat[0].stats.starttime - timeerrorsS[3]
        mpickS = edat[0].stats.starttime

    # create dictionary
    # for P phase
    ppick = dict(lpp=lpickP, epp=epickP, mpp=mpickP, spe=Perror, snr=SNRP,
                 snrdb=SNRPdB, weight=Pweight, fm=FM, w0=None, fc=None, Mo=None,
                 Mw=None, picker=picker, marked=Pmarker)
    # add S phase
    spick = dict(lpp=lpickS, epp=epickS, mpp=mpickS, spe=Serror, snr=SNRS,
                 snrdb=SNRSdB, weight=Sweight, fm=None, picker=picker, Ao=Ao)
    # merge picks into returning dictionary
    picks = dict(P=ppick, S=spick)
    return picks


def iteratepicker(wf, NLLocfile, picks, badpicks, pickparameter):
    '''
    Repicking of bad onsets. Uses theoretical onset times from NLLoc-location file.

    :param wf: waveform, obspy stream object

    :param NLLocfile: path/name of NLLoc-location file

    :param picks: dictionary of available onset times

    :param badpicks: picks to be repicked

    :param pickparameter: picking parameters from autoPyLoT-input file
    '''

    msg = '#######################################################\n' \
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
        msg = '#######################################################\n' \
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
        pickparameter.setParam(
            pstart=max([0, badpicks[i][1] - wf2pick[0].stats.starttime \
                        - pickparameter.get('tlta')]))
        pickparameter.setParam(pstop=pickparameter.get('pstart') + \
                                     (3 * pickparameter.get('tlta')))
        pickparameter.setParam(sstop=pickparameter.get('sstop') / 2)
        pickparameter.setParam(pickwinP=pickparameter.get('pickwinP') / 2)
        pickparameter.setParam(
            Precalcwin=pickparameter.get('Precalcwin') / 2)
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
        newpicks = autopickstation(wf2pick, pickparameter)

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
