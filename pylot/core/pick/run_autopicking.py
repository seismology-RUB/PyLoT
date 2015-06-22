#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Function to run automated picking algorithms using AIC,
HOS and AR prediction. Uses object CharFuns and Picker and
function conglomerate utils.

:author: MAGS2 EP3 working group / Ludger Kueperkoch
"""

import matplotlib.pyplot as plt
import numpy as np
from pylot.core.pick.Picker import *
from pylot.core.pick.CharFuns import *
import pdb
def run_autopicking(wfstream, pickparam):
    """
    :param: wfstream
    :type: `~obspy.core.stream.Stream`

    :param: pickparam
    :type: container of picking parameters from input file,
           usually autoPyLoT.in
    """

    # declaring pickparam variables (only for convenience)
    # read your autoPyLoT.in for details!

    # special parameters for P picking
    algoP = pickparam.getParam('algoP')
    iplot = pickparam.getParam('iplot')
    pstart = pickparam.getParam('pstart')
    pstop = pickparam.getParam('pstop')
    thosmw = pickparam.getParam('tlta')
    tsnrz = pickparam.getParam('tsnrz')
    hosorder = pickparam.getParam('hosorder')
    bpz1 = pickparam.getParam('bpz1')
    bpz2 = pickparam.getParam('bpz2')
    pickwinP = pickparam.getParam('pickwinP')
    tsmoothP = pickparam.getParam('tsmoothP')
    ausP = pickparam.getParam('ausP')
    nfacP = pickparam.getParam('nfacP')
    tpred1z = pickparam.getParam('tpred1z')
    tdet1z = pickparam.getParam('tdet1z')
    Parorder = pickparam.getParam('Parorder')
    addnoise = pickparam.getParam('addnoise')
    Precalcwin = pickparam.getParam('Precalcwin')
    minAICPslope = pickparam.getParam('minAICPslope')
    minAICPSNR = pickparam.getParam('minAICPSNR')
    timeerrorsP = pickparam.getParam('timeerrorsP')
    # special parameters for S picking
    algoS = pickparam.getParam('algoS')
    sstart = pickparam.getParam('sstart')
    sstop = pickparam.getParam('sstop')
    bph1 = pickparam.getParam('bph1')
    bph2 = pickparam.getParam('bph2')
    tsnrh = pickparam.getParam('tsnrh')
    pickwinS = pickparam.getParam('pickwinS')
    tpred1h = pickparam.getParam('tpred1h')
    tdet1h = pickparam.getParam('tdet1h')
    tpred2h = pickparam.getParam('tpred2h')
    tdet2h = pickparam.getParam('tdet2h')
    Sarorder = pickparam.getParam('Sarorder')
    aictsmoothS = pickparam.getParam('aictsmoothS')
    tsmoothS = pickparam.getParam('tsmoothS')
    ausS = pickparam.getParam('ausS')
    minAICSslope = pickparam.getParam('minAICSslope')
    minAICSSNR = pickparam.getParam('minAICSSNR')
    Srecalcwin = pickparam.getParam('Srecalcwin')
    nfacS = pickparam.getParam('nfacS')
    timeerrorsS = pickparam.getParam('timeerrorsS')
    # parameters for first-motion determination
    minFMSNR = pickparam.getParam('minFMSNR')
    fmpickwin = pickparam.getParam('fmpickwin')
    minfmweight = pickparam.getParam('minfmweight')

    # initialize output 
    Pweight = 4   # weight for P onset 
    Sweight = 4   # weight for S onset
    FM = 'N'      # first motion (polarity)
    SNRP = None   # signal-to-noise ratio of P onset
    SNRPdB = None # signal-to-noise ratio of P onset [dB]
    SNRS = None   # signal-to-noise ratio of S onset
    SNRSdB = None # signal-to-noise ratio of S onset [dB]
    mpickP = None # most likely P onset
    lpickP = None # latest possible P onset
    epickP = None # earliest possible P onset
    mpickS = None # most likely S onset
    lpickS = None # latest possible S onset
    epickS = None # earliest possible S onset
    Perror = None # symmetrized picking error P onset
    Serror = None # symmetrized picking error S onset

    aicSflag = 0
    aicPflag = 0
    Sflag = 0

    # split components
    zdat = wfstream.select(component="Z")
    edat = wfstream.select(component="E")
    if len(edat) == 0:  # check for other components
        edat = wfstream.select(component="2")
    ndat = wfstream.select(component="N")
    if len(ndat) == 0:  # check for other components
        ndat = wfstream.select(component="1")

    if algoP == 'HOS' or algoP == 'ARZ' and zdat is not None:
        print '##########################################'
        print 'run_autopicking: Working on P onset of station %s' % zdat[
            0].stats.station
        print 'Filtering vertical trace ...'
        print zdat
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
            print 'run_autopicking: Cutting times are too large for actual ' \
                  'waveform!'
            print 'Use entire waveform instead!'
            pstart = 0
            pstop = len(zdat[0].data) * zdat[0].stats.delta
        cuttimes = [pstart, pstop]
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
        tr_aic = tr_filt.copy()
        tr_aic.data = cf1.getCF()
        z_copy[0].data = tr_aic.data
        aiccf = AICcf(z_copy, cuttimes)  # instance of AICcf
        ##############################################################
        # get prelimenary onset time from AIC-HOS-CF using subclass AICPicker
        # of class AutoPicking
        aicpick = AICPicker(aiccf, tsnrz, pickwinP, iplot, None, tsmoothP)
        ##############################################################
        # go on with processing if AIC onset passes quality control
        if (aicpick.getSlope() >= minAICPslope and
                aicpick.getSNR() >= minAICPSNR):
            aicPflag = 1
            print 'AIC P-pick passes quality control: Slope: %f, SNR: %f' % \
                  (aicpick.getSlope(), aicpick.getSNR())
            print 'Go on with refined picking ...'
            # re-filter waveform with larger bandpass
            print 'run_autopicking: re-filtering vertical trace ...'
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
            refPpick = PragPicker(cf2, tsnrz, pickwinP, iplot, ausP, tsmoothP,
                                  aicpick.getpick())
            mpickP = refPpick.getpick()
            #############################################################
            if mpickP is not None:
            	# quality assessment
            	# get earliest and latest possible pick and symmetrized uncertainty
            	[lpickP, epickP, Perror] = earllatepicker(z_copy, nfacP, tsnrz, mpickP, iplot)

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

            	print 'run_autopicking: P-weight: %d, SNR: %f, SNR[dB]: %f, ' \
                  	'Polarity: %s' % (Pweight, SNRP, SNRPdB, FM)
            	Sflag = 1

        else:
            print 'Bad initial (AIC) P-pick, skip this onset!'
            print 'AIC-SNR=', aicpick.getSNR(), 'AIC-Slope=', aicpick.getSlope()
            Sflag = 0

    else:
        print 'run_autopicking: No vertical component data available, ' \
              'skipping station!'

    if edat is not None and ndat is not None and len(edat) > 0 and len(
            ndat) > 0 and Pweight < 4:
        print 'Go on picking S onset ...'
        print '##################################################'
        print 'Working on S onset of station %s' % edat[0].stats.station
        print 'Filtering horizontal traces ...'

        # determine time window for calculating CF after P onset
        cuttimesh = [round(max([mpickP + sstart, 0])),
                     round(min([mpickP + sstop, Lwf]))]

        if algoS == 'ARH':
            print edat, ndat
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
            print zdat, edat, ndat
            # re-create stream object including both horizontal components
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
                aicarhpick.getSNR() >= minAICSSNR):
            aicSflag = 1
            print 'AIC S-pick passes quality control: Slope: %f, SNR: %f' \
                  % (aicarhpick.getSlope(), aicarhpick.getSNR())
            print 'Go on with refined picking ...'
            # re-calculate CF from re-filtered trace in vicinity of initial
            # onset
            cuttimesh2 = [round(aicarhpick.getpick() - Srecalcwin),
                          round(aicarhpick.getpick() + Srecalcwin)]
            # re-filter waveform with larger bandpass
            print 'run_autopicking: re-filtering horizontal traces...'
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
            	# get earliest and latest possible pick and symmetrized uncertainty
            	h_copy[0].data = trH1_filt.data
            	[lpickS1, epickS1, Serror1] = earllatepicker(h_copy, nfacS, tsnrh,
               	                                          mpickS, iplot)
            	h_copy[0].data = trH2_filt.data
            	[lpickS2, epickS2, Serror2] = earllatepicker(h_copy, nfacS, tsnrh,
                                                         mpickS, iplot)
            	if algoS == 'ARH':
              		# get earliest pick of both earliest possible picks
                	epick = [epickS1, epickS2]
                	lpick = [lpickS1, lpickS2]
                	pickerr = [Serror1, Serror2]
                	ipick = np.argmin([epickS1, epickS2])
            	elif algoS == 'AR3':
                	[lpickS3, epickS3, Serror3] = earllatepicker(h_copy, nfacS,
                                                             tsnrh,
                                                             mpickS, iplot)
                	# get earliest pick of all three picks
                	epick = [epickS1, epickS2, epickS3]
                	lpick = [lpickS1, lpickS2, lpickS3]
                	pickerr = [Serror1, Serror2, Serror3]
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

            	print 'run_autopicking: S-weight: %d, SNR: %f, SNR[dB]: %f' % (
                	Sweight, SNRS, SNRSdB)

        else:
            print 'Bad initial (AIC) S-pick, skip this onset!'
            print 'AIC-SNR=', aicarhpick.getSNR(), \
                  'AIC-Slope=', aicarhpick.getSlope()

    else:
        print 'run_autopicking: No horizontal component data available or ' \
              'bad P onset, skipping S picking!'

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
        plt.yticks([])
        plt.ylim([-1.5, 1.5])
        plt.ylabel('Normalized Counts')
        plt.suptitle(tr_filt.stats.starttime)

        if len(edat) > 1 and len(ndat) > 1 and Sflag == 1:
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
                    p24, = plt.plot([aicarhpick.getpick(), aicarhpick.getpick()],
                                    [-1, 1], 'g')
                    plt.plot(
                        [aicarhpick.getpick() - 0.5, aicarhpick.getpick() + 0.5],
                        [1, 1], 'g')
                    plt.plot(
                        [aicarhpick.getpick() - 0.5, aicarhpick.getpick() + 0.5],
                        [-1, -1], 'g')
                    p25, = plt.plot([refSpick.getpick(), refSpick.getpick()],
                                    [-1.3, 1.3], 'g', linewidth=2)
                    plt.plot([refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                             [1.3, 1.3], 'g', linewidth=2)
                    plt.plot([refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
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
                    p24, = plt.plot([aicarhpick.getpick(), aicarhpick.getpick()],
                                    [-1, 1], 'g')
                    plt.plot(
                        [aicarhpick.getpick() - 0.5, aicarhpick.getpick() + 0.5],
                        [1, 1], 'g')
                    plt.plot(
                        [aicarhpick.getpick() - 0.5, aicarhpick.getpick() + 0.5],
                        [-1, -1], 'g')
                    p25, = plt.plot([refSpick.getpick(), refSpick.getpick()],
                             [-1.3, 1.3], 'g', linewidth=2)
                    plt.plot([refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
                             [1.3, 1.3], 'g', linewidth=2)
                    plt.plot([refSpick.getpick() - 0.5, refSpick.getpick() + 0.5],
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
    if mpickP is not None:
    	lpickP = zdat[0].stats.starttime + lpickP
    	epickP = zdat[0].stats.starttime + epickP
    	mpickP = zdat[0].stats.starttime + mpickP

    if mpickS is not None:
    	lpickS = edat[0].stats.starttime + lpickS
    	epickS = edat[0].stats.starttime + epickS
    	mpickS = edat[0].stats.starttime + mpickS
    
    # create dictionary
    # for P phase
    phase = 'P'
    phasepick = {'lpp': lpickP, 'epp': epickP, 'mpp': mpickP, 'spe': Perror, \
                 'snr': SNRP, 'snrdb': SNRPdB, 'weight': Pweight, 'fm': FM}
    picks = {phase: phasepick}
    # add S phase
    phase = 'S'
    phasepick = {'lpp': lpickS, 'epp': epickS, 'mpp': mpickS, 'spe': Serror, \
                 'snr': SNRS, 'snrdb': SNRSdB, 'weight': Sweight, 'fm': None}
    picks[phase] = phasepick
    
    return picks
