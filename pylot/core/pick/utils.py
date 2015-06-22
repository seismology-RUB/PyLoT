#!/usr/bin/env python
#
# -*- coding: utf-8 -*-
"""
   Created Mar/Apr 2015
   Collection of helpful functions for manual and automatic picking.

   :author: Ludger Kueperkoch / MAGS2 EP3 working group
"""

import numpy as np
import matplotlib.pyplot as plt
from obspy.core import Stream, UTCDateTime
import warnings


def earllatepicker(X, nfac, TSNR, Pick1, iplot=None):
    '''
    Function to derive earliest and latest possible pick after Diehl & Kissling (2009)
    as reasonable uncertainties. Latest possible pick is based on noise level,
    earliest possible pick is half a signal wavelength in front of most likely
    pick given by PragPicker or manually set by analyst. Most likely pick
    (initial pick Pick1) must be given.

    :param: X, time series (seismogram)
    :type:  `~obspy.core.stream.Stream`

    :param: nfac (noise factor), nfac times noise level to calculate latest possible pick
    :type: int

    :param: TSNR, length of time windows around pick used to determine SNR [s]
    :type: tuple (T_noise, T_gap, T_signal)

    :param: Pick1, initial (most likely) onset time, starting point for earllatepicker
    :type: float

    :param: iplot, if given, results are plotted in figure(iplot)
    :type: int
    '''

    assert isinstance(X, Stream), "%s is not a stream object" % str(X)

    LPick = None
    EPick = None
    PickError = None
    print 'earllatepicker: Get earliest and latest possible pick relative to most likely pick ...'

    x = X[0].data
    t = np.arange(0, X[0].stats.npts / X[0].stats.sampling_rate,
                  X[0].stats.delta)
    # get latest possible pick
    #get noise window
    inoise = getnoisewin(t, Pick1, TSNR[0], TSNR[1])
    #get signal window
    isignal = getsignalwin(t, Pick1, TSNR[2])
    #calculate noise level
    nlevel = np.sqrt(np.mean(np.square(x[inoise]))) * nfac
    #get time where signal exceeds nlevel
    ilup, = np.where(x[isignal] > nlevel)
    ildown, = np.where(x[isignal] < -nlevel)
    if not ilup.size and not ildown.size:
        raise ValueError('earllatepicker: Signal lower than noise level')
    il = min(np.min(ilup) if ilup.size else float('inf'),
             np.min(ildown) if ildown.size else float('inf'))
    LPick = t[isignal][il]

    #get earliest possible pick

    #determine all zero crossings in signal window
    zc = crossings_nonzero_all(x[isignal])
    #calculate mean half period T0 of signal as the average of the
    T0 = np.mean(np.diff(zc)) * X[0].stats.delta  #this is half wave length!
    #T0/4 is assumed as time difference between most likely and earliest possible pick!
    EPick = Pick1 - T0 / 2

    #get symmetric pick error as mean from earliest and latest possible pick
    #by weighting latest possible pick two times earliest possible pick
    diffti_tl = LPick - Pick1
    diffti_te = Pick1 - EPick
    PickError = (diffti_te + 2 * diffti_tl) / 3

    if iplot > 1:
        p = plt.figure(iplot)
        p1, = plt.plot(t, x, 'k')
        p2, = plt.plot(t[inoise], x[inoise])
        p3, = plt.plot(t[isignal], x[isignal], 'r')
        p4, = plt.plot([t[0], t[int(len(t)) - 1]], [nlevel, nlevel], '--k')
        p5, = plt.plot(t[isignal[0][zc]], np.zeros(len(zc)), '*g', markersize=14)
        plt.legend([p1, p2, p3, p4, p5],
                   ['Data', 'Noise Window', 'Signal Window', 'Noise Level',
                    'Zero Crossings'], \
                   loc='best')
        plt.plot([t[0], t[int(len(t)) - 1]], [-nlevel, -nlevel], '--k')
        plt.plot([Pick1, Pick1], [max(x), -max(x)], 'b', linewidth=2)
        plt.plot([LPick, LPick], [max(x) / 2, -max(x) / 2], '--k')
        plt.plot([EPick, EPick], [max(x) / 2, -max(x) / 2], '--k')
        plt.plot([Pick1 + PickError, Pick1 + PickError],
                 [max(x) / 2, -max(x) / 2], 'r--')
        plt.plot([Pick1 - PickError, Pick1 - PickError],
                 [max(x) / 2, -max(x) / 2], 'r--')
        plt.xlabel('Time [s] since %s' % X[0].stats.starttime)
        plt.yticks([])
        ax = plt.gca()
        ax.set_xlim([t[inoise[0][0]] - 2, t[isignal[0][len(isignal) - 1]] + 3])
        plt.title(
            'Earliest-/Latest Possible/Most Likely Pick & Symmetric Pick Error, %s' %
            X[0].stats.station)
        plt.show()
        raw_input()
        plt.close(p)

    return EPick, LPick, PickError


def fmpicker(Xraw, Xfilt, pickwin, Pick, iplot=None):
    '''
    Function to derive first motion (polarity) of given phase onset Pick.
    Calculation is based on zero crossings determined within time window pickwin
    after given onset time.

    :param: Xraw, unfiltered time series (seismogram)
    :type:  `~obspy.core.stream.Stream`

    :param: Xfilt, filtered time series (seismogram)
    :type:  `~obspy.core.stream.Stream`

    :param: pickwin, time window after onset Pick within zero crossings are calculated
    :type: float

    :param: Pick, initial (most likely) onset time, starting point for fmpicker
    :type: float

    :param: iplot, if given, results are plotted in figure(iplot)
    :type: int
    '''

    warnings.simplefilter('ignore', np.RankWarning)

    assert isinstance(Xraw, Stream), "%s is not a stream object" % str(Xraw)
    assert isinstance(Xfilt, Stream), "%s is not a stream object" % str(Xfilt)

    FM = None
    if Pick is not None:
        print 'fmpicker: Get first motion (polarity) of onset using unfiltered seismogram...'

        xraw = Xraw[0].data
        xfilt = Xfilt[0].data
        t = np.arange(0, Xraw[0].stats.npts / Xraw[0].stats.sampling_rate,
                      Xraw[0].stats.delta)
        # get pick window
        ipick = np.where(
            (t <= min([Pick + pickwin, len(Xraw[0])])) & (t >= Pick))
        #remove mean
        xraw[ipick] = xraw[ipick] - np.mean(xraw[ipick])
        xfilt[ipick] = xfilt[ipick] - np.mean(xfilt[ipick])

        #get next zero crossing after most likely pick
        #initial onset is assumed to be the first zero crossing
        #first from unfiltered trace
        zc1 = []
        zc1.append(Pick)
        index1 = []
        i = 0
        for j in range(ipick[0][1], ipick[0][len(t[ipick]) - 1]):
            i = i + 1
            if xraw[j - 1] <= 0 and xraw[j] >= 0:
                zc1.append(t[ipick][i])
                index1.append(i)
            elif xraw[j - 1] > 0 and xraw[j] <= 0:
                zc1.append(t[ipick][i])
                index1.append(i)
            if len(zc1) == 3:
                break

        #if time difference betweeen 1st and 2cnd zero crossing
        #is too short, get time difference between 1st and 3rd
        #to derive maximum
        if zc1[1] - zc1[0] <= Xraw[0].stats.delta:
            li1 = index1[1]
        else:
            li1 = index1[0]
        if np.size(xraw[ipick[0][1]:ipick[0][li1]]) == 0:
            print 'earllatepicker: Onset on unfiltered trace too emergent for first motion determination!'
            P1 = None
        else:
            imax1 = np.argmax(abs(xraw[ipick[0][1]:ipick[0][li1]]))
            islope1 = np.where((t >= Pick) & (t <= Pick + t[imax1]))
            #calculate slope as polynomal fit of order 1
            xslope1 = np.arange(0, len(xraw[islope1]), 1)
            P1 = np.polyfit(xslope1, xraw[islope1], 1)
            datafit1 = np.polyval(P1, xslope1)

        #now using filterd trace
        #next zero crossing after most likely pick
        zc2 = []
        zc2.append(Pick)
        index2 = []
        i = 0
        for j in range(ipick[0][1], ipick[0][len(t[ipick]) - 1]):
            i = i + 1
            if xfilt[j - 1] <= 0 and xfilt[j] >= 0:
                zc2.append(t[ipick][i])
                index2.append(i)
            elif xfilt[j - 1] > 0 and xfilt[j] <= 0:
                zc2.append(t[ipick][i])
                index2.append(i)
            if len(zc2) == 3:
                break

        #if time difference betweeen 1st and 2cnd zero crossing
        #is too short, get time difference between 1st and 3rd
        #to derive maximum
        if zc2[1] - zc2[0] <= Xfilt[0].stats.delta:
            li2 = index2[1]
        else:
            li2 = index2[0]
        if np.size(xfilt[ipick[0][1]:ipick[0][li2]]) == 0:
            print 'earllatepicker: Onset on filtered trace too emergent for first motion determination!'
            P2 = None
        else:
            imax2 = np.argmax(abs(xfilt[ipick[0][1]:ipick[0][li2]]))
            islope2 = np.where((t >= Pick) & (t <= Pick + t[imax2]))
            #calculate slope as polynomal fit of order 1
            xslope2 = np.arange(0, len(xfilt[islope2]), 1)
            P2 = np.polyfit(xslope2, xfilt[islope2], 1)
            datafit2 = np.polyval(P2, xslope2)

        #compare results
        if P1 is not None and P2 is not None:
            if P1[0] < 0 and P2[0] < 0:
                FM = 'D'
            elif P1[0] >= 0 and P2[0] < 0:
                FM = '-'
            elif P1[0] < 0 and P2[0] >= 0:
                FM = '-'
            elif P1[0] > 0 and P2[0] > 0:
                FM = 'U'
            elif P1[0] <= 0 and P2[0] > 0:
                FM = '+'
            elif P1[0] > 0 and P2[0] <= 0:
                FM = '+'

    if iplot > 1:
        plt.figure(iplot)
        plt.subplot(2, 1, 1)
        plt.plot(t, xraw, 'k')
        p1, = plt.plot([Pick, Pick], [max(xraw), -max(xraw)], 'b', linewidth=2)
        if P1 is not None:
            p2, = plt.plot(t[islope1], xraw[islope1])
            p3, = plt.plot(zc1, np.zeros(len(zc1)), '*g', markersize=14)
            p4, = plt.plot(t[islope1], datafit1, '--g', linewidth=2)
            plt.legend([p1, p2, p3, p4],
                       ['Pick', 'Slope Window', 'Zero Crossings', 'Slope'], \
                       loc='best')
            plt.text(Pick + 0.02, max(xraw) / 2, '%s' % FM, fontsize=14)
            ax = plt.gca()
            ax.set_xlim(
                [t[islope1[0][0]] - 0.1, t[islope1[0][len(islope1) - 1]] + 0.3])
        plt.yticks([])
        plt.title('First-Motion Determination, %s, Unfiltered Data' % Xraw[
            0].stats.station)

        plt.subplot(2, 1, 2)
        plt.title('First-Motion Determination, Filtered Data')
        plt.plot(t, xfilt, 'k')
        p1, = plt.plot([Pick, Pick], [max(xfilt), -max(xfilt)], 'b',
                       linewidth=2)
        if P2 is not None:
            p2, = plt.plot(t[islope2], xfilt[islope2])
            p3, = plt.plot(zc2, np.zeros(len(zc2)), '*g', markersize=14)
            p4, = plt.plot(t[islope2], datafit2, '--g', linewidth=2)
            plt.text(Pick + 0.02, max(xraw) / 2, '%s' % FM, fontsize=14)
            ax = plt.gca()
            ax.set_xlim(
                [t[islope2[0][0]] - 0.1, t[islope2[0][len(islope2) - 1]] + 0.3])
        plt.xlabel('Time [s] since %s' % Xraw[0].stats.starttime)
        plt.yticks([])
        plt.show()
        raw_input()
        plt.close(iplot)

    return FM

def crossings_nonzero_all(data):
    pos = data > 0
    npos = ~pos
    return ((pos[:-1] & npos[1:]) | (npos[:-1] & pos[1:])).nonzero()[0]

def getSNR(X, TSNR, t1):
    '''
    Function to calculate SNR of certain part of seismogram relative to
    given time (onset) out of given noise and signal windows. A safety gap
    between noise and signal part can be set. Returns SNR and SNR [dB] and
    noiselevel.

    :param: X, time series (seismogram)
    :type:  `~obspy.core.stream.Stream`

    :param: TSNR, length of time windows [s] around t1 (onset) used to determine SNR
    :type: tuple (T_noise, T_gap, T_signal)

    :param: t1, initial time (onset) from which noise and signal windows are calculated
    :type: float
    '''

    assert isinstance(X, Stream), "%s is not a stream object" % str(X)

    x = X[0].data
    t = np.arange(0, X[0].stats.npts / X[0].stats.sampling_rate,
                  X[0].stats.delta)

    # get noise window
    inoise = getnoisewin(t, t1, TSNR[0], TSNR[1])

    #get signal window
    isignal = getsignalwin(t, t1, TSNR[2])
    if np.size(inoise) < 1:
        print 'getSNR: Empty array inoise, check noise window!'
        return
    elif np.size(isignal) < 1:
        print 'getSNR: Empty array isignal, check signal window!'
        return

    #calculate ratios
    noiselevel = np.sqrt(np.mean(np.square(x[inoise])))
    signallevel = np.sqrt(np.mean(np.square(x[isignal])))
    SNR = signallevel / noiselevel
    SNRdB = 10 * np.log10(SNR)

    return SNR, SNRdB, noiselevel


def getnoisewin(t, t1, tnoise, tgap):
    '''
    Function to extract indeces of data out of time series for noise calculation.
    Returns an array of indeces.

    :param: t, array of time stamps
    :type:  numpy array

    :param: t1, time from which relativ to it noise window is extracted
    :type: float

    :param: tnoise, length of time window [s] for noise part extraction
    :type: float

    :param: tgap, safety gap between t1 (onset) and noise window to
            ensure, that noise window contains no signal
    :type: float
    '''

    inoise = None
    # get noise window
    inoise = np.where((t <= max([t1 - tgap, 0])) \
                      & (t >= max([t1 - tnoise - tgap, 0])))
    if np.size(inoise) < 1:
        print 'getnoisewin: Empty array inoise, check noise window!'

    return inoise


def getsignalwin(t, t1, tsignal):
    '''
    Function to extract data out of time series for signal level calculation.
    Returns an array of indeces.

    :param: t, array of time stamps
    :type:  numpy array

    :param: t1, time from which relativ to it signal window is extracted
    :type: float

    :param: tsignal, length of time window [s] for signal level calculation
    :type: float
    '''

    inoise = None
    # get signal window
    isignal = np.where((t <= min([t1 + tsignal, len(t)])) \
                       & (t >= t1))
    if np.size(isignal) < 1:
        print 'getsignalwin: Empty array isignal, check signal window!'

    return isignal

def wadaticheck(pickdic, dttolerance, iplot):
    '''
    Function to calculate Wadati-diagram from given P and S onsets in order
    to detect S pick outliers. If a certain S-P time deviates by dttolerance
    from regression of S-P time the S pick is marked and down graded.

    : param: pickdic, dictionary containing picks and quality parameters
    : type:  dictionary
    
    : param: dttolerance, maximum adjusted deviation of S-P time from 
             S-P time regression
    : type:  float

    : param: iplot, if iplot > 1, Wadati diagram is shown
    : type:  int
    '''

    checkedonsets = pickdic

    # search for good quality picks and calculate S-P time
    Ppicks = []
    Spicks = []
    SPtimes = []
    for key in pickdic:
        if pickdic[key]['P']['weight'] < 4 and pickdic[key]['S']['weight'] < 4:
           # calculate S-P time
           spt = pickdic[key]['S']['mpp'] - pickdic[key]['P']['mpp']
           # add S-P time to dictionary
           pickdic[key]['SPt'] = spt
           # add P onsets and corresponding S-P times to list
           UTCPpick = UTCDateTime(pickdic[key]['P']['mpp'])
           UTCSpick = UTCDateTime(pickdic[key]['S']['mpp'])
           Ppicks.append(UTCPpick.timestamp)
           Spicks.append(UTCSpick.timestamp)
           SPtimes.append(spt)


    if len(SPtimes) >= 3:
    	# calculate slope 
    	p1 = np.polyfit(Ppicks, SPtimes, 1)
    	wdfit = np.polyval(p1, Ppicks)
        wfitflag = 0
     
        # calculate vp/vs ratio before check
        vpvsr = p1[0] + 1
        print 'wadaticheck: Average Vp/Vs ratio before check:', vpvsr
 
        checkedPpicks = []
        checkedSpicks = []
        checkedSPtimes = []
        # calculate deviations from Wadati regression
        for key in pickdic:
            if pickdic[key].has_key('SPt'):
                ii = 0
            	wddiff = abs(pickdic[key]['SPt'] - wdfit[ii])    	
                ii += 1
                # check, if deviation is larger than adjusted
                if wddiff >= dttolerance:
                    # mark onset and downgrade S-weight to 9 
                    # (not used anymore)
                    marker = 'badWadatiCheck'
                    pickdic[key]['S']['weight'] = 9
                else:
                    marker = 'goodWadatiCheck'
                    checkedPpick =  UTCDateTime(pickdic[key]['P']['mpp'])
                    checkedPpicks.append(checkedPpick.timestamp)
                    checkedSpick = UTCDateTime(pickdic[key]['S']['mpp'])
                    checkedSpicks.append(checkedSpick.timestamp)
                    checkedSPtime = pickdic[key]['S']['mpp'] - pickdic[key]['P']['mpp']
                    checkedSPtimes.append(checkedSPtime)

                pickdic[key]['S']['marked'] = marker


    	# calculate new slope 
    	p2 = np.polyfit(checkedPpicks, checkedSPtimes, 1)
    	wdfit2 = np.polyval(p2, checkedPpicks)

        # calculate vp/vs ratio after check
        cvpvsr = p2[0] + 1
        print 'wadaticheck: Average Vp/Vs ratio after check:', cvpvsr

        checkedonsets = pickdic
 
    else:
    	print 'wadaticheck: Not enough S-P times available for reliable regression!'
        print 'Skip wadati check!'
        wfitflag = 1

    # plot results
    if iplot > 1:
    	plt.figure(iplot)
    	f1, = plt.plot(Ppicks, SPtimes, 'ro')
        if wfitflag == 0:
        	f2, = plt.plot(Ppicks, wdfit, 'k')
                f3, = plt.plot(checkedPpicks, checkedSPtimes, 'ko')
                f4, = plt.plot(checkedPpicks, wdfit2, 'g')
        plt.ylabel('S-P Times [s]')
        plt.xlabel('P Times [s]')
        plt.title('Wadati-Diagram, %d S-P Times, Vp/Vs(raw)=%5.2f, Vp/Vs(checked)=%5.2f' \
                                                     % (len(SPtimes), vpvsr, cvpvsr))
        plt.legend([f1, f2, f3, f4], ['Skipped S-Picks', 'Wadati 1', 'Reliable S-Picks', \
                                       'Wadati 2'], loc='best')
        plt.show()
        raw_input()
        plt.close(iplot)

    return checkedonsets
