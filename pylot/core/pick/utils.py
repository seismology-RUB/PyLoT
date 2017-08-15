#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
"""
   Created Mar/Apr 2015
   Collection of helpful functions for manual and automatic picking.

   :author: Ludger Kueperkoch, BESTEC GmbH
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
from obspy.core import Stream, UTCDateTime


def earllatepicker(X, nfac, TSNR, Pick1, iplot=0, verbosity=1, fig=None):
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

    if verbosity == 2:
        print('earllatepicker:')
        print('nfac:', nfac)
        print('Init pick:', Pick1)
        print('TSNR (T_noise, T_gap, T_signal):', TSNR)

    LPick = None
    EPick = None
    PickError = None
    plt_flag = 0
    try:
        iplot = int(iplot)
    except:
        if iplot == True or iplot == 'True':
           iplot = 2
        else:
           iplot = 0

    if verbosity:
        print('earllatepicker: Get earliest and latest possible pick'
              ' relative to most likely pick ...')

    x = X[0].data
    t = np.arange(0, X[0].stats.npts / X[0].stats.sampling_rate,
                  X[0].stats.delta)
    inoise = getnoisewin(t, Pick1, TSNR[0], TSNR[1])
    # get signal window
    isignal = getsignalwin(t, Pick1, TSNR[2])
    # remove mean
    x = x - np.mean(x[inoise])
    # calculate noise level
    nlevel = np.sqrt(np.mean(np.square(x[inoise]))) * nfac
    if verbosity == 2:
        print('x:', x)
        print('t:', t)
        print('x_inoise:', x[inoise])
        print('x_isignal:', x[isignal])
        print('nlevel:', nlevel)

    # get time where signal exceeds nlevel
    ilup, = np.where(x[isignal] > nlevel)
    ildown, = np.where(x[isignal] < -nlevel)
    if not ilup.size and not ildown.size:
        if verbosity:
            print("earllatepicker: Signal lower than noise level!\n"
                  "Skip this trace!")
        return LPick, EPick, PickError
    il = min(np.min(ilup) if ilup.size else float('inf'),
             np.min(ildown) if ildown.size else float('inf'))
    LPick = t[isignal][il]

    # get earliest possible pick

    EPick = np.nan;
    count = 0
    pis = isignal

    # if EPick stays NaN the signal window size will be doubled
    while np.isnan(EPick):
        if count > 0:
            if verbosity:
                print("\nearllatepicker: Doubled signal window size %s time(s) "
                      "because of NaN for earliest pick." % count)
            isigDoubleWinStart = pis[-1] + 1
            isignalDoubleWin = np.arange(isigDoubleWinStart,
                                         isigDoubleWinStart + len(pis))
            if (isigDoubleWinStart + len(pis)) < X[0].data.size:
                pis = np.concatenate((pis, isignalDoubleWin))
            else:
                if verbosity:
                    print("Could not double signal window. Index out of bounds.")
                break
        count += 1
        # determine all zero crossings in signal window (demeaned)
        zc = crossings_nonzero_all(x[pis] - x[pis].mean())
        # calculate mean half period T0 of signal as the average of the
        T0 = np.mean(np.diff(zc)) * X[0].stats.delta  # this is half wave length!
        EPick = Pick1 - T0  # half wavelength as suggested by Diehl et al.

    # get symmetric pick error as mean from earliest and latest possible pick
    # by weighting latest possible pick two times earliest possible pick
    diffti_tl = LPick - Pick1
    diffti_te = Pick1 - EPick
    PickError = symmetrize_error(diffti_te, diffti_tl)

    if iplot > 1:
        if fig == None or fig == 'None':
            fig = plt.figure()  # iplot)
            plt_flag = 1
        ax = fig.add_subplot(111)
        ax.plot(t, x, 'k', label='Data')
        ax.axvspan(t[inoise[0]], t[inoise[-1]], color='y', alpha=0.2, lw=0, label='Noise Window')
        ax.axvspan(t[isignal[0]], t[isignal[-1]], color='b', alpha=0.2, lw=0, label='Signal Window')
        ax.plot([t[0], t[int(len(t)) - 1]], [nlevel, nlevel], '--k', label='Noise Level')
        ax.plot(t[isignal[zc]], np.zeros(len(zc)), '*g',
                markersize=14, label='Zero Crossings')
        ax.plot([t[0], t[int(len(t)) - 1]], [-nlevel, -nlevel], '--k')
        ax.plot([Pick1, Pick1], [max(x), -max(x)], 'b', linewidth=2, label='mpp')
        ax.plot([LPick, LPick], [max(x) / 2, -max(x) / 2], '--k', label='lpp')
        ax.plot([EPick, EPick], [max(x) / 2, -max(x) / 2], '--k', label='epp')
        ax.plot([Pick1 + PickError, Pick1 + PickError],
                [max(x) / 2, -max(x) / 2], 'r--', label='spe')
        ax.plot([Pick1 - PickError, Pick1 - PickError],
                [max(x) / 2, -max(x) / 2], 'r--')
        ax.set_xlabel('Time [s] since %s' % X[0].stats.starttime)
        ax.set_yticks([])
        ax.set_title(
            'Earliest-/Latest Possible/Most Likely Pick & Symmetric Pick Error, %s' %
            X[0].stats.station)
        ax.legend()
        if plt_flag == 1:
            fig.show()
            try: input()
            except SyntaxError: pass
            plt.close(fig)

    return EPick, LPick, PickError


def fmpicker(Xraw, Xfilt, pickwin, Pick, iplot=0, fig=None):
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

    plt_flag = 0
    try:
        iplot = int(iplot)
    except:
        if iplot == True or iplot == 'True':
           iplot = 2
        else:
           iplot = 0

    warnings.simplefilter('ignore', np.RankWarning)

    assert isinstance(Xraw, Stream), "%s is not a stream object" % str(Xraw)
    assert isinstance(Xfilt, Stream), "%s is not a stream object" % str(Xfilt)

    FM = None
    if Pick is not None:
        print("fmpicker: Get first motion (polarity) of onset using unfiltered seismogram...")

        xraw = Xraw[0].data
        xfilt = Xfilt[0].data
        t = np.arange(0, Xraw[0].stats.npts / Xraw[0].stats.sampling_rate,
                      Xraw[0].stats.delta)
        # get pick window
        ipick = np.where(
            (t <= min([Pick + pickwin, len(Xraw[0])])) & (t >= Pick))
        # remove mean
        xraw[ipick] = xraw[ipick] - np.mean(xraw[ipick])
        xfilt[ipick] = xfilt[ipick] - np.mean(xfilt[ipick])

        # get zero crossings after most likely pick
        # initial onset is assumed to be the first zero crossing
        # first from unfiltered trace
        zc1 = []
        zc1.append(Pick)
        index1 = []
        i = 0
        for j in range(ipick[0][1], ipick[0][len(t[ipick]) - 1]):
            i = i + 1
            if xraw[j - 1] <= 0 <= xraw[j]:
                zc1.append(t[ipick][i])
                index1.append(i)
            elif xraw[j - 1] > 0 >= xraw[j]:
                zc1.append(t[ipick][i])
                index1.append(i)
            if len(zc1) == 3:
                break

        if len(zc1) < 3:
            print('fmpicker: Could not determine zero crossings!')
            return

        # if time difference betweeen 1st and 2cnd zero crossing
        # is too short, get time difference between 1st and 3rd
        # to derive maximum
        if zc1[1] - zc1[0] <= Xraw[0].stats.delta:
            li1 = index1[1]
        else:
            li1 = index1[0]
        if np.size(xraw[ipick[0][1]:ipick[0][li1]]) == 0 or len(index1) <= 1:
            print("fmpicker: Onset on unfiltered trace too emergent for first motion determination!")
            P1 = None
        else:
            imax1 = np.argmax(abs(xraw[ipick[0][1]:ipick[0][li1]]))
            if imax1 == 0:
                imax1 = np.argmax(abs(xraw[ipick[0][1]:ipick[0][index1[1]]]))
                if imax1 == 0:
                    print("fmpicker: Zero crossings too close!")
                    print("Skip first motion determination!")
                    return FM

            islope1 = np.where((t >= Pick) & (t <= Pick + t[imax1]))
            # calculate slope as polynomal fit of order 1
            xslope1 = np.arange(0, len(xraw[islope1]), 1)
            P1 = np.polyfit(xslope1, xraw[islope1], 1)
            datafit1 = np.polyval(P1, xslope1)

        # now using filterd trace
        # next zero crossings after most likely pick
        zc2 = []
        zc2.append(Pick)
        index2 = []
        i = 0
        for j in range(ipick[0][1], ipick[0][len(t[ipick]) - 1]):
            i = i + 1
            if xfilt[j - 1] <= 0 <= xfilt[j]:
                zc2.append(t[ipick][i])
                index2.append(i)
            elif xfilt[j - 1] > 0 >= xfilt[j]:
                zc2.append(t[ipick][i])
                index2.append(i)
            if len(zc2) == 3:
                break

        # if time difference betweeen 1st and 2cnd zero crossing
        # is too short, get time difference between 1st and 3rd
        # to derive maximum
        if zc2[1] - zc2[0] <= Xfilt[0].stats.delta:
            li2 = index2[1]
        else:
            li2 = index2[0]
        if np.size(xfilt[ipick[0][1]:ipick[0][li2]]) == 0 or len(index2) <= 1:
            print("fmpicker: Onset on filtered trace too emergent for first motion determination!")
            P2 = None
        else:
            imax2 = np.argmax(abs(xfilt[ipick[0][1]:ipick[0][li2]]))
            if imax2 == 0:
                imax2 = np.argmax(abs(xfilt[ipick[0][1]:ipick[0][index2[1]]]))
                if imax2 == 0:
                    print("fmpicker: Zero crossings too close!")
                    print("Skip first motion determination!")
                    return FM

            islope2 = np.where((t >= Pick) & (t <= Pick + t[imax2]))
            # calculate slope as polynomal fit of order 1
            xslope2 = np.arange(0, len(xfilt[islope2]), 1)
            P2 = np.polyfit(xslope2, xfilt[islope2], 1)
            datafit2 = np.polyval(P2, xslope2)

        # compare results
        if P1 is not None and P2 is not None:
            if P1[0] < 0 and P2[0] < 0:
                FM = 'D'
            elif P1[0] >= 0 > P2[0]:
                FM = '-'
            elif P1[0] < 0 <= P2[0]:
                FM = '-'
            elif P1[0] > 0 and P2[0] > 0:
                FM = 'U'
            elif P1[0] <= 0 < P2[0]:
                FM = '+'
            elif P1[0] > 0 >= P2[0]:
                FM = '+'

        print("fmpicker: Found polarity %s" % FM)

    if iplot > 1:
        if fig == None or fig == 'None':
            fig = plt.figure()  # iplot)
            plt_flag = 1
        ax1 = fig.add_subplot(211)
        ax1.plot(t, xraw, 'k')
        ax1.plot([Pick, Pick], [max(xraw), -max(xraw)], 'b', linewidth=2, label='Pick')
        if P1 is not None:
            ax1.plot(t[islope1], xraw[islope1], label='Slope Window')
            ax1.plot(zc1, np.zeros(len(zc1)), '*g', markersize=14, label='Zero Crossings')
            ax1.plot(t[islope1], datafit1, '--g', linewidth=2)
            ax1.legend()
            ax1.text(Pick + 0.02, max(xraw) / 2, '%s' % FM, fontsize=14)
        ax1.set_yticks([])
        ax1.set_title('First-Motion Determination, %s, Unfiltered Data' % Xraw[
            0].stats.station)

        ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
        ax2.set_title('First-Motion Determination, Filtered Data')
        ax2.plot(t, xfilt, 'k')
        ax2.plot([Pick, Pick], [max(xfilt), -max(xfilt)], 'b',
                 linewidth=2)
        if P2 is not None:
            ax2.plot(t[islope2], xfilt[islope2])
            ax2.plot(zc2, np.zeros(len(zc2)), '*g', markersize=14)
            ax2.plot(t[islope2], datafit2, '--g', linewidth=2)
            ax2.text(Pick + 0.02, max(xraw) / 2, '%s' % FM, fontsize=14)
        ax2.set_xlabel('Time [s] since %s' % Xraw[0].stats.starttime)
        ax2.set_yticks([])
        if plt_flag == 1:
            fig.show()
            try: input()
            except SyntaxError: pass
            plt.close(fig)

    return FM


def crossings_nonzero_all(data):
    pos = data > 0
    npos = ~pos
    return ((pos[:-1] & npos[1:]) | (npos[:-1] & pos[1:])).nonzero()[0]


def symmetrize_error(dte, dtl):
    """
    takes earliest and latest possible pick and returns the symmetrized pick
    uncertainty value
    :param dte: relative lower uncertainty
    :param dtl: relative upper uncertainty
    :return: symmetrized error
    """
    return (dte + 2 * dtl) / 3


def getSNR(X, TSNR, t1, tracenum=0):
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

    SNR = None
    SNRdB = None
    noiselevel = None

    x = X[tracenum].data
    npts = X[tracenum].stats.npts
    sr = X[tracenum].stats.sampling_rate
    dt = X[tracenum].stats.delta
    t = np.arange(0, npts / sr, dt)

    # get noise window
    inoise = getnoisewin(t, t1, TSNR[0], TSNR[1])

    # get signal window
    isignal = getsignalwin(t, t1, TSNR[2])
    if np.size(inoise) < 1:
        print("getSNR: Empty array inoise, check noise window!")
        return SNR, SNRdB, noiselevel

    # demean over entire waveform
    x = x - np.mean(x[inoise])

    # calculate ratios
    noiselevel = np.sqrt(np.mean(np.square(x[inoise])))
    # signallevel = np.sqrt(np.mean(np.square(x[isignal])))

    if np.size(isignal) < 1:
        print("getSNR: Empty array isignal, check signal window!")
        return SNR, SNRdB, noiselevel

    # noiselevel = np.abs(x[inoise]).max()
    signallevel = np.abs(x[isignal]).max()

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

    # get noise window
    inoise, = np.where((t <= max([t1 - tgap, 0])) \
                       & (t >= max([t1 - tnoise - tgap, 0])))
    if np.size(inoise) < 1:
        inoise, = np.where((t >= t[0]) & (t <= t1))
        if np.size(inoise) < 1:
            print("getnoisewin: Empty array inoise, check noise window!")

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

    # get signal window
    isignal, = np.where((t <= min([t1 + tsignal, t[-1]])) \
                        & (t >= t1))
    if np.size(isignal) < 1:
        print("getsignalwin: Empty array isignal, check signal window!")

    return isignal


def getResolutionWindow(snr, extent):
    """
    Number -> Float
    produce the half of the time resolution window width from given SNR
    value
          SNR >= 3    ->  2 sec    HRW
      3 > SNR >= 2    ->  5 sec    MRW
      2 > SNR >= 1.5  -> 10 sec    LRW
    1.5 > SNR         -> 15 sec   VLRW
    see also Diehl et al. 2009

    :parameter: extent, can be 'local', 'regional', 'global'

    >>> getResolutionWindow(0.5)
    7.5
    >>> getResolutionWindow(1.8)
    5.0
    >>> getResolutionWindow(2.3)
    2.5
    >>> getResolutionWindow(4)
    1.0
    >>> getResolutionWindow(2)
    2.5
    """

    res_wins = {
        'regional': {'HRW': 2., 'MRW': 5., 'LRW': 10., 'VLRW': 15.},
        'local': {'HRW': 2., 'MRW': 5., 'LRW': 10., 'VLRW': 15.},
        'global': {'HRW': 40., 'MRW': 100., 'LRW': 200., 'VLRW': 300.}
    }

    if snr:
        if snr < 1.5:
            time_resolution = res_wins[extent]['VLRW']
        elif snr < 2.:
            time_resolution = res_wins[extent]['LRW']
        elif snr < 3.:
            time_resolution = res_wins[extent]['MRW']
        elif snr > 3.:
            time_resolution = res_wins[extent]['HRW']
    else:
        time_resolution = res_wins[extent]['VLRW']

    return time_resolution / 2


def select_for_phase(st, phase):
    '''
    takes a STream object and a phase name and returns that particular component
    which presumably shows the chosen PHASE best

    :param st: stream object containing one or more component[s]
    :type st: `~obspy.core.stream.Stream`
    :param phase: label of the phase for which the stream selection is carried
        out; 'P' or 'S'
    :type phase: str
    :return:
    '''
    from pylot.core.util.defaults import SetChannelComponents

    sel_st = Stream()
    compclass = SetChannelComponents()
    if phase.upper() == 'P':
        comp = 'Z'
        alter_comp = compclass.getCompPosition(comp)
        alter_comp = str(alter_comp[0])
        sel_st += st.select(component=comp)
        sel_st += st.select(component=alter_comp)
    elif phase.upper() == 'S':
        comps = 'NE'
        for comp in comps:
            alter_comp = compclass.getCompPosition(comp)
            alter_comp = str(alter_comp[0])
            sel_st += st.select(component=comp)
            sel_st += st.select(component=alter_comp)
    else:
        raise TypeError('Unknown phase label: {0}'.format(phase))
    return sel_st


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
        print("###############################################")
        print("wadaticheck: Average Vp/Vs ratio before check: %f" % vpvsr)

        checkedPpicks = []
        checkedSpicks = []
        checkedSPtimes = []
        # calculate deviations from Wadati regression
        ii = 0
        ibad = 0
        for key in pickdic:
            if pickdic[key].has_key('SPt'):
                wddiff = abs(pickdic[key]['SPt'] - wdfit[ii])
                ii += 1
                # check, if deviation is larger than adjusted
                if wddiff > dttolerance:
                    # mark onset and downgrade S-weight to 9
                    # (not used anymore)
                    marker = 'badWadatiCheck'
                    pickdic[key]['S']['weight'] = 9
                    ibad += 1
                else:
                    marker = 'goodWadatiCheck'
                    checkedPpick = UTCDateTime(pickdic[key]['P']['mpp'])
                    checkedPpicks.append(checkedPpick.timestamp)
                    checkedSpick = UTCDateTime(pickdic[key]['S']['mpp'])
                    checkedSpicks.append(checkedSpick.timestamp)
                    checkedSPtime = pickdic[key]['S']['mpp'] - pickdic[key]['P']['mpp']
                    checkedSPtimes.append(checkedSPtime)

                pickdic[key]['S']['marked'] = marker

        if len(checkedPpicks) >= 3:
            # calculate new slope
            p2 = np.polyfit(checkedPpicks, checkedSPtimes, 1)
            wdfit2 = np.polyval(p2, checkedPpicks)

            # calculate vp/vs ratio after check
            cvpvsr = p2[0] + 1
            print("wadaticheck: Average Vp/Vs ratio after check: %f" % cvpvsr)
            print("wadatacheck: Skipped %d S pick(s)" % ibad)
        else:
            print("###############################################")
            print("wadatacheck: Not enough checked S-P times available!")
            print("Skip Wadati check!")

        checkedonsets = pickdic

    else:
        print("wadaticheck: Not enough S-P times available for reliable regression!")
        print("Skip wadati check!")
        wfitflag = 1

    # plot results
    if iplot > 0:
        plt.figure()  # iplot)
        f1, = plt.plot(Ppicks, SPtimes, 'ro')
        if wfitflag == 0:
            f2, = plt.plot(Ppicks, wdfit, 'k')
            f3, = plt.plot(checkedPpicks, checkedSPtimes, 'ko')
            f4, = plt.plot(checkedPpicks, wdfit2, 'g')
            plt.title('Wadati-Diagram, %d S-P Times, Vp/Vs(raw)=%5.2f,' \
                      'Vp/Vs(checked)=%5.2f' % (len(SPtimes), vpvsr, cvpvsr))
            plt.legend([f1, f2, f3, f4], ['Skipped S-Picks', 'Wadati 1',
                                          'Reliable S-Picks', 'Wadati 2'], loc='best')
        else:
            plt.title('Wadati-Diagram, %d S-P Times' % len(SPtimes))

        plt.ylabel('S-P Times [s]')
        plt.xlabel('P Times [s]')

    return checkedonsets


def RMS(X):
    '''
    Function returns root mean square of a given array X
    '''
    return np.sqrt(np.sum(np.power(X, 2)) / len(X))


def checksignallength(X, pick, TSNR, minsiglength, nfac, minpercent, iplot=0, fig=None):
    '''
    Function to detect spuriously picked noise peaks.
    Uses RMS trace of all 3 components (if available) to determine,
    how many samples [per cent] after P onset are below certain
    threshold, calculated from noise level times noise factor.

    : param: X, time series (seismogram)
    : type:  `~obspy.core.stream.Stream`

    : param: pick, initial (AIC) P onset time
    : type:  float

    : param: TSNR, length of time windows around initial pick [s]
    : type:  tuple (T_noise, T_gap, T_signal)

    : param: minsiglength, minium required signal length [s] to
             declare pick as P onset
    : type:  float

    : param: nfac, noise factor (nfac * noise level = threshold)
    : type:  float

    : param: minpercent, minimum required percentage of samples
             above calculated threshold
    : type:  float

    : param: iplot, if iplot > 1, results are shown in figure
    : type:  int
    '''

    plt_flag = 0
    try:
        iplot = int(iplot)
    except:
        if iplot == True or iplot == 'True':
           iplot = 2
        else:
           iplot = 0

    assert isinstance(X, Stream), "%s is not a stream object" % str(X)

    print("Checking signal length ...")

    if len(X) > 1:
        # all three components available
        # make sure, all components have equal lengths
        ilen = min([len(X[0].data), len(X[1].data), len(X[2].data)])
        x1 = X[0][0:ilen]
        x2 = X[1][0:ilen]
        x3 = X[2][0:ilen]
        # get RMS trace
        rms = np.sqrt((np.power(x1, 2) + np.power(x2, 2) + np.power(x3, 2)) / 3)
    else:
        x1 = X[0].data
        ilen = len(x1)
        rms = abs(x1)

    t = np.arange(0, ilen / X[0].stats.sampling_rate,
                  X[0].stats.delta)

    # get noise window in front of pick plus saftey gap
    inoise = getnoisewin(t, pick, TSNR[0], TSNR[1])
    # get signal window
    isignal = getsignalwin(t, pick, minsiglength)
    # calculate minimum adjusted signal level
    minsiglevel = max(rms[inoise]) * nfac
    # minimum adjusted number of samples over minimum signal level
    minnum = len(isignal) * minpercent / 100
    # get number of samples above minimum adjusted signal level
    numoverthr = len(np.where(rms[isignal] >= minsiglevel)[0])

    if numoverthr >= minnum:
        print("checksignallength: Signal reached required length.")
        returnflag = 1
    else:
        print("checksignallength: Signal shorter than required minimum signal length!")
        print("Presumably picked noise peak, pick is rejected!")
        print("(min. signal length required: %s s)" % minsiglength)
        returnflag = 0

    if iplot > 1:
        if fig == None or fig == 'None':
            fig = plt.figure()  # iplot)
            plt_flag = 1
        ax = fig.add_subplot(111)
        ax.plot(t, rms, 'k', label='RMS Data')
        ax.axvspan(t[inoise[0]], t[inoise[-1]], color='y', alpha=0.2, lw=0, label='Noise Window')
        ax.axvspan(t[isignal[0]], t[isignal[-1]], color='b', alpha=0.2, lw=0, label='Signal Window')
        ax.plot([t[isignal[0]], t[isignal[len(isignal) - 1]]],
                [minsiglevel, minsiglevel], 'g', linewidth=2, label='Minimum Signal Level')
        ax.plot([pick, pick], [min(rms), max(rms)], 'b', linewidth=2, label='Onset')
        ax.legend()
        ax.set_xlabel('Time [s] since %s' % X[0].stats.starttime)
        ax.set_ylabel('Counts')
        ax.set_title('Check for Signal Length, Station %s' % X[0].stats.station)
        ax.set_yticks([])
        if plt_flag == 1:
            fig.show()
            try: input()
            except SyntaxError: pass
            plt.close(fig)

    return returnflag


def checkPonsets(pickdic, dttolerance, iplot):
    '''
    Function to check statistics of P-onset times: Control deviation from
    median (maximum adjusted deviation = dttolerance) and apply pseudo-
    bootstrapping jackknife.

    : param: pickdic, dictionary containing picks and quality parameters
    : type:  dictionary

    : param: dttolerance, maximum adjusted deviation of P-onset time from
             median of all P onsets
    : type:  float

    : param: iplot, if iplot > 1, Wadati diagram is shown
    : type:  int
    '''

    checkedonsets = pickdic

    # search for good quality P picks
    Ppicks = []
    stations = []
    for key in pickdic:
        if pickdic[key]['P']['weight'] < 4:
            # add P onsets to list
            UTCPpick = UTCDateTime(pickdic[key]['P']['mpp'])
            Ppicks.append(UTCPpick.timestamp)
            stations.append(key)

    # apply jackknife bootstrapping on variance of P onsets
    print("###############################################")
    print("checkPonsets: Apply jackknife bootstrapping on P-onset times ...")
    [xjack, PHI_pseudo, PHI_sub] = jackknife(Ppicks, 'VAR', 1)
    # get pseudo variances smaller than average variances
    # (times safety factor), these picks passed jackknife test
    ij = np.where(PHI_pseudo <= 5 * xjack)
    # these picks did not pass jackknife test
    badjk = np.where(PHI_pseudo > 5 * xjack)
    badjkstations = np.array(stations)[badjk]
    print("checkPonsets: %d pick(s) did not pass jackknife test!" % len(badjkstations))
    print(badjkstations)

    # calculate median from these picks
    pmedian = np.median(np.array(Ppicks)[ij])
    # find picks that deviate less than dttolerance from median
    ii = np.where(abs(np.array(Ppicks)[ij] - pmedian) <= dttolerance)
    jj = np.where(abs(np.array(Ppicks)[ij] - pmedian) > dttolerance)
    igood = ij[0][ii]
    ibad = ij[0][jj]
    goodstations = np.array(stations)[igood]
    badstations = np.array(stations)[ibad]

    print("checkPonsets: %d pick(s) deviate too much from median!" % len(ibad))
    print("checkPonsets: Skipped %d P pick(s) out of %d" % (len(badstations) \
                                                            + len(badjkstations), len(stations)))

    goodmarker = 'goodPonsetcheck'
    badmarker = 'badPonsetcheck'
    badjkmarker = 'badjkcheck'
    for i in range(0, len(goodstations)):
        # mark P onset as checked and keep P weight
        pickdic[goodstations[i]]['P']['marked'] = goodmarker
    for i in range(0, len(badstations)):
        # mark P onset and downgrade P weight to 9
        # (not used anymore)
        pickdic[badstations[i]]['P']['marked'] = badmarker
        pickdic[badstations[i]]['P']['weight'] = 9
    for i in range(0, len(badjkstations)):
        # mark P onset and downgrade P weight to 9
        # (not used anymore)
        pickdic[badjkstations[i]]['P']['marked'] = badjkmarker
        pickdic[badjkstations[i]]['P']['weight'] = 9

    checkedonsets = pickdic

    if iplot > 0:
        p1, = plt.plot(np.arange(0, len(Ppicks)), Ppicks, 'ro', markersize=14)
        if len(badstations) < 1 and len(badjkstations) < 1:
            p2, = plt.plot(np.arange(0, len(Ppicks)), Ppicks, 'go', markersize=14)
        else:
            p2, = plt.plot(igood, np.array(Ppicks)[igood], 'go', markersize=14)
        p3, = plt.plot([0, len(Ppicks) - 1], [pmedian, pmedian], 'g',
                       linewidth=2)
        for i in range(0, len(Ppicks)):
            plt.text(i, Ppicks[i] + 0.01, '{0}'.format(stations[i]))

        plt.xlabel('Number of P Picks')
        plt.ylabel('Onset Time [s] from 1.1.1970')
        plt.legend([p1, p2, p3], ['Skipped P Picks', 'Good P Picks', 'Median'],
                   loc='best')
        plt.title('Jackknifing and Median Tests on P Onsets')

    return checkedonsets


def jackknife(X, phi, h):
    '''
    Function to calculate the Jackknife Estimator for a given quantity,
    special type of boot strapping. Returns the jackknife estimator PHI_jack
    the pseudo values PHI_pseudo and the subgroup parameters PHI_sub.

    : param: X, given quantity
    : type:  list

    : param: phi, chosen estimator, choose between:
             "MED" for median
             "MEA" for arithmetic mean
             "VAR" for variance
    : type:  string

    : param: h, size of subgroups, optinal, default = 1
    : type:  integer
    '''

    PHI_jack = None
    PHI_pseudo = None
    PHI_sub = None

    # determine number of subgroups
    g = len(X) / h

    if type(g) is not int:
        print("jackknife: Cannot divide quantity X in equal sized subgroups!")
        print("Choose another size for subgroups!")
        return PHI_jack, PHI_pseudo, PHI_sub
    else:
        # estimator of undisturbed spot check
        if phi == 'MEA':
            phi_sc = np.mean(X)
        elif phi == 'VAR':
            phi_sc = np.var(X)
        elif phi == 'MED':
            phi_sc = np.median(X)

        # estimators of subgroups
        PHI_pseudo = []
        PHI_sub = []
        for i in range(0, g - 1):
            # subgroup i, remove i-th sample
            xx = X[:]
            del xx[i]
            # calculate estimators of disturbed spot check
            if phi == 'MEA':
                phi_sub = np.mean(xx)
            elif phi == 'VAR':
                phi_sub = np.var(xx)
            elif phi == 'MED':
                phi_sub = np.median(xx)

            PHI_sub.append(phi_sub)
            # pseudo values
            phi_pseudo = g * phi_sc - ((g - 1) * phi_sub)
            PHI_pseudo.append(phi_pseudo)
        # jackknife estimator
        PHI_jack = np.mean(PHI_pseudo)

    return PHI_jack, PHI_pseudo, PHI_sub


def checkZ4S(X, pick, zfac, checkwin, iplot, fig=None):
    '''
    Function to compare energy content of vertical trace with
    energy content of horizontal traces to detect spuriously
    picked S onsets instead of P onsets. Usually, P coda shows
    larger longitudal energy on vertical trace than on horizontal
    traces, where the transversal energy is larger within S coda.
    Be careful: there are special circumstances, where this is not
    the case!

    : param: X, fitered(!) time series, three traces
    : type:  `~obspy.core.stream.Stream`

    : param: pick, initial (AIC) P onset time
    : type:  float

    : param: zfac, factor for threshold determination,
             vertical energy must exceed coda level times zfac
             to declare a pick as P onset
    : type:  float

    : param: checkwin, window length [s] for calculating P-coda
             energy content
    : type:  float

    : param: iplot, if iplot > 1, energy content and threshold
             are shown
    : type:  int
    '''
    
    plt_flag = 0
    try:
        iplot = int(iplot)
    except:
        if iplot == True or iplot == 'True':
           iplot = 2
        else:
           iplot = 0


    assert isinstance(X, Stream), "%s is not a stream object" % str(X)

    print("Check for spuriously picked S onset instead of P onset ...")

    returnflag = 0

    # split components
    zdat = X.select(component="Z")
    if len(zdat) == 0:  # check for other components
        zdat = X.select(component="3")
    edat = X.select(component="E")
    if len(edat) == 0:  # check for other components
        edat = X.select(component="2")
    ndat = X.select(component="N")
    if len(ndat) == 0:  # check for other components
        ndat = X.select(component="1")

    # get earliest time of all 3 traces
    min_t = min(zdat[0].stats.starttime, edat[0].stats.starttime, ndat[0].stats.starttime)

    # generate time arrays for all 3 traces
    tz = np.arange(0, zdat[0].stats.npts / zdat[0].stats.sampling_rate,
                   zdat[0].stats.delta)
    tn = np.arange(0, ndat[0].stats.npts / ndat[0].stats.sampling_rate,
                   ndat[0].stats.delta)
    te = np.arange(0, edat[0].stats.npts / edat[0].stats.sampling_rate,
                   edat[0].stats.delta)

    zdiff = (zdat[0].stats.starttime - min_t)
    ndiff = (ndat[0].stats.starttime - min_t)
    ediff = (edat[0].stats.starttime - min_t)

    # get signal windows
    isignalz = getsignalwin(tz, pick - zdiff, checkwin)
    isignaln = getsignalwin(tn, pick - ndiff, checkwin)
    isignale = getsignalwin(te, pick - ediff, checkwin)

    # calculate RMS of traces
    rmsz = RMS(zdat[0].data[isignalz])
    rmsn = RMS(ndat[0].data[isignaln])
    rmse = RMS(edat[0].data[isignale])

    # calculate threshold
    minsiglevel = (rmsn + rmse) / 2 * zfac

    # vertical P-coda level must exceed horizontal P-coda level
    # zfac times encodalevel
    if rmsz < minsiglevel:
        print("checkZ4S: Maybe S onset? Skip this P pick!")
    else:
        print("checkZ4S: P onset passes checkZ4S test!")
        returnflag = 1

    if iplot > 1:
        rms_dict = {'Z': rmsz,
                    'N': rmsn,
                    'E': rmse}

        traces_dict = {'Z': zdat[0],
                       'N': ndat[0],
                       'E': edat[0]}

        diff_dict = {'Z': zdiff,
                     'N': ndiff,
                     'E': ediff}

        signal_dict = {'Z': isignalz,
                       'N': isignaln,
                       'E': isignale}

        for i, key in enumerate(['Z', 'N', 'E']):
            rms = rms_dict[key]
            trace = traces_dict[key]
            t = np.arange(diff_dict[key], trace.stats.npts / trace.stats.sampling_rate + diff_dict[key],
                          trace.stats.delta)
            if i == 0:
                if fig == None or fig == 'None':
                    fig = plt.figure()  # self.iplot) ### WHY? MP MP
                    plt_flag = 1
                ax1 = fig.add_subplot(3, 1, i + 1)
                ax = ax1
                ax.set_title('CheckZ4S, Station %s' % zdat[0].stats.station)
            else:
                if fig == None or fig == 'None':
                    fig = plt.figure()  # self.iplot) ### WHY? MP MP
                    plt_flag = 1
                ax = fig.add_subplot(3, 1, i + 1, sharex=ax1)
            ax.plot(t, abs(trace.data), color='b', label='abs')
            ax.plot(t, trace.data, color='k')
            name = str(trace.stats.channel) + ': {}'.format(rms)
            ax.plot([pick, pick + checkwin], [rms, rms], 'r', label='RMS {}'.format(name))
            ax.plot([pick, pick], ax.get_ylim(), 'm', label='Pick')
            ax.set_ylabel('Normalized Counts')
            ax.axvspan(pick, pick + checkwin, color='c', alpha=0.2,
                       lw=0)
            ax.legend()
        ax.set_xlabel('Time [s] since %s' % zdat[0].stats.starttime)
        if plt_flag == 1:
            fig.show()
            try: input()
            except SyntaxError: pass
            plt.close(fig)
    return returnflag


def getQualityfromUncertainty(uncertainty, Errors):
    '''Script to transform uncertainty into quality classes 0-4
       regarding adjusted time errors Errors.
    '''

    # set initial quality to 4 (worst) and change only if one condition is hit
    quality = 4

    if uncertainty == None or uncertainty == 'None':
        return quality

    if uncertainty <= Errors[0]:
        quality = 0
    elif (uncertainty > Errors[0]) and \
         (uncertainty < Errors[1]):
        quality = 1
    elif (uncertainty > Errors[1]) and \
         (uncertainty < Errors[2]):
        quality = 2
    elif (uncertainty > Errors[2]) and \
         (uncertainty < Errors[3]):
        quality = 3
    elif uncertainty > Errors[3]:
        quality = 4

    return quality

if __name__ == '__main__':
    import doctest

    doctest.testmod()
