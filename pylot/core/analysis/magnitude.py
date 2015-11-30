#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created August/September 2015.

:author: Ludger Küperkoch / MAGS2 EP3 working group
"""

import matplotlib.pyplot as plt
import numpy as np
from obspy.core import Stream
from pylot.core.pick.utils import getsignalwin
from scipy.optimize import curve_fit

class Magnitude(object):
    '''
    Superclass for calculating Wood-Anderson peak-to-peak
    amplitudes, local magnitudes and moment magnitudes.
    '''

    def __init__(self, wfstream, To, pwin, iplot, w0=None, delta=None, rho=None, vp=None):
        '''
        :param: wfstream
        :type: `~obspy.core.stream.Stream

        :param: To, onset time, P- or S phase
        :type:  float

        :param: pwin, pick window [To To+pwin] to get maximum
                peak-to-peak amplitude (WApp) or to calculate
                source spectrum (DCfc)
        :type:  float

        :param: iplot, no. of figure window for plotting interims results
        :type: integer

        '''

        assert isinstance(wfstream, Stream), "%s is not a stream object" % str(wfstream)

        self.setwfstream(wfstream)
        self.setTo(To)
        self.setpwin(pwin)
        self.setiplot(iplot)
        self.setw0(w0)
        self.setrho(rho)
        self.setdelta(delta)
        self.setvp(vp)
        self.calcwapp()
        self.calcsourcespec()
        self.calcMoMw()


    def getwfstream(self):
        return self.wfstream

    def setwfstream(self, wfstream):
        self.wfstream = wfstream

    def getTo(self):
        return self.To

    def setTo(self, To):
        self.To = To

    def getpwin(self):
        return self.pwin

    def setpwin(self, pwin):
        self.pwin = pwin

    def getiplot(self):
        return self.iplot

    def setiplot(self, iplot):
        self.iplot = iplot

    def setw0(self, w0):
        self.w0 = w0

    def getw0(self):
        return self.w0

    def setrho(self, rho):
        self.rho = rho

    def getrho(self):
        return self.rho
 
    def setvp(self, vp):
        self.vp = vp

    def getvp(self):
        return self.vp

    def setdelta(self, delta):
        self.delta = delta

    def getdelta(self):
        return self.delta

    def getwapp(self):
        return self.wapp

    def getw0(self):
        return self.w0

    def getfc(self):
        return self.fc

    def getMo(self):
        return self.Mo

    def getMw(self):
        return self.Mw

    def calcwapp(self):
        self.wapp = None

    def calcsourcespec(self):
        self.sourcespek = None

    def calcMoMw(self):
        self.Mo = None
        self.Mw = None

class WApp(Magnitude):
    '''
    Method to derive peak-to-peak amplitude as seen on a Wood-Anderson-
    seismograph. Has to be derived from instrument corrected traces!
    '''

    def calcwapp(self):
        print ("Getting Wood-Anderson peak-to-peak amplitude ...")
        print ("Simulating Wood-Anderson seismograph ...")

        self.wapp = None
        stream = self.getwfstream()

        # poles, zeros and sensitivity of WA seismograph
        # (see Uhrhammer & Collins, 1990, BSSA, pp. 702-716)
        paz_wa = {
            'poles': [5.6089 - 5.4978j, -5.6089 - 5.4978j],
            'zeros': [0j, 0j],
            'gain': 2080,
            'sensitivity': 1}

        stream.simulate(paz_remove=None, paz_simulate=paz_wa)

        trH1 = stream[0].data
        trH2 = stream[1].data
        ilen = min([len(trH1), len(trH2)])
        # get RMS of both horizontal components
        sqH = np.sqrt(np.power(trH1[0:ilen], 2) + np.power(trH2[0:ilen], 2))
        # get time array
        th = np.arange(0, len(sqH) * stream[0].stats.delta, stream[0].stats.delta)
        # get maximum peak within pick window
        iwin = getsignalwin(th, self.getTo(), self.getpwin())
        self.wapp = np.max(sqH[iwin])
        print ("Determined Wood-Anderson peak-to-peak amplitude: %f mm") % self.wapp

        if self.getiplot() > 1:
            stream.plot()
            f = plt.figure(2)
            plt.plot(th, sqH)
            plt.plot(th[iwin], sqH[iwin], 'g')
            plt.plot([self.getTo(), self.getTo()], [0, max(sqH)], 'r', linewidth=2)
            plt.title('Station %s, RMS Horizontal Traces, WA-peak-to-peak=%4.1f mm' \
                      % (stream[0].stats.station, self.wapp))
            plt.xlabel('Time [s]')
            plt.ylabel('Displacement [mm]')
            plt.show()
            raw_input()
            plt.close(f)


class M0Mw(Magnitude):
    '''
    Method to calculate seismic moment Mo and moment magnitude Mw.
    Uses class w0fc for calculating plateau wo and corner frequency 
    fc of source spectrum, respectively.
    '''

    def calcMoMw(self):

        stream = self.getwfstream()
        tr = stream[0]

        print("Calculating seismic moment Mo and moment magnitude Mw for station %s ..." \
               % tr.stats.station)

        # additional common parameters for calculating Mo
        rP = 0.52      # average radiation pattern of P waves (Aki & Richards, 1980)
        freesurf = 2.0 # free surface correction, assuming vertical incidence 

        self.Mo = (self.getw0() * 4 * np.pi * self.getrho() * np.power(self.getvp(), 3) * \
                  self.getdelta()) / (rP * freesurf) 

        self.Mw = 2/3 * np.log10(self.Mo) - 6
        print("MoMw: Calculated seismic moment Mo = %e Nm => Mw = %3.1f " % (self.Mo, self.Mw))

        
class w0fc(Magnitude):
    '''
    Method to calculate the source spectrum and to derive from that the plateau
    (usually called omega0) and the corner frequency assuming Aki's omega-square
    source model. Has to be derived from instrument corrected displacement traces!
    '''

    def calcsourcespec(self):
        print ("Calculating source spectrum ....")

        self.w0 = None # DC-value
        self.fc = None # corner frequency

        stream = self.getwfstream()
        tr = stream[0]

        # get time array
        t = np.arange(0, len(tr) * tr.stats.delta, tr.stats.delta)
        iwin = getsignalwin(t, self.getTo(), self.getpwin())
        xdat = tr.data[iwin]

        # fft
        fny = tr.stats.sampling_rate / 2
        l = len(xdat) / tr.stats.sampling_rate
        n = tr.stats.sampling_rate * l # number of fft bins after Bath
        # find next power of 2 of data length
        m = pow(2, np.ceil(np.log(len(xdat)) / np.log(2)))
        N = int(np.power(m, 2))
        y = tr.stats.delta * np.fft.fft(xdat, N)
        Y = abs(y[: N/2])
        L = (N - 1) / tr.stats.sampling_rate
        f = np.arange(0, fny, 1/L)

        # remove zero-frequency and frequencies above
        # corner frequency of seismometer (assumed
        # to be 100 Hz)
        fi = np.where((f >= 1) & (f < 100))
        F = f[fi]
        YY = Y[fi]
        # get plateau (DC value) and corner frequency
        # initial guess of plateau
        w0in = np.mean(YY[0:100])
        # initial guess of corner frequency
        # where spectral level reached 50% of flat level
        iin = np.where(YY >= 0.5 * w0in)
        Fcin = F[iin[0][np.size(iin) - 1]]

        # use of implicit scipy otimization function
        fit = synthsourcespec(F, w0in, Fcin)
        [optspecfit, pcov] = curve_fit(synthsourcespec, F, YY.real, [w0in, Fcin])
        w01 = optspecfit[0]
        fc1 = optspecfit[1]
        print ("w0fc: Determined w0-value: %e m/Hz, \n"
               "Determined corner frequency: %f Hz" % (w01, fc1))
        
        # use of conventional fitting 
        [w02, fc2] = fitSourceModel(F, YY.real, Fcin, self.getiplot())
 
        # get w0 and fc as median 
        self.w0 = np.median([w01, w02])
        self.fc = np.median([fc1, fc2])
        print("w0fc: Using w0-value = %e m/Hz and fc = %f Hz" % (self.w0, self.fc))

	if self.getiplot() > 1:
            f1 = plt.figure()
            plt.subplot(2,1,1)
            # show displacement in mm
            plt.plot(t, np.multiply(tr, 1000), 'k')
            plt.plot(t[iwin], np.multiply(xdat, 1000), 'g')
            plt.title('Seismogram and P pulse, station %s' % tr.stats.station)
            plt.xlabel('Time since %s' % tr.stats.starttime)
            plt.ylabel('Displacement [mm]')

            plt.subplot(2,1,2)
            plt.loglog(f, Y.real, 'k')
            plt.loglog(F, YY.real)
            plt.loglog(F, fit, 'g')
            plt.loglog([self.fc, self.fc], [self.w0/100, self.w0], 'g') 
            plt.title('Source Spectrum from P Pulse, w0=%e m/Hz, fc=%6.2f Hz' \
                       % (self.w0, self.fc))
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Amplitude [m/Hz]')
            plt.grid()
            plt.show()
            raw_input()
            plt.close(f1)



def synthsourcespec(f, omega0, fcorner):
    '''
    Calculates synthetic source spectrum from given plateau and corner
    frequency assuming Akis omega-square model.

    :param: f, frequencies
    :type:  array

    :param: omega0, DC-value (plateau) of source spectrum
    :type:  float

    :param: fcorner, corner frequency of source spectrum
    :type:  float
    '''

    #ssp = omega0 / (pow(2, (1 + f / fcorner)))
    ssp = omega0 / (1 + pow(2, (f / fcorner)))

    return ssp


def fitSourceModel(f, S, fc0, iplot):
    '''
    Calculates synthetic source spectrum by varying corner frequency fc.
    Returns best approximated plateau omega0 and corner frequency, i.e. with least 
    common standard deviations. 

    :param:  f, frequencies
    :type:   array

    :param:  S, observed source spectrum
    :type:   array

    :param:  fc0, initial corner frequency
    :type:   float
    '''

    w0 =  []
    stdw0 = []
    fc = []
    stdfc = []
    STD = []

    # get window around initial corner frequency for trials
    fcstopl = fc0 - max(1, len(f) / 10)
    il = np.argmin(abs(f-fcstopl))
    fcstopl = f[il]
    fcstopr = fc0 + min(len(f), len(f) /10) 
    ir = np.argmin(abs(f-fcstopr))
    fcstopr = f[ir]
    iF = np.where((f >= fcstopl) & (f <= fcstopr))

    # vary corner frequency around initial point
    for i in range(il, ir): 
        FC = f[i]
        indexdc = np.where((f > 0 ) & (f <= FC))
        dc = np.mean(S[indexdc])
        stddc = np.std(dc - S[indexdc])
        w0.append(dc)
        stdw0.append(stddc)
        fc.append(FC)
        # slope
        indexfc = np.where((f >= FC) & (f <= fcstopr))
        yi = dc/(1+(f[indexfc]/FC)**2)
        stdFC = np.std(yi - S[indexfc])
        stdfc.append(stdFC)
        STD.append(stddc + stdFC)

    # get best found w0 anf fc from minimum
    fc = fc[np.argmin(STD)]
    w0 = w0[np.argmin(STD)]
    print("fitSourceModel: best fc: %fHz, best w0: %e m/Hz" \
           % (fc, w0))

    if iplot > 1:
        plt.figure(iplot)
        plt.loglog(f, S, 'k')
        plt.loglog([f[0], fc], [w0, w0], 'g')
        plt.loglog([fc, fc], [w0/100, w0], 'g') 
        plt.title('Calculated Source Spectrum, Omega0=%e m/Hz, fc=%6.2f Hz' \
                   % (w0, fc)) 
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [m/Hz]')
        plt.grid()
        plt.figure(iplot+1)
        plt.subplot(311)
        plt.plot(f[il:ir], STD,'*')
        plt.title('Common Standard Deviations')
        plt.xticks([])
        plt.subplot(312)
        plt.plot(f[il:ir], stdw0,'*')
        plt.title('Standard Deviations of w0-Values')
        plt.xticks([])
        plt.subplot(313)
        plt.plot(f[il:ir],stdfc,'*')
        plt.title('Standard Deviations of Corner Frequencies')
        plt.xlabel('Corner Frequencies [Hz]')
        plt.show()
        raw_input()
        plt.close()

    return w0, fc







