#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created August/September 2015.

:author: Ludger Küperkoch / MAGS2 EP3 working group
"""

import matplotlib.pyplot as plt
import numpy as np
from obspy.core import Stream, UTCDateTime
from pylot.core.pick.utils import getsignalwin, crossings_nonzero_all
from pylot.core.util.utils import getPatternLine
from scipy.optimize import curve_fit
from scipy import integrate
from pylot.core.read.data import Data

class Magnitude(object):
    '''
    Superclass for calculating Wood-Anderson peak-to-peak
    amplitudes, local magnitudes, seismic moments
    and moment magnitudes.
    '''

    def __init__(self, wfstream, To, pwin, iplot, NLLocfile=None, \
                 picks=None, rho=None, vp=None, invdir=None):
        '''
        :param: wfstream
        :type: `~obspy.core.stream.Stream

        :param: To, onset time, P- or S phase
        :type:  float

        :param: pwin, pick window [To To+pwin] to get maximum
                peak-to-peak amplitude (WApp) or to calculate
                source spectrum (DCfc) around P onset
        :type:  float

        :param: iplot, no. of figure window for plotting interims results
        :type:  integer

        :param: NLLocfile, name and full path to NLLoc-location file
                needed when calling class MoMw
        :type:  string

        :param: picks, dictionary containing picking results
        :type:  dictionary

        :param: rho [kg/m³], rock density, parameter from autoPyLoT.in
        :type:  integer

        :param: vp [m/s], P-velocity
        :param: integer

        :param: invdir, path to inventory or dataless-SEED file
        :type:  string
        '''

        assert isinstance(wfstream, Stream), "%s is not a stream object" % str(wfstream)

        self.setwfstream(wfstream)
        self.setTo(To)
        self.setpwin(pwin)
        self.setiplot(iplot)
        self.setNLLocfile(NLLocfile)
        self.setrho(rho)
        self.setpicks(picks)
        self.setvp(vp)
        self.setinvdir(invdir)
        self.calcwapp()
        self.calcsourcespec()
        self.run_calcMoMw()


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

    def setNLLocfile(self, NLLocfile):
        self.NLLocfile = NLLocfile

    def getNLLocfile(self):
        return self.NLLocfile

    def setrho(self, rho):
        self.rho = rho

    def getrho(self):
        return self.rho
 
    def setvp(self, vp):
        self.vp = vp

    def getvp(self):
        return self.vp

    def setpicks(self, picks):
        self.picks = picks

    def getpicks(self):
        return self.picks

    def getwapp(self):
        return self.wapp

    def getw0(self):
        return self.w0

    def getfc(self):
        return self.fc

    def setinvdir(self, invdir):
        self.invdir = invdir

    def getinvdir(self):
        return self.invdir

    def getpicdic(self):
        return self.picdic

    def calcwapp(self):
        self.wapp = None

    def calcsourcespec(self):
        self.sourcespek = None

    def run_calcMoMw(self):
        self.pickdic = None

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
    Requires results of class w0fc for calculating plateau w0 
    and corner frequency fc of source spectrum, respectively. Uses 
    subfunction calcMoMw.py. Returns modified dictionary of picks including 
    Dc-value, corner frequency fc, seismic moment Mo and 
    corresponding moment magntiude Mw.
    '''

    def run_calcMoMw(self):

        picks = self.getpicks()
        nllocfile = self.getNLLocfile()
        wfdat = self.getwfstream()
        # get vertical component data only
        zdat = wfdat.select(component="Z")
        if len(zdat) == 0:  # check for other components
            zdat = wfdat.select(component="3")

        for key in picks:
             if picks[key]['P']['weight'] < 4:
                 # select waveform
                 selwf = zdat.select(station=key)
                 # get hypocentral distance of station
                 # from NLLoc-location file
                 if len(key) > 4:
                     Ppattern = '%s  ?    ?    ? P' % key
                 elif len(key) == 4:
                     Ppattern = '%s   ?    ?    ? P' % key
                 elif len(key) < 4:
                     Ppattern = '%s    ?    ?    ? P' % key
                 nllocline = getPatternLine(nllocfile, Ppattern)
                 delta = float(nllocline.split(None)[21])
                 # call subfunction to estimate source spectrum
                 # and to derive w0 and fc 
                 [w0, fc] = calcsourcespec(selwf, picks[key]['P']['mpp'], \
                             self.getiplot(), self.getinvdir())

                 if w0 is not None:
                     # call subfunction to calculate Mo and Mw
                     [Mo, Mw] = calcMoMw(selwf, w0, self.getrho(), self.getvp(), \
                                 delta, self.getinvdir())
                 else:
                     Mo = None
                     Mw = None

                 # add w0, fc, Mo and Mw to dictionary
                 picks[key]['P']['w0'] = w0
                 picks[key]['P']['fc'] = fc
                 picks[key]['P']['Mo'] = Mo
                 picks[key]['P']['Mw'] = Mw
                 self.picdic = picks

def calcMoMw(wfstream, w0, rho, vp, delta, inv):
    '''
    Subfunction of run_calcMoMw to calculate individual
    seismic moments and corresponding moment magnitudes.
    '''

    tr = wfstream[0]

    print("calcMoMw: Calculating seismic moment Mo and moment magnitude Mw for station %s ..." \
           % tr.stats.station)

    # additional common parameters for calculating Mo
    rP = 2 / np.sqrt(15) # average radiation pattern of P waves (Aki & Richards, 1980)
    freesurf = 2.0       # free surface correction, assuming vertical incidence 

    Mo = w0 * 4 * np.pi * rho * np.power(vp, 3) * delta / (rP * freesurf) 

    Mw = np.log10(Mo * 1e07) * 2 / 3 - 10.7 #after Hanks & Kanamori (1979), defined for [dyn*cm]!

    print("calcMoMw: Calculated seismic moment Mo = %e Nm => Mw = %3.1f " % (Mo, Mw))

    return Mo, Mw

        

def calcsourcespec(wfstream, onset, iplot, inventory):
    '''
    Subfunction to calculate the source spectrum and to derive from that the plateau
    (usually called omega0) and the corner frequency assuming Aki's omega-square
    source model. Has to be derived from instrument corrected displacement traces,
    thus restitution and integration necessary!
    '''
    print ("Calculating source spectrum ....")

    fc = None
    w0 = None
    data = Data()
    z_copy = wfstream.copy()

    [corzdat, restflag] = data.restituteWFData(inventory, z_copy)

    if restflag == 1:
        # integrate to displacment
        corintzdat = integrate.cumtrapz(corzdat[0], None, corzdat[0].stats.delta)
        z_copy[0].data = corintzdat
        tr = z_copy[0]
        # get window after P pulse for 
        # calculating source spectrum
        if tr.stats.sampling_rate <= 100:
            winzc = tr.stats.sampling_rate
        elif tr.stats.sampling_rate > 100 and \
              tr.stats.sampling_rate <= 200:
             winzc = 0.5 * tr.stats.sampling_rate
        elif tr.stats.sampling_rate > 200 and \
              tr.stats.sampling_rate <= 400:
             winzc = 0.2 * tr.stats.sampling_rate
        elif tr.stats.sampling_rate > 400:
             winzc = tr.stats.sampling_rate
        tstart = UTCDateTime(tr.stats.starttime)
        tonset = onset.timestamp -tstart.timestamp
        impickP = tonset * tr.stats.sampling_rate
        wfzc = tr.data[impickP : impickP + winzc]
        # get time array
        t = np.arange(0, len(tr) * tr.stats.delta, tr.stats.delta)
        # calculate spectrum using only first cycles of
        # waveform after P onset!
        zc = crossings_nonzero_all(wfzc)
        if np.size(zc) == 0 or len(zc) <= 3:
            print ("Something is wrong with the waveform, "
                   "no zero crossings derived!")
            print ("No calculation of source spectrum possible!")
            plotflag = 0
        else:
            plotflag = 1
            index = min([3, len(zc) - 1])
            calcwin = (zc[index] - zc[0]) * z_copy[0].stats.delta
            iwin = getsignalwin(t, tonset, calcwin)
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
            [w02, fc2] = fitSourceModel(F, YY.real, Fcin, iplot)
 
            # get w0 and fc as median 
            w0 = np.median([w01, w02])
            fc = np.median([fc1, fc2])
            print("w0fc: Using w0-value = %e m/Hz and fc = %f Hz" % (w0, fc))

    if iplot > 1:
        f1 = plt.figure()
        plt.subplot(2,1,1)
        # show displacement in mm
        plt.plot(t, np.multiply(tr, 1000), 'k')
        if plotflag == 1:
            plt.plot(t[iwin], np.multiply(xdat, 1000), 'g')
            plt.title('Seismogram and P Pulse, Station %s-%s' \
                        % (tr.stats.station, tr.stats.channel))
        else:
            plt.title('Seismogram, Station %s-%s' \
                        % (tr.stats.station, tr.stats.channel))
        plt.xlabel('Time since %s' % tr.stats.starttime)
        plt.ylabel('Displacement [mm]')

        if plotflag == 1:
            plt.subplot(2,1,2)
            plt.loglog(f, Y.real, 'k')
            plt.loglog(F, YY.real)
            plt.loglog(F, fit, 'g')
            plt.loglog([fc, fc], [w0/100, w0], 'g') 
            plt.title('Source Spectrum from P Pulse, w0=%e m/Hz, fc=%6.2f Hz' \
                       % (w0, fc))
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Amplitude [m/Hz]')
            plt.grid()
        plt.show()
        raw_input()
        plt.close(f1)

    return w0, fc
     

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







