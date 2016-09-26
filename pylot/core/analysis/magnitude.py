#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created autumn/winter 2015.

:author: Ludger Küperkoch / MAGS2 EP3 working group
"""
import os
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import Stream, UTCDateTime
from obspy.core.event import Magnitude as ObspyMagnitude
from obspy.geodetics import degrees2kilometers

from pylot.core.pick.utils import getsignalwin, crossings_nonzero_all, select_for_phase
from pylot.core.util.utils import getPatternLine
from scipy.optimize import curve_fit
from scipy import integrate, signal
from pylot.core.util.dataprocessing import restitute_data, read_metadata
from pylot.core.util.utils import common_range, fit_curve

def richter_magnitude_scaling(delta):
    relation = np.loadtxt(os.path.join(os.path.expanduser('~'),
                            '.pylot', 'richter_scaling.data'))
    # prepare spline interpolation to calculate return value
    func, params = fit_curve(relation[:,0], relation[:, 1])
    return func(delta, params)


class Magnitude(object):
    """
    Base class object for Magnitude calculation within PyLoT.
    """

    def __init__(self, stream):
        self._stream = stream
        self._magnitudes = dict()

    @property
    def stream(self):
        return self._stream

    @stream.setter
    def stream(self, value):
        self._stream = value

    @property
    def magnitudes(self):
        return self._magnitudes

    @magnitudes.setter
    def magnitudes(self, value):
        """
        takes a tuple and saves the key value pair to private
        attribute _magnitudes
        :param value: station, magnitude value pair
        :type value: tuple or list
        :return:
        """
        station, magnitude = value
        self._magnitudes[station] = magnitude

    def get(self):
        return self.magnitudes

class Magnitude(object):
    '''
    Superclass for calculating Wood-Anderson peak-to-peak
    amplitudes, local magnitudes, source spectra, seismic moments
    and moment magnitudes.
    '''

    def __init__(self, wfstream, t0, pwin, iplot, NLLocfile=None, \
                 picks=None, rho=None, vp=None, Qp=None, invdir=None):
        '''
        :param: wfstream
        :type: `~obspy.core.stream.Stream

        :param: t0, onset time, P- or S phase
        :type:  float

        :param: pwin, pick window [t0 t0+pwin] to get maximum
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

        :param: invdir, name and path to inventory or dataless-SEED file
        :type:  string
        '''

        assert isinstance(wfstream, Stream), "%s is not a stream object" % str(wfstream)

        self._stream = wfstream
        self._invdir = invdir
        self._t0 = t0
        self._pwin = pwin
        self._iplot = iplot
        self.setNLLocfile(NLLocfile)
        self.setrho(rho)
        self.setpicks(picks)
        self.setvp(vp)
        self.setQp(Qp)
        self.calcwapp()
        self.calcsourcespec()
        self.run_calcMoMw()

    @property
    def stream(self):
        return self._stream

    @stream.setter
    def stream(self, wfstream):
        self._stream = wfstream

    @property
    def t0(self):
        return self._t0

    @t0.setter
    def t0(self, value):
        self._t0 = value

    @property
    def invdir(self):
        return self._invdir

    @invdir.setter
    def invdir(self, value):
        self._invdir = value

    @property
    def pwin(self):
        return self._pwin

    @pwin.setter
    def pwin(self, value):
        self._pwin = value

    @property
    def plot_flag(self):
        return self.iplot

    @plot_flag.setter
    def plot_flag(self, value):
        self._iplot = value

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

    def setQp(self, Qp):
        self.Qp = Qp

    def getQp(self):
        return self.Qp

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

    def get_metadata(self):
        return read_metadata(self.invdir)

    def plot(self):
        pass

    def getpicdic(self):
        return self.picdic

    def calcwapp(self):
        self.wapp = None

    def calcsourcespec(self):
        self.sourcespek = None

    def run_calcMoMw(self):
        self.pickdic = None


class RichterMagnitude(Magnitude):
    '''
    Method to derive peak-to-peak amplitude as seen on a Wood-Anderson-
    seismograph. Has to be derived from instrument corrected traces!
    '''

    # poles, zeros and sensitivity of WA seismograph
    # (see Uhrhammer & Collins, 1990, BSSA, pp. 702-716)
    _paz = {
        'poles': [5.6089 - 5.4978j, -5.6089 - 5.4978j],
        'zeros': [0j, 0j],
        'gain': 2080,
        'sensitivity': 1
    }

    def __init__(self, stream, t0, calc_win):
        super(RichterMagnitude, self).__init__(stream)

        self._t0 = t0
        self._calc_win = calc_win

    @property
    def t0(self):
        return self._t0

    @t0.setter
    def t0(self, value):
        self._t0 = value

    @property
    def calc_win(self):
        return self._calc_win

    @calc_win.setter
    def calc_win(self, value):
        self._calc_win = value

    def get(self):

        st = self.stream

        # simulate Wood-Anderson response
        st.simulate(paz_remove=None, paz_simulate=self._paz)

        # trim waveform to common range
        stime, etime = common_range(st)
        st.trim(stime, etime)

        # get time delta from waveform data
        dt = st[0].stats.delta

        power = [np.power(tr.data, 2) for tr in st if tr.stats.channel[-1] not
                 in 'Z3']
        if len(power) != 2:
            raise ValueError('Wood-Anderson amplitude defintion only valid for '
                             'two horizontals: {0} given'.format(len(power)))
        power_sum = power[0] + power[1]
        #
        sqH = np.sqrt(power_sum)



    def calcwapp(self):
        print ("Getting Wood-Anderson peak-to-peak amplitude ...")
        print ("Simulating Wood-Anderson seismograph ...")

        self.wapp = None
        stream = self.stream

        stream.simulate(paz_remove=None, paz_simulate=paz_wa)

        trH1 = stream[0].data
        trH2 = stream[1].data
        ilen = min([len(trH1), len(trH2)])
        # get RMS of both horizontal components
        sqH = np.sqrt(np.power(trH1[0:ilen], 2) + np.power(trH2[0:ilen], 2))
        # get time array
        th = np.arange(0, len(sqH) * stream[0].stats.delta, stream[0].stats.delta)
        # get maximum peak within pick window
        iwin = getsignalwin(th, self.t0, self.pwin)
        self.wapp = np.max(sqH[iwin])
        print ("Determined Wood-Anderson peak-to-peak amplitude: %f mm") % self.wapp

        if self.plot_flag > 1:
            stream.plot()
            f = plt.figure(2)
            plt.plot(th, sqH)
            plt.plot(th[iwin], sqH[iwin], 'g')
            plt.plot([self.t0, self.t0], [0, max(sqH)], 'r', linewidth=2)
            plt.title('Station %s, RMS Horizontal Traces, WA-peak-to-peak=%4.1f mm' \
                      % (stream[0].stats.station, self.wapp))
            plt.xlabel('Time [s]')
            plt.ylabel('Displacement [mm]')
            plt.show()
            raw_input()
            plt.close(f)


class MomentMagnitude(Magnitude):
    '''
    Method to calculate seismic moment Mo and moment magnitude Mw.
    Requires results of class calcsourcespec for calculating plateau w0
    and corner frequency fc of source spectrum, respectively. Uses
    subfunction calcMoMw.py. Returns modified dictionary of picks including
    Dc-value, corner frequency fc, seismic moment Mo and
    corresponding moment magntiude Mw.
    '''

    def run_calcMoMw(self):

        picks = self.getpicks()
        nllocfile = self.getNLLocfile()
        wfdat = self.stream
        self.picdic = None

        for key in picks:
            if picks[key]['P']['weight'] < 4:
                # select waveform
                selwf = wfdat.select(station=key)
                if len(key) > 4:
                    Ppattern = '%s  ?    ?    ? P' % key
                elif len(key) == 4:
                    Ppattern = '%s   ?    ?    ? P' % key
                elif len(key) < 4:
                    Ppattern = '%s    ?    ?    ? P' % key
                nllocline = getPatternLine(nllocfile, Ppattern)
                # get hypocentral distance, station azimuth and
                # angle of incidence from NLLoc-location file
                delta = float(nllocline.split(None)[21])
                az = float(nllocline.split(None)[22])
                inc = float(nllocline.split(None)[24])
                # call subfunction to estimate source spectrum
                # and to derive w0 and fc
                [w0, fc] = calcsourcespec(selwf, picks[key]['P']['mpp'], \
                                          self.get_metadata(), self.getvp(), delta, az, \
                                          inc, self.getQp(), self.plot_flag)

                if w0 is not None:
                    # call subfunction to calculate Mo and Mw
                    zdat = selwf.select(component="Z")
                    if len(zdat) == 0:  # check for other components
                        zdat = selwf.select(component="3")
                    [Mo, Mw] = calcMoMw(zdat, w0, self.getrho(), self.getvp(),
                                        delta)
                else:
                    Mo = None
                    Mw = None

                # add w0, fc, Mo and Mw to dictionary
                picks[key]['P']['w0'] = w0
                picks[key]['P']['fc'] = fc
                picks[key]['P']['Mo'] = Mo
                picks[key]['P']['Mw'] = Mw
                self.picdic = picks


def calc_woodanderson_pp(st, metadata, T0, win=10., verbosity=False):
    if verbosity:
        print ("Getting Wood-Anderson peak-to-peak amplitude ...")
        print ("Simulating Wood-Anderson seismograph ...")

    # poles, zeros and sensitivity of WA seismograph
    # (see Uhrhammer & Collins, 1990, BSSA, pp. 702-716)
    paz_wa = {
        'poles': [5.6089 - 5.4978j, -5.6089 - 5.4978j],
        'zeros': [0j, 0j],
        'gain': 2080,
        'sensitivity': 1
    }

    stime, etime = common_range(st)
    st.trim(stime, etime)

    invtype, inventory = metadata

    st = restitute_data(st, invtype, inventory)

    dt = st[0].stats.delta

    st.simulate(paz_remove=None, paz_simulate=paz_wa)

    # get RMS of both horizontal components
    power = [np.power(tr.data, 2) for tr in st if tr.stats.channel[-1] not
              in 'Z3']
    if len(power) != 2:
        raise ValueError('Wood-Anderson amplitude defintion only valid for '
                         'two horizontals: {0} given'.format(len(power)))
    power_sum = power[0] + power[1]
    sqH = np.sqrt(power_sum)

    # get time array
    th = np.arange(0, len(sqH) * dt, dt)
    # get maximum peak within pick window
    iwin = getsignalwin(th, T0 - stime, win)
    wapp = np.max(sqH[iwin])
    if verbosity:
        print("Determined Wood-Anderson peak-to-peak amplitude: {0} "
              "mm".format(wapp))
    return wapp


def calcMoMw(wfstream, w0, rho, vp, delta):
    '''
    Subfunction of run_calcMoMw to calculate individual
    seismic moments and corresponding moment magnitudes.

    :param: wfstream
    :type:  `~obspy.core.stream.Stream`

    :param: w0, height of plateau of source spectrum
    :type:  float

    :param: rho, rock density [kg/m³]
    :type:  integer

    :param: delta, hypocentral distance [km]
    :type:  integer

    :param: inv, name/path of inventory or dataless-SEED file
    :type:  string
    '''

    tr = wfstream[0]
    delta = delta * 1000  # hypocentral distance in [m]

    print("calcMoMw: Calculating seismic moment Mo and moment magnitude Mw for station %s ..." \
          % tr.stats.station)

    # additional common parameters for calculating Mo
    rP = 2 / np.sqrt(15)  # average radiation pattern of P waves (Aki & Richards, 1980)
    freesurf = 2.0  # free surface correction, assuming vertical incidence

    Mo = w0 * 4 * np.pi * rho * np.power(vp, 3) * delta / (rP * freesurf)

    # Mw = np.log10(Mo * 1e07) * 2 / 3 - 10.7 # after Hanks & Kanamori (1979), defined for [dyn*cm]!
    Mw = np.log10(Mo) * 2 / 3 - 6.7  # for metric units

    print("calcMoMw: Calculated seismic moment Mo = %e Nm => Mw = %3.1f " % (Mo, Mw))

    return Mo, Mw


def calcsourcespec(wfstream, onset, metadata, vp, delta, azimuth, incidence,
                   qp, iplot):
    '''
    Subfunction to calculate the source spectrum and to derive from that the plateau
    (usually called omega0) and the corner frequency assuming Aki's omega-square
    source model. Has to be derived from instrument corrected displacement traces,
    thus restitution and integration necessary! Integrated traces are rotated
    into ray-coordinate system ZNE => LQT using Obspy's rotate modul!

    :param: wfstream
    :type:  `~obspy.core.stream.Stream`

    :param: onset, P-phase onset time
    :type: float

    :param: metadata, tuple or list containing type of inventory and either
    list of files or inventory object
    :type:  tuple or list

    :param: vp, Vp-wave velocity
    :type:  float

    :param: delta, hypocentral distance [km]
    :type:  integer

    :param: azimuth
    :type:  integer

    :param: incidence
    :type:  integer

    :param: Qp, quality factor for P-waves
    :type:  integer

    :param: iplot, show results (iplot>1) or not (iplot(<2)
    :type:  integer
    '''
    print ("Calculating source spectrum ....")

    # get Q value
    Q, A = qp

    delta = delta * 1000  # hypocentral distance in [m]

    fc = None
    w0 = None
    wf_copy = wfstream.copy()

    invtype, inventory = metadata

    [cordat, restflag] = restitute_data(wf_copy, invtype, inventory)
    if restflag is True:
        zdat = cordat.select(component="Z")
        if len(zdat) == 0:
            zdat = cordat.select(component="3")
        cordat_copy = cordat.copy()
        # get equal time stamps and lengths of traces
        # necessary for rotation of traces
        trstart, trend = common_range(cordat_copy)
        cordat_copy.trim(trstart, trend)
        try:
            # rotate into LQT (ray-coordindate-) system using Obspy's rotate
            # L: P-wave direction
            # Q: SV-wave direction
            # T: SH-wave direction
            LQT = cordat_copy.rotate('ZNE->LQT', azimuth, incidence)
            ldat = LQT.select(component="L")
            if len(ldat) == 0:
                # if horizontal channels are 2 and 3
                # no azimuth information is available and thus no
                # rotation is possible!
                print("calcsourcespec: Azimuth information is missing, "
                      "no rotation of components possible!")
                ldat = LQT.select(component="Z")

            # integrate to displacement
            # unrotated vertical component (for copmarison)
            inttrz = signal.detrend(integrate.cumtrapz(zdat[0].data, None,
                                                       zdat[0].stats.delta))
            # rotated component Z => L
            Ldat = signal.detrend(integrate.cumtrapz(ldat[0].data, None,
                                                     ldat[0].stats.delta))

            # get window after P pulse for
            # calculating source spectrum
            tstart = UTCDateTime(zdat[0].stats.starttime)
            tonset = onset.timestamp - tstart.timestamp
            impickP = tonset * zdat[0].stats.sampling_rate
            wfzc = Ldat[impickP: len(Ldat) - 1]
            # get time array
            t = np.arange(0, len(inttrz) * zdat[0].stats.delta, \
                          zdat[0].stats.delta)
            # calculate spectrum using only first cycles of
            # waveform after P onset!
            zc = crossings_nonzero_all(wfzc)
            if np.size(zc) == 0 or len(zc) <= 3:
                print ("calcsourcespec: Something is wrong with the waveform, "
                       "no zero crossings derived!")
                print ("No calculation of source spectrum possible!")
                plotflag = 0
            else:
                plotflag = 1
                index = min([3, len(zc) - 1])
                calcwin = (zc[index] - zc[0]) * zdat[0].stats.delta
                iwin = getsignalwin(t, tonset, calcwin)
                xdat = Ldat[iwin]

                # fft
                fny = zdat[0].stats.sampling_rate / 2
                l = len(xdat) / zdat[0].stats.sampling_rate
                # number of fft bins after Bath
                n = zdat[0].stats.sampling_rate * l
                # find next power of 2 of data length
                m = pow(2, np.ceil(np.log(len(xdat)) / np.log(2)))
                N = int(np.power(m, 2))
                y = zdat[0].stats.delta * np.fft.fft(xdat, N)
                Y = abs(y[: N / 2])
                L = (N - 1) / zdat[0].stats.sampling_rate
                f = np.arange(0, fny, 1 / L)

                # remove zero-frequency and frequencies above
                # corner frequency of seismometer (assumed
                # to be 100 Hz)
                fi = np.where((f >= 1) & (f < 100))
                F = f[fi]
                YY = Y[fi]

                # correction for attenuation
                wa = 2 * np.pi * F  # angular frequency
                D = np.exp((wa * delta) / (2 * vp * Q * F ** A))
                YYcor = YY.real * D

                # get plateau (DC value) and corner frequency
                # initial guess of plateau
                w0in = np.mean(YYcor[0:100])
                # initial guess of corner frequency
                # where spectral level reached 50% of flat level
                iin = np.where(YYcor >= 0.5 * w0in)
                Fcin = F[iin[0][np.size(iin) - 1]]

                # use of implicit scipy otimization function
                fit = synthsourcespec(F, w0in, Fcin)
                [optspecfit, _] = curve_fit(synthsourcespec, F, YYcor, [w0in,
                                                                        Fcin])
                w01 = optspecfit[0]
                fc1 = optspecfit[1]
                print ("calcsourcespec: Determined w0-value: %e m/Hz, \n"
                       "Determined corner frequency: %f Hz" % (w01, fc1))

                # use of conventional fitting
                [w02, fc2] = fitSourceModel(F, YYcor, Fcin, iplot)

                # get w0 and fc as median of both
                # source spectrum fits
                w0 = np.median([w01, w02])
                fc = np.median([fc1, fc2])
                print("calcsourcespec: Using w0-value = %e m/Hz and fc = %f Hz" % (w0, fc))

        except TypeError as er:
            raise TypeError('''{0}'''.format(er))

        if iplot > 1:
            f1 = plt.figure()
            tLdat = np.arange(0, len(Ldat) * zdat[0].stats.delta, \
                              zdat[0].stats.delta)
            plt.subplot(2, 1, 1)
            # show displacement in mm
            p1, = plt.plot(t, np.multiply(inttrz, 1000), 'k')
            p2, = plt.plot(tLdat, np.multiply(Ldat, 1000))
            plt.legend([p1, p2], ['Displacement', 'Rotated Displacement'])
            if plotflag == 1:
                plt.plot(t[iwin], np.multiply(xdat, 1000), 'g')
                plt.title('Seismogram and P Pulse, Station %s-%s' \
                          % (zdat[0].stats.station, zdat[0].stats.channel))
            else:
                plt.title('Seismogram, Station %s-%s' \
                          % (zdat[0].stats.station, zdat[0].stats.channel))
            plt.xlabel('Time since %s' % zdat[0].stats.starttime)
            plt.ylabel('Displacement [mm]')

            if plotflag == 1:
                plt.subplot(2, 1, 2)
                p1, = plt.loglog(f, Y.real, 'k')
                p2, = plt.loglog(F, YY.real)
                p3, = plt.loglog(F, YYcor, 'r')
                p4, = plt.loglog(F, fit, 'g')
                plt.loglog([fc, fc], [w0 / 100, w0], 'g')
                plt.legend([p1, p2, p3, p4], ['Raw Spectrum', \
                                              'Used Raw Spectrum', \
                                              'Q-Corrected Spectrum', \
                                              'Fit to Spectrum'])
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

    # ssp = omega0 / (pow(2, (1 + f / fcorner)))
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

    w0 = []
    stdw0 = []
    fc = []
    stdfc = []
    STD = []

    # get window around initial corner frequency for trials
    fcstopl = fc0 - max(1, len(f) / 10)
    il = np.argmin(abs(f - fcstopl))
    fcstopl = f[il]
    fcstopr = fc0 + min(len(f), len(f) / 10)
    ir = np.argmin(abs(f - fcstopr))
    fcstopr = f[ir]
    iF = np.where((f >= fcstopl) & (f <= fcstopr))

    # vary corner frequency around initial point
    for i in range(il, ir):
        FC = f[i]
        indexdc = np.where((f > 0) & (f <= FC))
        dc = np.mean(S[indexdc])
        stddc = np.std(dc - S[indexdc])
        w0.append(dc)
        stdw0.append(stddc)
        fc.append(FC)
        # slope
        indexfc = np.where((f >= FC) & (f <= fcstopr))
        yi = dc / (1 + (f[indexfc] / FC) ** 2)
        stdFC = np.std(yi - S[indexfc])
        stdfc.append(stdFC)
        STD.append(stddc + stdFC)

    # get best found w0 anf fc from minimum
    if len(STD) > 0:
        fc = fc[np.argmin(STD)]
        w0 = w0[np.argmin(STD)]
    elif len(STD) == 0:
        fc = fc0
        w0 = max(S)

    print("fitSourceModel: best fc: %fHz, best w0: %e m/Hz" \
          % (fc, w0))

    if iplot > 1:
        plt.figure(iplot)
        plt.loglog(f, S, 'k')
        plt.loglog([f[0], fc], [w0, w0], 'g')
        plt.loglog([fc, fc], [w0 / 100, w0], 'g')
        plt.title('Calculated Source Spectrum, Omega0=%e m/Hz, fc=%6.2f Hz' \
                  % (w0, fc))
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [m/Hz]')
        plt.grid()
        plt.figure(iplot + 1)
        plt.subplot(311)
        plt.plot(f[il:ir], STD, '*')
        plt.title('Common Standard Deviations')
        plt.xticks([])
        plt.subplot(312)
        plt.plot(f[il:ir], stdw0, '*')
        plt.title('Standard Deviations of w0-Values')
        plt.xticks([])
        plt.subplot(313)
        plt.plot(f[il:ir], stdfc, '*')
        plt.title('Standard Deviations of Corner Frequencies')
        plt.xlabel('Corner Frequencies [Hz]')
        plt.show()
        raw_input()
        plt.close()

    return w0, fc


def calc_richter_magnitude(e, wf, metadata, amp_win):
    if e.origins:
        o = e.origins[0]
        mags = dict()
        for a in o.arrivals:
            if a.phase not in 'sS':
                continue
            pick = a.pick_id.get_referred_object()
            station = pick.waveform_id.station_code
            wf = select_for_phase(wf.select(
                station=station), a.phase)
            delta = degrees2kilometers(a.distance)
            wapp = calc_woodanderson_pp(wf, metadata, pick.time,
                                        amp_win)
            # using standard Gutenberg-Richter relation
            # TODO make the ML calculation more flexible by allowing
            # use of custom relation functions
            mag = np.log10(wapp) + richter_magnitude_scaling(delta)
            mags[station] = mag
        mag = np.median([M for M in mags.values()])
        # give some information on the processing
        print('number of stations used: {0}\n'.format(len(mags.values())))
        print('stations used:\n')
        for s in mags.keys(): print('\t{0}'.format(s))

        return ObspyMagnitude(mag=mag, magnitude_type='ML')
    else:
        return None


def calc_moment_magnitude(e, wf, metadata, vp, Qp, rho):
    if e.origins:
        o = e.origins[0]
        mags = dict()
        for a in o.arrivals:
            if a.phase in 'sS':
                continue
            pick = a.pick_id.get_referred_object()
            station = pick.waveform_id.station_code
            wf = wf.select(station=station)
            if not wf:
                continue
            onset = pick.time
            dist = degrees2kilometers(a.distance)
            w0, fc = calcsourcespec(wf, onset, metadata, vp, dist, a.azimuth, a.takeoff_angle, Qp, 0)
            if w0 is None or fc is None:
                continue
            station_mag = calcMoMw(wf, w0, rho, vp, dist)
            mags[station] = station_mag
        mag = np.median([M[1] for M in mags.values()])
        # give some information on the processing
        print('number of stations used: {0}\n'.format(len(mags.values())))
        print('stations used:\n')
        for s in mags.keys(): print('\t{0}'.format(s))
        return ObspyMagnitude(mag=mag, magnitude_type='Mw')
    else:
        return None
