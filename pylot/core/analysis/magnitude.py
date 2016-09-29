#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created autumn/winter 2015.

:author: Ludger Küperkoch / MAGS2 EP3 working group
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import obspy.core.event as ope
from obspy.geodetics import degrees2kilometers
from scipy import integrate, signal
from scipy.optimize import curve_fit

from pylot.core.pick.utils import getsignalwin, crossings_nonzero_all, \
    select_for_phase
from pylot.core.util.utils import common_range, fit_curve


def richter_magnitude_scaling(delta):
    relation = np.loadtxt(os.path.join(os.path.expanduser('~'),
                                       '.pylot', 'richter_scaling.data'))
    # prepare spline interpolation to calculate return value
    func, params = fit_curve(relation[:, 0], relation[:, 1])
    return func(delta, params)


class Magnitude(object):
    """
    Base class object for Magnitude calculation within PyLoT.
    """

    def __init__(self, stream, event, verbosity=False, iplot=0):
        self._type = "M"
        self._plot_flag = iplot
        self._verbosity = verbosity
        self._event = event
        self._stream = stream
        self._magnitudes = dict()

    def __str__(self):
        print(
        'number of stations used: {0}\n'.format(len(self.magnitudes.values())))
        print('\tstation\tmagnitude')
        for s, m in self.magnitudes.items(): print('\t{0}\t{1}'.format(s, m))

    def __nonzero__(self):
        return bool(self.magnitudes)

    @property
    def type(self):
        return self._type

    @property
    def plot_flag(self):
        return self._plot_flag

    @plot_flag.setter
    def plot_flag(self, value):
        self._plot_flag = value

    @property
    def verbose(self):
        return self._verbosity

    @verbose.setter
    def verbose(self, value):
        if not isinstance(value, bool):
            print('WARNING: only boolean values accepted...\n')
            value = bool(value)
        self._verbosity = value

    @property
    def stream(self):
        return self._stream

    @stream.setter
    def stream(self, value):
        self._stream = value

    @property
    def event(self):
        return self._event

    @property
    def origin_id(self):
        return self._event.origins[0].resource_id

    @property
    def arrivals(self):
        return self._event.origins[0].arrivals

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

    def calc(self):
        pass

    def updated_event(self):
        self.event.magnitudes.append(self.net_magnitude())
        return self.event

    def net_magnitude(self):
        if self:
            # TODO if an average Magnitude instead of the median is calculated
            # StationMagnitudeContributions should be added to the returned
            # Magnitude object
            # mag_error => weights (magnitude error estimate from peak_to_peak, calcsourcespec?)
            # weights => StationMagnitdeContribution
            mag = ope.Magnitude(
                mag=np.median([M.mag for M in self.magnitudes.values()]),
                magnitude_type=self.type,
                origin_id=self.origin_id,
                station_count=len(self.magnitudes),
                azimuthal_gap=self.origin_id.get_referred_object().quality.azimuthal_gap)
            return mag
        return None


class RichterMagnitude(Magnitude):
    """
    Method to derive peak-to-peak amplitude as seen on a Wood-Anderson-
    seismograph. Has to be derived from instrument corrected traces!
    """

    # poles, zeros and sensitivity of WA seismograph
    # (see Uhrhammer & Collins, 1990, BSSA, pp. 702-716)
    _paz = {
        'poles': [5.6089 - 5.4978j, -5.6089 - 5.4978j],
        'zeros': [0j, 0j],
        'gain': 2080,
        'sensitivity': 1
    }

    _amplitudes = dict()

    def __init__(self, stream, event, calc_win, verbosity=False, iplot=0):
        super(RichterMagnitude, self).__init__(stream, event, verbosity, iplot)

        self._calc_win = calc_win
        self._type = 'ML'
        self.calc()

    @property
    def calc_win(self):
        return self._calc_win

    @calc_win.setter
    def calc_win(self, value):
        self._calc_win = value

    @property
    def amplitudes(self):
        return self._amplitudes

    @amplitudes.setter
    def amplitudes(self, value):
        station, a0 = value
        self._amplitudes[station] = a0

    def peak_to_peak(self, st, t0):

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

        # get time array
        th = np.arange(0, len(sqH) * dt, dt)
        # get maximum peak within pick window
        iwin = getsignalwin(th, t0 - stime, self.calc_win)
        wapp = np.max(sqH[iwin])
        if self.verbose:
            print("Determined Wood-Anderson peak-to-peak amplitude: {0} "
                  "mm".format(wapp))

        # check for plot flag (for debugging only)
        if self.plot_flag > 1:
            st.plot()
            f = plt.figure(2)
            plt.plot(th, sqH)
            plt.plot(th[iwin], sqH[iwin], 'g')
            plt.plot([t0, t0], [0, max(sqH)], 'r', linewidth=2)
            plt.title(
                'Station %s, RMS Horizontal Traces, WA-peak-to-peak=%4.1f mm' \
                % (st[0].stats.station, wapp))
            plt.xlabel('Time [s]')
            plt.ylabel('Displacement [mm]')
            plt.show()
            raw_input()
            plt.close(f)

        return wapp

    def calc(self):
        for a in self.arrivals:
            if a.phase not in 'sS':
                continue
            pick = a.pick_id.get_referred_object()
            station = pick.waveform_id.station_code
            wf = select_for_phase(self.stream.select(
                station=station), a.phase)
            if not wf:
                if self.verbose:
                    print(
                    'WARNING: no waveform data found for station {0}'.format(
                        station))
                continue
            delta = degrees2kilometers(a.distance)
            onset = pick.time
            a0 = self.peak_to_peak(wf, onset)
            amplitude = ope.Amplitude(generic_amplitude=a0 * 1e-3)
            amplitude.unit = 'm'
            amplitude.category = 'point'
            amplitude.waveform_id = pick.waveform_id
            amplitude.magnitude_hint = self.type
            amplitude.pick_id = pick.resource_id
            amplitude.type = 'AML'
            self.event.amplitudes.append(amplitude)
            self.amplitudes = (station, amplitude)
            # using standard Gutenberg-Richter relation
            # TODO make the ML calculation more flexible by allowing
            # use of custom relation functions
            magnitude = ope.StationMagnitude(
                mag=np.log10(a0) + richter_magnitude_scaling(delta))
            magnitude.origin_id = self.origin_id
            magnitude.waveform_id = pick.waveform_id
            magnitude.amplitude_id = amplitude.resource_id
            magnitude.station_magnitude_type = self.type
            self.event.station_magnitudes.append(magnitude)
            self.magnitudes = (station, magnitude)


class MomentMagnitude(Magnitude):
    '''
    Method to calculate seismic moment Mo and moment magnitude Mw.
    Requires results of class calcsourcespec for calculating plateau w0
    and corner frequency fc of source spectrum, respectively. Uses
    subfunction calcMoMw.py. Returns modified dictionary of picks including
    Dc-value, corner frequency fc, seismic moment Mo and
    corresponding moment magntiude Mw.
    '''

    _props = dict()

    def __init__(self, stream, event, vp, Qp, density, verbosity=False,
                 iplot=False):
        super(MomentMagnitude, self).__init__(stream, event, verbosity, iplot)

        self._vp = vp
        self._Qp = Qp
        self._density = density
        self._type = 'Mw'
        self.calc()

    @property
    def p_velocity(self):
        return self._vp

    @property
    def p_attenuation(self):
        return self._Qp

    @property
    def rock_density(self):
        return self._density

    @property
    def moment_props(self):
        return self._props

    @moment_props.setter
    def moment_props(self, value):
        station, props = value
        self._props[station] = props

    @property
    def seismic_moment(self):
        return self._m0

    @seismic_moment.setter
    def seismic_moment(self, value):
        self._m0 = value

    def calc(self):
        for a in self.arrivals:
            if a.phase not in 'pP':
                continue
            pick = a.pick_id.get_referred_object()
            station = pick.waveform_id.station_code
            wf = select_for_phase(self.stream.select(
                station=station), a.phase)
            if not wf:
                continue
            onset = pick.time
            distance = degrees2kilometers(a.distance)
            azimuth = a.azimuth
            incidence = a.takeoff_angle
            w0, fc = calcsourcespec(wf, onset, self.p_velocity, distance,
                                    azimuth,
                                    incidence, self.p_attenuation,
                                    self.plot_flag, self.verbose)
            if w0 is None or fc is None:
                if self.verbose:
                    print("WARNING: insufficient frequency information")
                continue
            wf = select_for_phase(wf, "P")
            m0, mw = calcMoMw(wf, w0, self.rock_density, self.p_velocity,
                              distance, self.verbose)
            self.moment_props = (station, dict(w0=w0, fc=fc, Mo=m0))
            magnitude = ope.StationMagnitude(mag=mw)
            magnitude.origin_id = self.origin_id
            magnitude.waveform_id = pick.waveform_id
            magnitude.station_magnitude_type = self.type
            self.event.station_magnitudes.append(magnitude)
            self.magnitudes = (station, magnitude)


def calcMoMw(wfstream, w0, rho, vp, delta, verbosity=False):
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

    if verbosity:
        print(
        "calcMoMw: Calculating seismic moment Mo and moment magnitude Mw for station {0} ...".format(
            tr.stats.station))

    # additional common parameters for calculating Mo
    rP = 2 / np.sqrt(
        15)  # average radiation pattern of P waves (Aki & Richards, 1980)
    freesurf = 2.0  # free surface correction, assuming vertical incidence

    Mo = w0 * 4 * np.pi * rho * np.power(vp, 3) * delta / (rP * freesurf)

    # Mw = np.log10(Mo * 1e07) * 2 / 3 - 10.7 # after Hanks & Kanamori (1979), defined for [dyn*cm]!
    Mw = np.log10(Mo) * 2 / 3 - 6.7  # for metric units

    if verbosity:
        print(
        "calcMoMw: Calculated seismic moment Mo = {0} Nm => Mw = {1:3.1f} ".format(
            Mo, Mw))

    return Mo, Mw


def calcsourcespec(wfstream, onset, vp, delta, azimuth, incidence,
                   qp, iplot=0, verbosity=False):
    '''
    Subfunction to calculate the source spectrum and to derive from that the plateau
    (usually called omega0) and the corner frequency assuming Aki's omega-square
    source model. Has to be derived from instrument corrected displacement traces,
    thus restitution and integration necessary! Integrated traces are rotated
    into ray-coordinate system ZNE => LQT using Obspy's rotate modul!

    :param: wfstream (corrected for instrument)
    :type:  `~obspy.core.stream.Stream`

    :param: onset, P-phase onset time
    :type: float

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

    :param: iplot, show results (iplot>1) or not (iplot<1)
    :type:  integer
    '''
    if verbosity:
        print ("Calculating source spectrum ....")

    # get Q value
    Q, A = qp

    dist = delta * 1000  # hypocentral distance in [m]

    fc = None
    w0 = None

    zdat = select_for_phase(wfstream, "P")

    dt = zdat[0].stats.delta

    freq = zdat[0].stats.sampling_rate

    # trim traces to common range (for rotation)
    trstart, trend = common_range(wfstream)
    wfstream.trim(trstart, trend)

    # rotate into LQT (ray-coordindate-) system using Obspy's rotate
    # L: P-wave direction
    # Q: SV-wave direction
    # T: SH-wave direction
    LQT = wfstream.rotate('ZNE->LQT', azimuth, incidence)
    ldat = LQT.select(component="L")
    if len(ldat) == 0:
        # if horizontal channels are 2 and 3
        # no azimuth information is available and thus no
        # rotation is possible!
        if verbosity:
            print("calcsourcespec: Azimuth information is missing, "
                  "no rotation of components possible!")
        ldat = LQT.select(component="Z")

    # integrate to displacement
    # unrotated vertical component (for comparison)
    inttrz = signal.detrend(integrate.cumtrapz(zdat[0].data, None, dt))

    # rotated component Z => L
    Ldat = signal.detrend(integrate.cumtrapz(ldat[0].data, None, dt))

    # get window after P pulse for
    # calculating source spectrum
    rel_onset = onset - trstart
    impickP = int(rel_onset * freq)
    wfzc = Ldat[impickP: len(Ldat) - 1]
    # get time array
    t = np.arange(0, len(inttrz) * dt, dt)
    # calculate spectrum using only first cycles of
    # waveform after P onset!
    zc = crossings_nonzero_all(wfzc)
    if np.size(zc) == 0 or len(zc) <= 3:
        if verbosity:
            print ("calcsourcespec: Something is wrong with the waveform, "
                   "no zero crossings derived!\n")
            print ("No calculation of source spectrum possible!")
        plotflag = 0
    else:
        plotflag = 1
        index = min([3, len(zc) - 1])
        calcwin = (zc[index] - zc[0]) * dt
        iwin = getsignalwin(t, rel_onset, calcwin)
        xdat = Ldat[iwin]

        # fft
        fny = freq / 2
        l = len(xdat) / freq
        # number of fft bins after Bath
        n = freq * l
        # find next power of 2 of data length
        m = pow(2, np.ceil(np.log(len(xdat)) / np.log(2)))
        N = int(np.power(m, 2))
        y = dt * np.fft.fft(xdat, N)
        Y = abs(y[: N / 2])
        L = (N - 1) / freq
        f = np.arange(0, fny, 1 / L)

        # remove zero-frequency and frequencies above
        # corner frequency of seismometer (assumed
        # to be 100 Hz)
        fi = np.where((f >= 1) & (f < 100))
        F = f[fi]
        YY = Y[fi]

        # correction for attenuation
        wa = 2 * np.pi * F  # angular frequency
        D = np.exp((wa * dist) / (2 * vp * Q * F ** A))
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
        [optspecfit, _] = curve_fit(synthsourcespec, F, YYcor, [w0in, Fcin])
        w01 = optspecfit[0]
        fc1 = optspecfit[1]
        if verbosity:
            print ("calcsourcespec: Determined w0-value: %e m/Hz, \n"
                   "Determined corner frequency: %f Hz" % (w01, fc1))

        # use of conventional fitting
        [w02, fc2] = fitSourceModel(F, YYcor, Fcin, iplot, verbosity)

        # get w0 and fc as median of both
        # source spectrum fits
        w0 = np.median([w01, w02])
        fc = np.median([fc1, fc2])
        if verbosity:
            print("calcsourcespec: Using w0-value = %e m/Hz and fc = %f Hz" % (
            w0, fc))

    if iplot > 1:
        f1 = plt.figure()
        tLdat = np.arange(0, len(Ldat) * dt, dt)
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


def fitSourceModel(f, S, fc0, iplot, verbosity=False):
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
    if verbosity:
        print(
        "fitSourceModel: best fc: {0} Hz, best w0: {1} m/Hz".format(fc, w0))

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
