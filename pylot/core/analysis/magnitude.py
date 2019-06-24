#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created autumn/winter 2015.
Revised/extended summer 2017.

:author: Ludger Küperkoch / MAGS2 EP3 working group
"""

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
    distance = np.array([0, 10, 20, 25, 30, 35, 40, 45, 50, 60, 70, 75, 85, 90, 100, 110,
                         120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 230, 240, 250,
                         260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380,
                         390, 400, 430, 470, 510, 560, 600, 700, 800, 900, 1000])
    richter_scaling = np.array([1.4, 1.5, 1.7, 1.9, 2.1, 2.3, 2.4, 2.5, 2.6, 2.8, 2.8, 2.9,
                                2.9, 3.0, 3.1, 3.1, 3.2, 3.2, 3.3, 3.3, 3.4, 3.4, 3.5, 3.5,
                                3.6, 3.7, 3.7, 3.8, 3.8, 3.9, 3.9, 4.0, 4.0, 4.1, 4.2, 4.2,
                                4.2, 4.2, 4.3, 4.3, 4.3, 4.4, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,
                                5.1, 5.2, 5.4, 5.5, 5.7])
    # prepare spline interpolation to calculate return value
    func, params = fit_curve(distance, richter_scaling)
    return func(delta, params)


class Magnitude(object):
    """
    Base class object for Magnitude calculation within PyLoT.
    """

    def __init__(self, stream, event, verbosity=False, iplot=0):
        self._type = "M"
        self._stream = stream
        self._plot_flag = iplot
        self._verbosity = verbosity
        self._event = event
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

    def updated_event(self, magscaling=None):
        net_ml = self.net_magnitude(magscaling)
        if net_ml:
            self.event.magnitudes.append(net_ml)
        return self.event

    def net_magnitude(self, magscaling=None):
        if self:
            if magscaling is None:
                scaling = False
            elif magscaling[0] != 0 and magscaling[1] != 0:
                scaling = False
            else:
                scaling = True
            if scaling:
                # scaling necessary
                print("Scaling network magnitude ...")
                mag = ope.Magnitude(
                    mag=np.median([M.mag for M in self.magnitudes.values()]) * \
                        magscaling[0] + magscaling[1],
                    magnitude_type=self.type,
                    origin_id=self.origin_id,
                    station_count=len(self.magnitudes),
                    azimuthal_gap=self.origin_id.get_referred_object().quality.azimuthal_gap)
            else:
                # no scaling necessary
                mag = ope.Magnitude(
                    mag=np.median([M.mag for M in self.magnitudes.values()]),
                    magnitude_type=self.type,
                    origin_id=self.origin_id,
                    station_count=len(self.magnitudes),
                    azimuthal_gap=self.origin_id.get_referred_object().quality.azimuthal_gap)
            return mag


class LocalMagnitude(Magnitude):
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

    def __init__(self, stream, event, calc_win, wascaling, verbosity=False, iplot=0):
        super(LocalMagnitude, self).__init__(stream, event, verbosity, iplot)

        self._calc_win = calc_win
        self._wascaling = wascaling
        self._type = 'ML'
        self.calc()

    @property
    def calc_win(self):
        return self._calc_win

    @calc_win.setter
    def calc_win(self, value):
        self._calc_win = value

    @property
    def wascaling(self):
        return self._wascaling

    @property
    def amplitudes(self):
        return self._amplitudes

    @amplitudes.setter
    def amplitudes(self, value):
        station, a0 = value
        self._amplitudes[station] = a0

    def peak_to_peak(self, st, t0):

        try:
            iplot = int(self.plot_flag)
        except:
            if self.plot_flag == True or self.plot_flag == 'True':
                iplot = 2
            else:
                iplot = 0

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
        th = np.arange(0, st[0].stats.npts / st[0].stats.sampling_rate, st[0].stats.delta)
        # get maximum peak within pick window
        iwin = getsignalwin(th, t0 - stime, self.calc_win)
        ii = min([iwin[len(iwin) - 1], len(th)])
        iwin = iwin[0:ii]
        wapp = np.max(sqH[iwin])
        if self.verbose:
            print("Determined Wood-Anderson peak-to-peak amplitude for station {0}: {1} "
                  "mm".format(st[0].stats.station, wapp))

        # check for plot flag (for debugging only)
        fig = None
        if iplot > 1:
            fig = plt.figure()
            ax = fig.add_subplot(211)
            ax.plot(th, st[0].data, 'k')
            ax.plot(th, sqH)
            ax.plot(th[iwin], sqH[iwin], 'g')
            ax.plot([t0 - stime, t0 - stime], [0, max(sqH)], 'r', linewidth=2)
            ax.set_title('Station %s, Channel %s, RMS Horizontal Trace, '
                         'WA-peak-to-peak=%6.3f mm' % (st[0].stats.station,
                                                       st[0].stats.channel,
                                                       wapp))
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Displacement [mm]')
            ax = fig.add_subplot(212)
            ax.plot(th, st[1].data, 'k')
            ax.plot(th, sqH)
            ax.plot(th[iwin], sqH[iwin], 'g')
            ax.plot([t0 - stime, t0 - stime], [0, max(sqH)], 'r', linewidth=2)
            ax.set_title('Channel %s, RMS Horizontal Trace, '
                         'WA-peak-to-peak=%6.3f mm' % (st[1].stats.channel,
                                                       wapp))
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Displacement [mm]')
            fig.show()
            try:
                input()
            except SyntaxError:
                pass
            plt.close(fig)

        return wapp, fig

    def calc(self):
        for a in self.arrivals:
            if a.phase not in 'sS':
                continue
            # make sure calculating Ml only from reliable onsets
            # NLLoc: time_weight = 0 => do not use onset!
            if a.time_weight == 0:
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
            a0, self.p2p_fig = self.peak_to_peak(wf, onset)
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
            # or scale WA amplitude with given scaling relation
            if str(self.wascaling) == '[0.0, 0.0, 0.0]':
                print("Calculating original Richter magnitude ...")
                magnitude = ope.StationMagnitude(mag=np.log10(a0) \
                                                     + richter_magnitude_scaling(delta))
            else:
                print("Calculating scaled local magnitude ...")
                a0 = a0 * 1e03  # mm to nm (see Havskov & Ottemöller, 2010)
                magnitude = ope.StationMagnitude(mag=np.log10(a0) \
                                                     + self.wascaling[0] * np.log10(delta) + self.wascaling[1]
                                                     * delta + self.wascaling[
                                                         2])
            if self.verbose:
                print(
                    "Local Magnitude for station {0}: ML = {1:3.1f}".format(
                       station, magnitude.mag))
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
            # make sure calculating Mo only from reliable onsets
            # NLLoc: time_weight = 0 => do not use onset!
            if a.time_weight == 0:
                continue
            pick = a.pick_id.get_referred_object()
            station = pick.waveform_id.station_code
            if len(self.stream) <= 2:
                print("Station:" '{0}'.format(station))
                print("WARNING: No instrument corrected data available,"
                      " no magnitude calculation possible! Go on.")
                continue
            scopy = self.stream.copy()
            wf = scopy.select(station=station)
            if not wf:
                continue
            onset = pick.time
            distance = degrees2kilometers(a.distance)
            azimuth = a.azimuth
            incidence = a.takeoff_angle
            w0, fc = calcsourcespec(wf, onset, self.p_velocity, distance,
                                    azimuth, incidence, self.p_attenuation,
                                    self.plot_flag, self.verbose)
            if w0 is None or fc is None:
                if self.verbose:
                    print("WARNING: insufficient frequency information")
                continue
            WF = select_for_phase(self.stream.select(
                station=station), a.phase)
            WF = select_for_phase(WF, "P")
            m0, mw = calcMoMw(WF, w0, self.rock_density, self.p_velocity,
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
        print("Calculating source spectrum for station %s ...." % wfstream[0].stats.station)

    try:
        iplot = int(iplot)
    except:
        if iplot == True or iplot == 'True':
            iplot = 2
        else:
            iplot = 0

    # get Q value
    Q, A = qp

    dist = delta * 1000  # hypocentral distance in [m]

    fc = None
    w0 = None

    zdat = select_for_phase(wfstream, "P")

    if len(zdat) == 0:
        raise IOError('No vertical component found in stream:\n{}'.format(wfstream))

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
        # if horizontal channels are 1 and 2
        # no azimuth information is available and thus no
        # rotation is possible!
        if verbosity:
            print("calcsourcespec: Azimuth information is missing, "
                  "no rotation of components possible!")
        # instead, use component 3
        ldat = LQT.select(component="3")
        if len(ldat) == 0:
            # maybe component z available
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
            print("calcsourcespec: Something is wrong with the waveform, "
                  "no zero crossings derived!\n")
            print("No calculation of source spectrum possible!")
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
        w0 = optspecfit[0]
        fc = optspecfit[1]
        # w01 = optspecfit[0]
        # fc1 = optspecfit[1]
        if verbosity:
            print("calcsourcespec: Determined w0-value: %e m/Hz, \n"
                  "calcsourcespec: Determined corner frequency: %f Hz" % (w0, fc))

        # use of conventional fitting
        # [w02, fc2] = fitSourceModel(F, YYcor, Fcin, iplot, verbosity)

        # get w0 and fc as median of both
        # source spectrum fits
        # w0 = np.median([w01, w02])
        # fc = np.median([fc1, fc2])
        # if verbosity:
        #    print("calcsourcespec: Using w0-value = %e m/Hz and fc = %f Hz" % (
        #    w0, fc))
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
            plt.legend([p1, p2, p3, p4], ['Raw Spectrum',
                                          'Used Raw Spectrum',
                                          'Q-Corrected Spectrum',
                                          'Fit to Spectrum'])
            plt.title('Source Spectrum from P Pulse, w0=%e m/Hz, fc=%6.2f Hz' \
                      % (w0, fc))
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Amplitude [m/Hz]')
            plt.grid()
            plt.show()
            try:
                input()
            except SyntaxError:
                pass
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

    try:
        iplot = int(iplot)
    except:
        if iplot == True or iplot == 'True':
            iplot = 2
        else:
            iplot = 0

    w0 = []
    stdw0 = []
    fc = []
    stdfc = []
    STD = []
    # get window around initial corner frequency for trials
    # left side of initial corner frequency
    fcstopl = max(f[0], fc0 - max(1, fc0 / 2))
    il = np.where(f <= fcstopl)
    il = il[0][np.size(il) - 1]
    # right side of initial corner frequency
    fcstopr = min(fc0 + (fc0 / 2), f[len(f) - 1])
    ir = np.where(f >= fcstopr)
    # check, if fcstopr is available
    if np.size(ir) == 0:
        fcstopr = fc0
        ir = len(f) - 1
    else:
        ir = ir[0][0]

    # vary corner frequency around initial point
    print("fitSourceModel: Varying corner frequency "
          "around initial corner frequency ...")
    # check difference of il and ir in order to
    # keep calculation time acceptable
    idiff = ir - il
    if idiff > 10000:
        increment = 100
    elif idiff <= 20:
        increment = 1
    else:
        increment = 10

    for i in range(il, ir, increment):
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
        plt.figure()  # iplot)
        plt.loglog(f, S, 'k')
        plt.loglog([f[0], fc], [w0, w0], 'g')
        plt.loglog([fc, fc], [w0 / 100, w0], 'g')
        plt.title('Calculated Source Spectrum, Omega0=%e m/Hz, fc=%6.2f Hz' \
                  % (w0, fc))
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [m/Hz]')
        plt.grid()
        plt.figure()  # iplot + 1)
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
        try:
            input()
        except SyntaxError:
            pass
        plt.close()

    return w0, fc
