#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np
from obspy import UTCDateTime
from pylot.core.util.utils import fit_curve, find_nearest, clims
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()
__author__ = 'sebastianw'

def create_axis(x0, incr, npts):
    ax = np.zeros(npts)
    for i in range(npts):
        ax[i] = x0 + incr * i
    return ax

def gauss_parameter(te, tm, tl, eta):
    '''
    takes three onset times and returns the parameters sig1, sig2, a1 and a2
    to represent the pick as a probability density funtion (PDF) with two
    Gauss branches
    :param te:
    :param tm:
    :param tl:
    :param eta:
    :return:
    '''

    sig1 = (tm - te) / np.sqrt(2 * np.log(1 / eta))
    sig2 = (tl - tm) / np.sqrt(2 * np.log(1 / eta))

    a1 = 2 / (1 + sig2 / sig1)
    a2 = 2 / (1 + sig1 / sig2)

    return tm, sig1, sig2, a1, a2


def exp_parameter(te, tm, tl, eta):
    '''
    takes three onset times te, tm and tl and returns the parameters sig1,
    sig2 and a to represent the pick as a probability density function (PDF)
    with two exponential decay branches
    :param te:
    :param tm:
    :param tl:
    :param eta:
    :return:
    '''

    sig1 = np.log(eta) / (te - tm)
    sig2 = np.log(eta) / (tm - tl)
    if np.isinf(sig1) == True:
        sig1 = np.log(eta) / (tm - tl)
    if np.isinf(sig2) == True:
        sig2 = np.log(eta) / (te - tm)
    a = 1 / (1 / sig1 + 1 / sig2)
    return tm, sig1, sig2, a


def gauss_branches(k, param_tuple):
    '''
    function gauss_branches takes an axes x, a center value mu, two sigma
    values sig1 and sig2 and two scaling factors a1 and a2 and return a
    list containing the values of a probability density function (PDF)
    consisting of gauss branches
    :param x:
    :type x:
    :param mu:
    :type mu:
    :param sig1:
    :type sig1:
    :param sig2:
    :type sig2:
    :param a1:
    :type a1:
    :param a2:
    :returns fun_vals: list with function values along axes x
    '''

    #python 3 workaround
    mu, sig1, sig2, a1, a2 = param_tuple
    
    def _func(k, mu, sig1, sig2, a1, a2):
        if k < mu:
            rval = a1 * 1 / (np.sqrt(2 * np.pi) * sig1) * np.exp(-((k - mu) / sig1) ** 2 / 2)
        else:
            rval = a2 * 1 / (np.sqrt(2 * np.pi) * sig2) * np.exp(-((k - mu) /
                                                                   sig2) ** 2 / 2)
        return rval

    try:
        return [_func(x, mu, sig1, sig2, a1, a2) for x in iter(k)]
    except TypeError:
        return _func(k, mu, sig1, sig2, a1, a2)


def exp_branches(k, param_tuple):
    '''
    function exp_branches takes an axes x, a center value mu, two sigma
    values sig1 and sig2 and a scaling factor a and return a
    list containing the values of a probability density function (PDF)
    consisting of exponential decay branches
    :param x:
    :param mu:
    :param sig1:
    :param sig2:
    :param a:
    :returns fun_vals: list with function values along axes x:
    '''

    #python 3 workaround
    mu, sig1, sig2, a = param_tuple
    
    def _func(k, mu, sig1, sig2, a):
        mu = float(mu)
        if k < mu:
            rval = a * np.exp(sig1 * (k - mu))
        else:
            rval = a * np.exp(-sig2 * (k - mu))
        return rval

    try:
        return [_func(x, mu, sig1, sig2, a) for x in iter(k)]
    except TypeError:
        return _func(k, mu, sig1, sig2, a)


# define container dictionaries for different types of pdfs
parameter = dict(gauss=gauss_parameter, exp=exp_parameter)
branches = dict(gauss=gauss_branches, exp=exp_branches)


class ProbabilityDensityFunction(object):
    '''
    A probability density function toolkit.
    '''

    version = __version__

    def __init__(self, x0, incr, npts, pdf, mu, params, eta=0.01):
        self.x0 = x0
        self.incr = incr
        self.npts = npts
        self.axis = create_axis(x0, incr, npts)
        self.mu = mu
        self.eta = eta
        self._pdf = pdf
        self.params = params

    def __add__(self, other):

        assert self.eta == other.eta, 'decline factors differ please use equally defined pdfs for comparison'

        eta = self.eta

        x0, incr, npts = self.commonparameter(other)

        axis = create_axis(x0, incr, npts)
        pdf_self = np.array(self.data(axis))
        pdf_other = np.array(other.data(axis))

        pdf = np.convolve(pdf_self, pdf_other, 'full') * incr

        # shift axis values for correct plotting
        npts = pdf.size
        x0 *= 2
        axis = create_axis(x0, incr, npts)
        mu = axis[np.where(pdf == max(pdf))][0]

        func, params = fit_curve(axis, pdf)

        return ProbabilityDensityFunction(x0, incr, npts, func, mu,
                                          params, eta)

    def __sub__(self, other):

        assert self.eta == other.eta, 'decline factors differ please use equally defined pdfs for comparison'

        eta = self.eta

        x0, incr, npts = self.commonparameter(other)

        axis = create_axis(x0, incr, npts)
        pdf_self = np.array(self.data(axis))
        pdf_other = np.array(other.data(axis))

        pdf = np.correlate(pdf_self, pdf_other, 'full') * incr

        # shift axis values for correct plotting
        npts = len(pdf)
        midpoint = npts / 2
        x0 = -incr * midpoint
        axis = create_axis(x0, incr, npts)
        mu = axis[np.where(pdf == max(pdf))][0]

        func, params = fit_curve(axis, pdf)

        return ProbabilityDensityFunction(x0, incr, npts, func, mu,
                                          params, eta)

    def __nonzero__(self):
        prec = self.precision(self.incr)
        data = np.array(self.data())
        gtzero = np.all(data >= 0)
        probone = bool(np.round(self.prob_gt_val(self.axis[0]), prec) == 1.)
        return bool(gtzero and probone)

    def __str__(self):
        return str([self.data()])

    @staticmethod
    def precision(incr):
        prec = int(np.ceil(np.abs(np.log10(incr)))) - 2
        return prec if prec >= 0 else 0

    def data(self, value=None):
        if value is None:
            return self._pdf(self.axis, self.params)
        return self._pdf(value, self.params)

    @property
    def eta(self):
        return self._eta

    @eta.setter
    def eta(self, value):
        self._eta = value

    @property
    def mu(self):
        return self._mu

    @mu.setter
    def mu(self, mu):
        self._mu = mu

    @property
    def axis(self):
        return self._x

    @axis.setter
    def axis(self, x):
        self._x = np.array(x)

    @classmethod
    def from_pick(self, lbound, barycentre, rbound, incr=0.001, decfact=0.01,
                  type='exp'):
        '''
        Initialize a new ProbabilityDensityFunction object.
        Takes incr, lbound, barycentre and rbound to derive x0 and the number
        of points npts for the axis vector.
        Maximum density
        is given at the barycentre and on the boundaries the function has
        declined to decfact times the maximum value. Integration of the
        function over a particular interval gives the probability for the
        variable value to be in that interval.
        '''

        # derive adequate window of definition
        margin = 2. * np.max([barycentre - lbound, rbound - barycentre])

        # find midpoint accounting also for `~obspy.UTCDateTime` object usage
        try:
            midpoint = (rbound + lbound) / 2
        except TypeError:
            try:
                midpoint = (rbound + float(lbound)) / 2
            except TypeError:
                midpoint = float(rbound + float(lbound)) / 2

        # find x0 on a grid point and sufficient npts
        was_datetime = None
        if isinstance(barycentre, UTCDateTime):
            barycentre = float(barycentre)
            was_datetime = True
        n = int(np.ceil((barycentre - midpoint) / incr))
        m = int(np.ceil((margin / incr)))
        midpoint = barycentre - n * incr
        margin = m * incr
        x0 = midpoint - margin
        npts = 2 * m

        if was_datetime:
            barycentre = UTCDateTime(barycentre)

        # calculate parameter for pdf representing function
        params = parameter[type](lbound, barycentre, rbound, decfact)

        # select pdf type
        pdf = branches[type]

        # return the object
        return ProbabilityDensityFunction(x0, incr, npts, pdf, barycentre,
                                          params, decfact)

    def broadcast(self, pdf, si, ei, data):
        try:
            pdf[si:ei] = data
        except ValueError as e:
            warnings.warn(str(e), Warning)
            return self.broadcast(pdf, si, ei, data[:-1])
        return pdf

    def expectation(self):
        '''
        returns the expectation value of the actual pdf object

        ..formula::
                mu_{\Delta t} = \int\limits_{-\infty}^\infty x \cdot f(x)dx

        :return float: rval
        '''

        #rval = 0
        #for x in self.axis:
        #    rval += x * self.data(x)
        rval = self.mu
        # Not sure about this! That might not be the barycentre.
        # However, for std calculation (next function)
        # self.mu is also used!! (LK, 02/2017) 
        return rval 

    def standard_deviation(self):
        mu = self.mu
        rval = 0
        for x in self.axis:
            rval += (x - float(mu)) ** 2 * self.data(x)
        return rval * self.incr

    def prob_lt_val(self, value):
        if value <= self.axis[0] or value > self.axis[-1]:
            raise ValueError('value out of bounds: {0}'.format(value))
        return self.prob_limits((self.axis[0], value))

    def prob_gt_val(self, value):
        if value < self.axis[0] or value >= self.axis[-1]:
            raise ValueError('value out of bounds: {0}'.format(value))
        return self.prob_limits((value, self.axis[-1]))

    def prob_limits(self, limits, oversampling=1.):
        sampling = self.incr / oversampling
        lim = np.arange(limits[0], limits[1], sampling)
        data = self.data(lim)
        min_est, max_est = 0., 0.
        for n in range(len(data) - 1):
            min_est += min(data[n], data[n + 1])
            max_est += max(data[n], data[n + 1])
        return (min_est + max_est) / 2. * sampling

    def prob_val(self, value):
        if not (self.axis[0] <= value <= self.axis[-1]):
            Warning('{0} not on axis'.format(value))
            return None
        return self.data(value) * self.incr

    def quantile(self, prob_value, eps=0.01):
        '''

        :param prob_value:
        :param eps:
        :return:
        '''
        l = self.axis[0]
        r = self.axis[-1]
        m = (r + l) / 2
        diff = prob_value - self.prob_lt_val(m)
        while abs(diff) > eps and ((r - l) > self.incr):
            if diff > 0:
                l = m
            else:
                r = m
            m = (r + l) / 2
            diff = prob_value - self.prob_lt_val(m)
        return m

    def quantile_distance(self, prob_value):
        """
        takes a probability value and and returns the distance
        between two complementary quantiles

        .. math::

            QA_\alpha = Q(1 - \alpha) - Q(\alpha)

        :param value: probability value :math:\alpha
        :type  value: float
        :return: quantile distance
        """
        if 0 >= prob_value or prob_value >= 0.5:
            raise ValueError('Value out of range.')
        ql = self.quantile(prob_value)
        qu = self.quantile(1 - prob_value)
        return qu - ql


    def quantile_dist_frac(self, x):
        """
        takes a probability value and returns the fraction of two
        corresponding quantile distances (
        :func:`pylot.core.util.pdf.ProbabilityDensityFunction
        #quantile_distance`)

        .. math::

            Q\Theta_\alpha = \frac{QA(0.5 - \alpha)}{QA(\alpha)}

        :param value: probability value :math:\alpha
        :return: quantile distance fraction
        """
        if x <= 0 or x >= 0.25:
            raise ValueError('Value out of range.')
        return self.quantile_distance(0.5-x)/self.quantile_distance(x)


    def plot(self, label=None):
        import matplotlib.pyplot as plt

        plt.plot(self.axis, self.data())
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.autoscale(axis='x', tight=True)
        if self:
            title_str = 'Probability density function '
            if label:
                title_str += label
            title_str.strip()
        else:
            title_str = 'Function not suitable as probability density function'
        plt.title(title_str)
        plt.show()

    def limits(self):
        l1 = self.x0
        r1 = l1 + self.incr * self.npts

        return l1, r1

    def cincr(self, other):
        if not self.incr == other.incr:
            raise NotImplementedError(
                'Upsampling of the lower sampled PDF not implemented yet!')
        else:
            return self.incr

    def commonlimits(self, incr, other, max_npts=1e5):
        '''
        Takes an increment incr and two left and two right limits and returns
        the left most limit and the minimum number of points needed to cover
        the whole given interval.
        :param incr:
        :param l1:
        :param l2:
        :param r1:
        :param r2:
        :param max_npts:
        :return:
        '''

        x0, r = clims(self.limits(), other.limits())

        # calculate index for rounding
        ri = self.precision(incr)

        npts = int(round(r - x0, ri) // incr)

        if npts > max_npts:
            raise ValueError('Maximum number of points exceeded:\n'
                             'max_npts - %d\n'
                             'npts - %d\n' % (max_npts, npts))

        npts = np.max([npts, self.npts, other.npts])

        if npts < self.npts or npts < other.npts:
            raise ValueError('new npts is to small')

        return x0, npts

    def commonparameter(self, other):
        assert isinstance(other, ProbabilityDensityFunction), \
            'both operands must be of type ProbabilityDensityFunction'

        incr = self.cincr(other)

        x0, npts = self.commonlimits(incr, other)

        return x0, incr, npts

