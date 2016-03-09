#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np
from obspy import UTCDateTime
from pylot.core.util.utils import find_nearest
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

    return sig1, sig2, a1, a2


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
    a = 1 / (1 / sig1 + 1 / sig2)

    return sig1, sig2, a


def gauss_branches(x, mu, sig1, sig2, a1, a2):
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
    fun_vals = []
    for k in x:
        if k < mu:
            fun_vals.append(a1 * 1 / (np.sqrt(2 * np.pi) * sig1) * np.exp(-((k - mu) / sig1) ** 2 / 2))
        else:
            fun_vals.append(a2 * 1 / (np.sqrt(2 * np.pi) * sig2) * np.exp(-((k - mu) / sig2) ** 2 / 2))
    return np.array(fun_vals)


def exp_branches(x, mu, sig1, sig2, a):
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
    fun_vals = []
    for k in x:
        if k < mu:
            fun_vals.append(a * np.exp(sig1 * (k - mu)))
        else:
            fun_vals.append(a * np.exp(-sig2 * (k - mu)))
    return np.array(fun_vals)

# define container dictionaries for different types of pdfs
parameter = dict(gauss=gauss_parameter, exp=exp_parameter)
branches = dict(gauss=gauss_branches, exp=exp_branches)


class ProbabilityDensityFunction(object):
    '''
    A probability density function toolkit.
    '''

    version = __version__

    def __init__(self, x0, incr, npts, pdf):
        self.x0 = x0
        self.incr = incr
        self.npts = npts
        self.axis = create_axis(x0, incr, npts)
        self.data = pdf

    def __add__(self, other):
        assert isinstance(other, ProbabilityDensityFunction), \
            'both operands must be of type ProbabilityDensityFunction'

        x0, incr, npts, pdf_self, pdf_other = self.rearrange(other)
        pdf = np.convolve(pdf_self, pdf_other, 'full') * incr

        # shift axis values for correct plotting
        npts = pdf.size
        x0 *= 2
        return ProbabilityDensityFunction(x0, incr, npts, pdf)

    def __sub__(self, other):
        assert isinstance(other, ProbabilityDensityFunction), \
            'both operands must be of type ProbabilityDensityFunction'

        x0, incr, npts, pdf_self, pdf_other = self.rearrange(other)

        pdf = np.correlate(pdf_self, pdf_other, 'same') * incr

        # shift axis values for correct plotting
        midpoint = npts / 2
        x0 = -incr * midpoint

        return ProbabilityDensityFunction(x0, incr, npts, pdf)

    def __nonzero__(self):
        prec = self.precision(self.incr)
        gtzero = np.any(self.data >= 0)
        probone = bool(np.round(self.prob_gt_val(self.axis[0]), prec) == 1.)
        return (gtzero and probone)

    def __str__(self):
        return str(self.data)

    @staticmethod
    def precision(incr):
        return int(np.ceil(np.abs(np.log10(incr))))

    @property
    def data(self):
        return self._pdf

    @data.setter
    def data(self, pdf):
        self._pdf = np.array(pdf)

    @property
    def axis(self):
        return self._x

    @axis.setter
    def axis(self, x):
        self._x = np.array(x)

    @classmethod
    def fromPick(self, incr, lbound, barycentre, rbound, decfact=0.01, type='gauss'):
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
        margin = 1.5 * np.max([barycentre - lbound, rbound - barycentre])

        # find midpoint accounting also for `~obspy.UTCDateTime` object usage
        try:
            midpoint = (rbound + lbound) / 2
        except TypeError:
            midpoint = (rbound + float(lbound)) / 2

        # find x0 on a grid point and sufficient npts
        n = int(np.ceil((barycentre - midpoint) / incr))
        m = int(np.ceil((margin / incr)))
        midpoint = barycentre - n * incr
        margin = m * incr
        x0 = midpoint - margin
        npts = 2 * m

        # calculate parameter for pdf representing function
        params = parameter[type](lbound, barycentre, rbound, decfact)

        # calculate pdf values
        try:
            pdf = branches[type](create_axis(x0, incr, npts), barycentre, *params)
        except TypeError as e:
            print('Warning:\n' + e.message + '\n' + 'trying timestamp instead')
            assert isinstance(barycentre, UTCDateTime), 'object not capable of' \
                                                    ' timestamp representation'
            pdf = branches[type](create_axis(x0, incr, npts),
                                 barycentre.timestamp, *params)

        # return the object
        return ProbabilityDensityFunction(x0, incr, npts, pdf)

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

        rval = 0
        for n, x in enumerate(self.axis):
            rval += x * self.data[n]
        return rval * self.incr

    def standard_deviation(self):
        mu = self.expectation()
        rval = 0
        for n, x in enumerate(self.axis):
            rval += (x - mu) ** 2 * self.data[n]
        return rval * self.incr

    def prob_lt_val(self, value):
        if value <= self.axis[0] or value > self.axis[-1]:
            raise ValueError('value out of bounds: {0}'.format(value))
        return self.data[self.axis <= value].sum() * self.incr

    def prob_gt_val(self, value):
        if value < self.axis[0] or value >= self.axis[-1]:
            raise ValueError('value out of bounds: {0}'.format(value))
        return self.data[self.axis >= value].sum() * self.incr

    def prob_val(self, value):
        if not (self.axis[0] <= value <= self.axis[-1]):
            Warning('{0} not on axis'.format(value))
            return None
        return self.data[find_nearest(self.axis, value)] * self.incr

    def plot(self):
        import matplotlib.pyplot as plt

        plt.plot(self.axis, self.data)
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.autoscale(axis='x', tight=True)
        if self:
            plt.title('Probability density function')
        else:
            plt.title('Function not suitable as probability density function')
        plt.show()

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
        # >>> manu = ProbabilityDensityFunction.fromPick(0.01, 0.3, 0.5, 0.54)
        # >>> auto = ProbabilityDensityFunction.fromPick(0.01, 0.3, 0.34, 0.54)
        # >>> manu.commonlimits(0.01, auto)
        # (

        l1 = self.x0
        r1 = l1 + self.incr * self.npts
        l2 = other.x0
        r2 = l2 + other.incr * other.npts

        if l1 < l2:
            x0 = l1
        else:
            x0 = l2

        # calculate index for rounding
        ri = self.precision(incr)

        if r1 < r2:
            npts = int(round(r2 - x0, ri) // incr)
        else:
            npts = int(round(r1 - x0, ri) // incr)

        if npts > max_npts:
            raise ValueError('Maximum number of points exceeded:\n'
                             'max_npts - %d\n'
                             'npts - %d\n' % (max_npts, npts))

        npts = np.max([npts, self.npts, other.npts])

        if npts < self.npts or npts < other.npts:
            raise ValueError('new npts is to small')

        return x0, npts


    def rearrange(self, other):
        '''
        Method rearrange takes another Probability Density Function and returns
        a new axis with mid-point 0 and covering positive and negative range
        of axis values, either containing the maximum value of both axis or
        the sum of the maxima
        :param other:
        :return:
        '''

        assert isinstance(other, ProbabilityDensityFunction), \
            'both operands must be of type ProbabilityDensityFunction'

        if not self.incr == other.incr:
            raise NotImplementedError('Upsampling of the lower sampled PDF not implemented yet!')
        else:
            incr = self.incr

        x0, npts = self.commonlimits(incr, other)

        pdf_self = np.zeros(npts)
        pdf_other = np.zeros(npts)

        x = create_axis(x0, incr, npts)

        sstart = find_nearest(x, self.x0)
        s_end = sstart + self.data.size
        ostart = find_nearest(x, other.x0)
        o_end = ostart + other.data.size

        pdf_self = self.broadcast(pdf_self, sstart, s_end, self.data)
        pdf_other = self.broadcast(pdf_other, ostart, o_end, other.data)

        return x0, incr, npts, pdf_self, pdf_other
