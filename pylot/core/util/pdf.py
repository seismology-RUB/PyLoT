#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()
__author__ = 'sebastianw'


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
            fun_vals.append(a1 * 1 / (np.sqrt(2 * np.pi) * sig1) * np.exp(-((k - mu) / sig1)**2 / 2 ))
        else:
            fun_vals.append(a2 * 1 / (np.sqrt(2 * np.pi) * sig2) * np.exp(-((k - mu) / sig2)**2 / 2))
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

    def __init__(self, x, pdf):
        self.axis = x
        self.data = pdf

    def __add__(self, other):
        assert isinstance(other, ProbabilityDensityFunction), \
            'both operands must be of type ProbabilityDensityFunction'

        x, pdf_self, pdf_other = self.rearrange(other)

        pdf = np.convolve(pdf_self, pdf_other, 'same') * self.delta()

        return ProbabilityDensityFunction(x, pdf)

    def __sub__(self, other):
        assert isinstance(other, ProbabilityDensityFunction), \
            'both operands must be of type ProbabilityDensityFunction'

        x, pdf_self, pdf_other = self.rearrange(other, plus=False)

        pdf = np.convolve(pdf_self, pdf_other[::-1], 'same') * self.delta()

        return ProbabilityDensityFunction(x, pdf)

    def __nonzero__(self):
        return True

    @property
    def data(self):
        return self.data

    @data.setter
    def data(self, pdf):
        self.data = np.array(pdf)

    @property
    def axis(self):
        return self.axis

    @axis.setter
    def axis(self, x):
        self.axis = np.array(x)

    def delta(self):
        return self.axis[1] - self.axis[0]

    def rearrange(self, other, plus=True):
        '''
        Method rearrange takes another Probability Density Function and returns
        a new axis with mid-point 0 and covering positive and negative range
        of axis values, either containing the maximum value of both axis or
        the sum of the maxima
        :param other:
        :param plus:
        :return:
        '''

        assert isinstance(other, ProbabilityDensityFunction), \
            'both operands must be of type ProbabilityDensityFunction'

        sd = self.delta()
        od = other.delta()

        samin = np.min(self.axis)
        oamin = np.min(other.axis)

        # test if 0 is a sampling node
        nodes_test = (not samin % sd and not oamin % od)

        # test if sampling rates match and if 0 is a sampling node
        if sd == od and nodes_test:
            if plus:
                max = np.max(self.axis) + np.max(other.axis)
            else:
                max = np.max(np.max(self.axis), np.max(other.axis))
            x = np.arange(-max, max + sd, sd)
        else:
            raise ValueError('Sampling rates do not match or nodes shifted!')

        pdf_self = np.zeros(x.size)
        pdf_other = np.zeros(x.size)

        sstart = np.where(x == samin)
        s_end = sstart + self.data.size
        ostart = np.where(x == oamin)
        o_end = ostart + other.data.size

        pdf_self[sstart:s_end] = self.data
        pdf_other[ostart:o_end] = other.data

        return x, pdf_self, pdf_other


class PickPDF(ProbabilityDensityFunction):

    def __init__(self, x, lbound, midpoint, rbound, decfact=0.01, type='gauss'):
        '''
        Initialize a new ProbabilityDensityFunction object. Takes arguments x,
        lbound, midpoint and rbound to define a probability density function
        defined on the interval of x. Maximum density is given at the midpoint
        and on the boundaries the function has declined to decfact times the
        maximum value. Integration of the function over a particular interval
        gives the probability for the variable value to lie in that interval.
        :param x: interval on which the pdf is defined
        :param lbound: left boundary
        :param midpoint: point of maximum probability density
        :param rbound: right boundary
        :param decfact: boundary decline factor
        :param type: determines the type of the probability density function's
         branches
        '''

        self.nodes = dict(te=lbound, tm=midpoint, tl=rbound, eta=decfact)
        self.type = type
        super(PickPDF, self).__init__(x, self.pdf())

    @property
    def type(self):
        return self.type

    @type.setter
    def type(self, type):
        self.type = type

    def params(self):
        return parameter[self.type](**self.nodes)

    def get(self, key):
        return self.nodes[key]

    def pdf(self):
        return branches[self.type](self.axis, self.get('tm'), *self.params())
