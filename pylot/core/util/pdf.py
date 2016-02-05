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

        self.axis = np.array(x)
        self.nodes = dict(lbound=lbound, midpoint=midpoint, rbound=rbound, eta=decfact)
        self.type = type

    def __add__(self, other):
        pass

    def __sub__(self, other):
        pass

    @property
    def type(self):
        return self.type

    @type.setter
    def type(self, type):
        self.type = type

    def params(self):
        return parameter[self.type](**self.nodes)

    def data(self):
        return branches[self.type](self.axis, self.nodes['midpoint'], *self.params())
