#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def crosscorrsingle(wf1, wf2, taumax):
    '''
    Calculates the crosscorrelation between two waveforms with a defined maximum timedifference.
    :param wf1: first waveformdata
    :type wf1: list
    :param wf2: second waveformdata
    :type wf2: list
    :param taumax: maximum time difference between waveforms
    :type taumax: positive integer
    :return: returns the crosscorrelation funktion 'c' and the lagvector 'l'
    :rtype: c and l are lists
    '''
    N = len(wf1)
    c = np.zeros(2 * taumax - 1)
    l = np.zeros(2 * taumax - 1)
    for tau in range(taumax):
        Cxyplus = 0
        Cxyminus = 0
        for n in range(N - tau):
            Cxy1plus = wf1[n] * wf2[n + tau]
            Cxy1minus = wf1[n + tau] * wf2[n]
            Cxyplus = Cxyplus + Cxy1plus
            Cxyminus = Cxyminus + Cxy1minus

        c[(taumax - 1) - tau] = Cxyminus
        c[(taumax - 1) + tau] = Cxyplus
        l[(taumax - 1) - tau] = -tau
        l[(taumax - 1) + tau] = tau
    return c, l


def crosscorrnormcalc(weights, wfs):
    '''
    crosscorrnormcalc - function that calculates the normalization for the
    cross correlation carried out by 'wfscrosscorr'
    :param weights: weighting factors for the single components
    :type weights: tuple
    :param wfs: tuple of `~numpy.array` object containing waveform data
    :type wfs: tuple
    :return: a floating point number yielding the by 'weights' weighted energy
    of the waveforms in 'wfs'
    :rtype: float
    '''

    # check if the parameters are of the right type
    if not isinstance(weights, tuple):
        raise TypeError("type of 'weight' should be 'tuple', but is {0}".format(
            type(weights)))
    if not isinstance(wfs, tuple):
        raise TypeError(
            "type of parameter 'wfs' should be 'tuple', but is {0}".format(
                type(wfs)))
    sqrsumwfs = 0.
    for n, wf in enumerate(wfs):
        sqrsumwf = np.sum(weights[n] ** 2. * wf ** 2.)
        sqrsumwfs += sqrsumwf
    return np.sqrt(sqrsumwfs)


def wfscrosscorr(weights, wfs, taumax):
    '''
    wfscrosscorr - function that calculates successive cross-correlations from a set of waveforms stored in a matrix

        base formula is:
            C(i)=SUM[p=1:nComponent](eP(p)*(SUM[n=1:N]APp(x,n)*APp(y,n+i)))/(SQRT(SUM[p=1:nComponent]eP(p)^2*(SUM[n=1:N](APp(x,n)^2)))*SQRT(SUM[p=1:nComponent]eP(p)^2*(SUM[n=1:N]APp(y,n)^2)))
                whereas
                    nComponent is the number of components
                    N is the number of samples
                    i is the lag-index

    input:
        APp		rowvectors containing the waveforms of each component p for which the cross-correlation is calculated
        tPp		rowvectros containing times
        eP		vector containing the weighting factors for the components (maxsize = [1x3])

    output:
        C		cross-correlation function
        L		lag-vector

    author(s):

    SWB 26.01.2010 as arranged with Thomas Meier and Monika Bischoff

    :param weights: weighting factors for the single components
    :type weights: tuple
    :param wfs: tuple of `~numpy.array` object containing waveform data
    :type wfs: tuple
    :param taumax: maximum time difference
    :type taumax: positive integer
    :return: returns cross correlation function normalized by the waveform energy
    '''

    ccnorm = 0.
    ccnorm = crosscorrnormcalc(weights, wfs[0])
    ccnorm *= crosscorrnormcalc(weights, wfs[1])

    c = 0.
    for n in range(len(wfs)):
        cc, l = crosscorrsingle(wfs[0][n], wfs[1][n], taumax)
        c += cc
    return c / ccnorm, l
