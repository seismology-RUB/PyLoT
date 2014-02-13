#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 09:58:45 2014

@author: sebastianw
"""


class AutoPickParameters(object):
    '''
    AutoPickParameters is a parameter type object capable to read and/or write
    parameter ASCII and binary files. Parameters are given by example:
    phl         S           #phaselabel
    ff1         0.1         #freqmin
    ff2         0.5         #freqmax
    tdet        6.875       #det-window_(s)_for_ar
    tpred       2.5         #pred-window_(s)_for_ar
    order       4           #order_of_ar
    fnoise      0           #noise_level_for_ar
    suppp       7           #envelopecoeff
    tolt        300         #(s)time_int_around_the_arrival_time_of_the_phase_to_analyze
    f1tpwt      4           #propfact_minfreq_secondtaper
    pickwindow  9           #length_of_pick_window
    w1          1           #length_of_smoothing_window
    w2          0.37        #cf(i-1)*(1+peps)_for_local_min
    w3          0.25        #cf(i-1)*(1+peps)_for_local_min
    tslope      0.8;2       #slope_det_window_loc_glob
    aerr        30;60       #adjusted_error_slope_fitting_loc_glob
    tsn         20;5;20;10  #length_signal_window_S/N
    proPh       Sn          #nextprominentphase
    '''

    def __init__(self, fn):
        '''
        Initialize parameter object:

        read content of an ASCII file an form a type consistent dictionary
        contain all parameters.
        '''

        parFileCont = {}
        try:
            inputFile = open(fn, 'r')
            lines = inputFile.readlines()
            for line in lines:
                parspl = line.split('\t')[:2]
                parFileCont[parspl[0]] = parspl[1]
            for key, value in parFileCont.iteritems():
                try:
                    val = int(value)
                except:
                    try:
                        val = float(value)
                    except:
                        if value.find(';') > 0:
                            vallist = value.strip().split(';')
                            val = []
                            for val0 in vallist:
                                val0 = float(val0)
                                val.append(val0)
                        else:
                            val = str(value.strip())
                parFileCont[key] = val
        except e:
            raise 'ParameterError:\n %s' % e
        finally:
            parFileCont = None
        self.__parameter = parFileCont

    def __str__(self):
        pass

    def __repr__(self):
        reprString = ''
        reprString += 'Automated picking parameter:\n\n'
        if self.__parameter is not None:
            for key, value in self.__parameter.iteritems():
                reprString += '%s:\t\t%s' % (key, value)
        else:
            reprString += 'Empty parameter dictionary.'
