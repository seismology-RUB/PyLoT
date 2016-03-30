#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 12:31:25 2014

@author: sebastianw
"""

import os
from pylot.core.loc import nll
from pylot.core.loc import hsat
from pylot.core.loc import velest


def readFilterInformation(fname):
    def convert2FreqRange(*args):
        if len(args) > 1:
            return [float(arg) for arg in args]
        elif len(args) == 1:
            return float(args[0])
        return None

    filter_file = open(fname, 'r')
    filter_information = dict()
    for filter_line in filter_file.readlines():
        filter_line = filter_line.split(' ')
        for n, pos in enumerate(filter_line):
            if pos == '\n':
                filter_line[n] = ''
        filter_information[filter_line[0]] = {'filtertype': filter_line[1]
        if filter_line[1]
        else None,
                                              'order': int(filter_line[2])
                                              if filter_line[1]
                                              else None,
                                              'freq': convert2FreqRange(*filter_line[3:])
                                              if filter_line[1]
                                              else None}
    return filter_information


FILTERDEFAULTS = readFilterInformation(os.path.join(os.path.expanduser('~'),
                                                    '.pylot',
                                                    'filter.in'))

OUTPUTFORMATS = {'.xml': 'QUAKEML',
                 '.cnv': 'CNV',
                 '.obs': 'NLLOC_OBS'}

LOCTOOLS = dict(nll=nll, hsat=hsat, velest=velest)

COMPPOSITION_MAP = dict(Z=2, N=1, E=0)
COMPPOSITION_MAP['1'] = 1
COMPPOSITION_MAP['2'] = 0
COMPPOSITION_MAP['3'] = 2

COMPNAME_MAP = dict(Z='3', N='1', E='2')
