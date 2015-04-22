#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import argparse

from pylot.core.util import _getVersionString
from pylot.core.read import Data, AutoPickParameter


__version__ = _getVersionString()


def autoPyLoT(fnames, inputfile):
    parameter = AutoPickParameter(inputfile)
    data = Data()
    data.setWFData(fnames)

    print parameter.getParam('algoP')

    print(parameter)

    print(data)


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='''This program visualizes Probabilistic Power Spectral Densities.''')

    parser.add_argument('-d', '-D', '--data', type=str, action='store',
                        help='''path or url to the data to be picked''',
                        default='http://examples.obspy.org/BW.KW1..EHZ.D.2011.037'
                        )
    parser.add_argument('-i', '-I', '--inputfile', type=str,
                        action='store',
                        help='''full path to the file containing the input
                        parameters for autoPyLoT''',
                        default=os.path.join([os.path.expanduser('~'),
                                              'autoPyLoT.in'])
                        )
    parser.add_argument('-v', '-V', '--version', action='version',
                        version='autoPyLoT ' + __version__,
                        help='show version information and exit')

    cla = parser.parse_args()

    autoPyLoT(cla.data, cla.input)
