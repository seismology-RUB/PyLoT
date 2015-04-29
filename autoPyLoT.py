#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import argparse

from pylot.core.util import _getVersionString
from pylot.core.read import Data, AutoPickParameter
from pylot.core.pick.CharFuns import HOScf, AICcf
from pylot.core.util.structure import DATASTRUCTURE


__version__ = _getVersionString()

METHOD = {'HOS':HOScf, 'AIC':AICcf}

def autoPyLoT(fnames, inputfile):
    '''
    Determine phase onsets automatically utilizing the automatic picking
    algorithm by Kueperkoch et al. 2011.

    :param fnames: list of strings containing the paths or urls to the
    waveform data to be picked
    :type fnames: list
    :param inputfile: path to the input file containing all parameter
    information for automatic picking (for formatting details, see.
    `~pylot.core.read.input.AutoPickParameter`
    :type inputfile: str
    :return:

    .. rubric:: Example

    '''
    parameter = AutoPickParameter(inputfile)

    data = Data()

    cfP = METHOD[parameter.getParam('algoP')]()

    if parameter.hasParam('datastructure'):
        datastructure = DATASTRUCTURE[parameter.getParam('datastructure')]()

    def expandSDS(datastructure):
        return datastructure.expandDataPath()

    def expandPDS(datastructure):
        return os.path.join(datastructure.expandDataPath(),'*')

    dsem = {'PDS':expandPDS, 'SDS':expandSDS}

    expandMethod = dsem[datastructure.getName()]

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
