#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import argparse
import glob

from pylot.core.util import _getVersionString
from pylot.core.read import Data, AutoPickParameter
from pylot.core.pick.CharFuns import HOScf, AICcf
from pylot.core.util.structure import DATASTRUCTURE


__version__ = _getVersionString()

METHOD = {'HOS':HOScf, 'AIC':AICcf}

def autoPyLoT(inputfile):
    '''
    Determine phase onsets automatically utilizing the automatic picking
    algorithm by Kueperkoch et al. 2011.

    :param inputfile: path to the input file containing all parameter
    information for automatic picking (for formatting details, see.
    `~pylot.core.read.input.AutoPickParameter`
    :type inputfile: str
    :return:

    .. rubric:: Example

    '''

    # reading parameter file

    parameter = AutoPickParameter(inputfile)

    data = Data()

    # declaring parameter variables (only for convenience)

    meth = parameter.getParam('algoP')
    tsnr1 = parameter.getParam('tsnr1')
    tsnr2 = parameter.getParam('tsnr2')
    tnoise = parameter.getParam('pnoiselen')
    tsignal = parameter.getParam('tlim')
    order = parameter.getParam('hosorder')
    thosmw = parameter.getParam('tlta')

    # getting information on data structure

    if parameter.hasParam('datastructure'):
        datastructure = DATASTRUCTURE[parameter.getParam('datastructure')]()
        dsfields = {'root':parameter.getParam('rootpath'),
                    'dpath':parameter.getParam('datapath'),
                    'dbase':parameter.getParam('database')}

        datastructure.modifyFields(**dsfields)

        datastructure.setExpandFields(['root', 'dpath', 'dbase'])

        # process each event in database
        datapath = datastructure.expandDataPath()
        for event in glob.glob(datapath):
            data.setWFData(os.path.join(datapath, event, '*'))
            print data

            cfP = METHOD[meth](data.getWFData(), (tnoise, tsignal), thosmw, order)

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
                        default=os.path.join(os.path.expanduser('~'), '.pylot',
                                             'autoPyLoT.in')
                        )
    parser.add_argument('-v', '-V', '--version', action='version',
                        version='autoPyLoT ' + __version__,
                        help='show version information and exit')

    cla = parser.parse_args()

    autoPyLoT(cla.data, str(cla.inputfile))
