#!/usr/bin/python
# -*- coding: utf-8 -*-


import os
import argparse
import glob

import matplotlib.pyplot as plt
from obspy.core import read
from pylot.core.util import _getVersionString
from pylot.core.read.data import Data
from pylot.core.read.inputs import AutoPickParameter
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.pick.autopick import autopickevent

__version__ = _getVersionString()


def autoPyLoT(inputfile):
    '''
    Determine phase onsets automatically utilizing the automatic picking
    algorithms by Kueperkoch et al. 2010/2012.

    :param inputfile: path to the input file containing all parameter
    information for automatic picking (for formatting details, see.
    `~pylot.core.read.input.AutoPickParameter`
    :type inputfile: str
    :return:

    .. rubric:: Example

    '''
    print '************************************'
    print '*********autoPyLoT starting*********'
    print 'The Python picking and Location Tool'
    print '      Version ', _getVersionString(), '2015'
    print ' '
    print 'Authors:'
    print 'S. Wehling-Benatelli'
    print '   Ruhr-Universität Bochum'
    print 'L. Küperkoch'
    print '   BESTEC GmbH, Landau (Pfalz)'
    print 'K. Olbert'
    print '   Christian-Albrechts Universität Kiel'
    print '************************************'

    # reading parameter file

    parameter = AutoPickParameter(inputfile)

    data = Data()

    # getting information on data structure

    if parameter.hasParam('datastructure'):
        datastructure = DATASTRUCTURE[parameter.getParam('datastructure')]()
        dsfields = {'root':parameter.getParam('rootpath'),
                    'dpath':parameter.getParam('datapath'),
                    'dbase':parameter.getParam('database')}

        exf = ['root', 'dpath', 'dbase']

        if parameter.hasParam('eventID'):
            dsfields['eventID'] = parameter.getParam('eventID')
            exf.append('eventID')

        datastructure.modifyFields(**dsfields)
        datastructure.setExpandFields(exf)

        # get path to inventory or dataless-seed file with station meta data
        invdir = parameter.getParam('invdir')

        # multiple event processing
        # read each event in database
        datapath = datastructure.expandDataPath()
        if not parameter.hasParam('eventID'):
            for event in [events for events in glob.glob(os.path.join(datapath, '*')) if os.path.isdir(events)]:
                data.setWFData(glob.glob(os.path.join(datapath, event, '*')))
                print 'Working on event %s' %event
                print data

                wfdat = data.getWFData() # all available streams
                ##########################################################
                # !automated picking starts here!
                picks = autopickevent(wfdat, parameter)

                print '------------------------------------------'
                print '-----Finished event %s!-----' % event
                print '------------------------------------------'

        # single event processing
        else:
            data.setWFData(glob.glob(os.path.join(datapath, parameter.getParam('eventID'), '*')))
            print 'Working on event ', parameter.getParam('eventID')
            print data

            wfdat = data.getWFData() # all available streams
            ##########################################################
            # !automated picking starts here!
            picks = autopickevent(wfdat, parameter)

            print '------------------------------------------'
            print '-------Finished event %s!-------' % parameter.getParam('eventID')
            print '------------------------------------------'

    print '####################################'
    print '************************************'
    print '*********autoPyLoT terminates*******'
    print 'The Python picking and Location Tool'
    print '************************************'

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='''autoPyLoT automatically picks phase onset times using higher order statistics,
                       autoregressive prediction and AIC''')

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

    autoPyLoT(str(cla.inputfile))
