#!/usr/bin/python
# -*- coding: utf-8 -*-


import os
import argparse
import glob
import subprocess

from obspy.core import read
from pylot.core.read.data import Data
from pylot.core.read.inputs import AutoPickParameter
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.pick.autopick import autopickevent
from pylot.core.pick.utils import writephases
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()


def autoPyLoT(inputfile):
    """
    Determine phase onsets automatically utilizing the automatic picking
    algorithms by Kueperkoch et al. 2010/2012.

    :param inputfile: path to the input file containing all parameter
    information for automatic picking (for formatting details, see.
    `~pylot.core.read.input.AutoPickParameter`
    :type inputfile: str
    :return:

    .. rubric:: Example

    """
    splash = '''************************************\n
                *********autoPyLoT starting*********\n
                The Python picking and Location Tool\n
                Version {version} 2015\n
                \n
                Authors:\n
                S. Wehling-Benatelli (Ruhr-Universität Bochum)\n
                L. Küperkoch (BESTEC GmbH, Landau i. d. Pfalz)\n
                K. Olbert (Christian-Albrechts Universität zu Kiel)\n
                ***********************************'''.format(version=_getVersionString())
    print(splash)

    # reading parameter file

    parameter = AutoPickParameter(inputfile)

    data = Data()

    # getting information on data structure

    if parameter.hasParam('datastructure'):
        datastructure = DATASTRUCTURE[parameter.getParam('datastructure')]()
        dsfields = {'root' :parameter.getParam('rootpath'),
                    'dpath' :parameter.getParam('datapath'),
                    'dbase' :parameter.getParam('database')}

        exf = ['root', 'dpath', 'dbase']

        if parameter.hasParam('eventID'):
            dsfields['eventID'] = parameter.getParam('eventID')
            exf.append('eventID')

        datastructure.modifyFields(**dsfields)
        datastructure.setExpandFields(exf)

        # check if default location routine NLLoc is available
        if parameter.hasParam('nllocbin'):
            locflag = 1
            # get NLLoc-root path
            nllocroot = parameter.getParam('nllocroot')
            # get path to NLLoc executable
            nllocbin = parameter.getParam('nllocbin')
            nlloccall = '%s/NLLoc' % nllocbin
            # get name of phase file
            phasef = parameter.getParam('phasefile')
            phasefile = '%s/obs/%s' % (nllocroot, phasef)
            # get name of NLLoc-control file
            locf = parameter.getParam('locfile')
            locfile = '%s/run/%s' % (nllocroot, locf)
            # patter of NLLoc ttimes from location grid
            ttpat = parameter.getParam('ttpatter')
            ttpatter = '%s/time/%s' % (nllocroot, ttpat)
            # patter of NLLoc-output file
            nllocoutpatter = parameter.getParam('outpatter')
        else:
            locflag = 0
            print ("!!No location routine available, autoPyLoT just picks the events without locating them!!")


        # multiple event processing
        # read each event in database
        datapath = datastructure.expandDataPath()
        if not parameter.hasParam('eventID'):
            for event in [events for events in glob.glob(os.path.join(datapath, '*')) if os.path.isdir(events)]:
                data.setWFData(glob.glob(os.path.join(datapath, event, '*')))
                print('Working on event %s' % event)
                print(data)

                wfdat = data.getWFData()  # all available streams
                ##########################################################
                # !automated picking starts here!
                picks = autopickevent(wfdat, parameter)

                ##########################################################
                # locating
                if locflag == 1:
                    # write phases to NLLoc-phase file
                    writephases(picks, 'NLLoc', phasefile)

                    # For locating the event the NLLoc-control file has to be modified!
                    # create comment line for NLLoc-control file
                    # NLLoc-output file
                    nllocout = '%s/loc/%s_%s' % (nllocroot, event, nllocoutpatter)
                    locfiles = 'LOCFILES %s NLLOC_OBS %s %s 0' % (phasefile, ttpatter, nllocout)
                    print ("Modifying  NLLoc-control file %s ..." % locfile)
                    # modification of NLLoc-control file
                    filedata = None
                    nllfile = open(locfile, 'r')
                    filedata = nllfile.read()
                    if filedata.find(locfiles) < 0:
                        # replace old command
                        filedata = filedata.replace('LOCFILES', locfiles)
                        nllfile = open(locfile, 'w')
                        nllfile.write(filedata)
                        nllfile.close()

                    # locate the event
                    subprocess.call([nlloccall, locfile])
                ##########################################################

                print '------------------------------------------'
                print '-----Finished event %s!-----' % event
                print '------------------------------------------'

        # single event processing
        else:
            data.setWFData(glob.glob(os.path.join(datapath, parameter.getParam('eventID'), '*')))
            print 'Working on event ', parameter.getParam('eventID')
            print data

            wfdat = data.getWFData()  # all available streams
            ##########################################################
            # !automated picking starts here!
            picks = autopickevent(wfdat, parameter)

            ##########################################################
            # locating
            if locflag == 1:
                # write phases to NLLoc-phase file
                writephases(picks, 'NLLoc', phasefile)

                # For locating the event the NLLoc-control file has to be modified!
                # create comment line for NLLoc-control file NLLoc-output file
                nllocout = '%s/loc/%s_%s' % (nllocroot, parameter.getParam('eventID'), nllocoutpatter)
                locfiles = 'LOCFILES %s NLLOC_OBS %s %s 0' % (phasefile, ttpatter, nllocout)
                print ("Modifying  NLLoc-control file %s ..." % locfile)
                # modification of NLLoc-control file
                filedata = None
                nllfile = open(locfile, 'r')
                filedata = nllfile.read()
                if filedata.find(locfiles) < 0:
                    # replace old command
                    filedata = filedata.replace('LOCFILES', locfiles)
                    nllfile = open(locfile, 'w')
                    nllfile.write(filedata)
                    nllfile.close()

                # locate the event
                subprocess.call([nlloccall, locfile])
            ##########################################################
            # write phase files for various location routines
            # HYPO71
            hypo71file = '%s/%s/autoPyLoT_HYPO71.pha' % (datapath, parameter.getParam('eventID'))
            writephases(picks, 'HYPO71', hypo71file)
           
                   
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
