#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import argparse
import glob
import subprocess
import string
from obspy.core import read, UTCDateTime
from pylot.core.read.data import Data
from pylot.core.read.inputs import AutoPickParameter
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.pick.autopick import autopickevent, iteratepicker
from pylot.core.loc.nll import *
from pylot.core.util.version import get_git_version as _getVersionString
from pylot.core.analysis.magnitude import M0Mw

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
            ctrf = parameter.getParam('ctrfile')
            ctrfile = '%s/run/%s' % (nllocroot, ctrf)
            # pattern of NLLoc ttimes from location grid
            ttpat = parameter.getParam('ttpatter')
            # pattern of NLLoc-output file
            nllocoutpatter = parameter.getParam('outpatter')
            maxnumit = 3 # maximum number of iterations for re-picking
        else:
            locflag = 0
            print ("                 !!!              ")
            print ("!!No location routine available, autoPyLoT is running in non-location mode!!")
            print ("                 !!!              ")


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
                    picksExport(picks, 'NLLoc', phasefile)

                    # For locating the event the NLLoc-control file has to be modified!
                    evID = event[string.rfind(event, "/") + 1 : len(events) - 1]
                    nllocout = '%s_%s' % (evID, nllocoutpatter)
                    # create comment line for NLLoc-control file
                    modifyInputFile(ctrf, nllocroot, nllocout, phasef, ttpat)

                    # locate the event
                    locate(nlloccall, ctrfile)

                    # !iterative picking if traces remained unpicked or occupied with bad picks!
                    # get theoretical onset times for picks with weights >= 4
                    # in order to reprocess them using smaller time windows around theoretical onset
                    # get stations with bad onsets
                    badpicks = []
                    for key in picks:
                         if picks[key]['P']['weight'] >= 4 or picks[key]['S']['weight'] >= 4:
                             badpicks.append([key, picks[key]['P']['mpp']])

                    if len(badpicks) == 0:
                        print("autoPyLoT: No bad onsets found, thus no iterative picking necessary!") 
                    else:
                        # get theoretical P-onset times from NLLoc-location file
                        locsearch = '%s/loc/%s.????????.??????.grid?.loc.hyp' % (nllocroot, nllocout)
                        if len(glob.glob(locsearch)) > 0:
                            # get latest file if several are available
                            nllocfile = max(glob.glob(locsearch), key=os.path.getctime) 
                            nlloccounter = 0
                            while len(badpicks) > 0 and nlloccounter <= maxnumit: 
                                nlloccounter += 1
                                if nlloccounter > maxnumit:
                                    print("autoPyLoT: Number of maximum iterations reached, stop iterative picking!")
                                    break
                                print("autoPyLoT: Starting with iteration No. %d ..." % nlloccounter)
                                picks = iteratepicker(wfdat, nllocfile, picks, badpicks, parameter)
                                # write phases to NLLoc-phase file
                                picksExport(picks, 'NLLoc', phasefile)
                                # locate the event
                                locate(nlloccall, ctrfile)
                                print("autoPyLoT: Iteration No. %d finished." % nlloccounter)
                                badpicks = []
                                for key in picks:
                                     if picks[key]['P']['weight'] >= 4 or picks[key]['S']['weight'] >= 4:
                                         badpicks.append([key, picks[key]['P']['mpp']])
                                print("autoPyLoT: After iteration No. %d: %d bad onsets found ..." % (nlloccounter, \
                                       len(badpicks)))
                                if len(badpicks) == 0:
                                    print("autoPyLoT: No more bad onsets found, stop iterative picking!")
                                    break
                            # calculating seismic moment Mo and corresponding moment
                            # magnitude Mw after Hanks and Kanamori (1979) from reliable
                            # picks/waveforms 
                            for key in picks:
                                 if picks[key]['P']['weight'] < 4 and picks[key]['P']['w0'] is not None:
                                     selwf = wfdat.select(station=key)
                                     w0 = picks[key]['P']['w0']
                                     sourcepara = M0Mw(selwf, None, None,  None, w0, 5, \
                                     parameter.getParam('rho'), parameter.getParam('vp'))
                        else:
                            print("autoPyLoT: No NLLoc-location file available! Stop iteration!")
                ##########################################################
                # write phase files for various location routines
                # HYPO71
                hypo71file = '%s/%s/autoPyLoT_HYPO71.pha' % (datapath, evID)
                writephases(picks, 'HYPO71', hypo71file)

                endsplash = '''------------------------------------------\n'
                               -----Finished event %s!-----\n' 
                               ------------------------------------------'''.format \
                                (version=_getVersionString()) % evID
                print(endsplash)
                if locflag == 0:
                    print("autoPyLoT was running in non-location mode!")

        # single event processing
        else:
            data.setWFData(glob.glob(os.path.join(datapath, parameter.getParam('eventID'), '*')))
            print("Working on event "), parameter.getParam('eventID')
            print data

            wfdat = data.getWFData()  # all available streams
            ##########################################################
            # !automated picking starts here!
            picks = autopickevent(wfdat, parameter)

            ##########################################################
            # locating
            if locflag == 1:
                # write phases to NLLoc-phase file
                picksExport(picks, 'NLLoc', phasefile)

                # For locating the event the NLLoc-control file has to be modified!
                nllocout = '%s_%s' % (parameter.getParam('eventID'), nllocoutpatter)
                # create comment line for NLLoc-control file
                modifyInputFile(ctrf, nllocroot, nllocout, phasef, ttpat)

                # locate the event
                locate(nlloccall, ctrfile)
                # !iterative picking if traces remained unpicked or occupied with bad picks!
                # get theoretical onset times for picks with weights >= 4
                # in order to reprocess them using smaller time windows around theoretical onset
                # get stations with bad onsets
                badpicks = []
                for key in picks:
                     if picks[key]['P']['weight'] >= 4 or picks[key]['S']['weight'] >= 4:
                         badpicks.append([key, picks[key]['P']['mpp']])

                if len(badpicks) == 0:
                    print("autoPyLoT: No bad onsets found, thus no iterative picking necessary!") 
                else:
                    # get theoretical P-onset times from NLLoc-location file
                    locsearch = '%s/loc/%s.????????.??????.grid?.loc.hyp' % (nllocroot, nllocout)
                    if len(glob.glob(locsearch)) > 0:
                        # get latest file if several are available
                        nllocfile = max(glob.glob(locsearch), key=os.path.getctime) 
                        nlloccounter = 0
                        while len(badpicks) > 0 and nlloccounter <= maxnumit: 
                            nlloccounter += 1
                            if nlloccounter > maxnumit:
                                print("autoPyLoT: Number of maximum iterations reached, stop iterative picking!")
                                break
                            print("autoPyLoT: Starting with iteration No. %d ..." % nlloccounter)
                            picks = iteratepicker(wfdat, nllocfile, picks, badpicks, parameter)
                            # write phases to NLLoc-phase file
                            picksExport(picks, 'NLLoc', phasefile)
                            # locate the event
                            locate(nlloccall, ctrfile)
                            print("autoPyLoT: Iteration No. %d finished." % nlloccounter)
                            badpicks = []
                            for key in picks:
                                 if picks[key]['P']['weight'] >= 4 or picks[key]['S']['weight'] >= 4:
                                     badpicks.append([key, picks[key]['P']['mpp']])
                            print("autoPyLoT: After iteration No. %d: %d bad onsets found ..." % (nlloccounter, \
                                   len(badpicks)))
                            if len(badpicks) == 0:
                                print("autoPyLoT: No more bad onsets found, stop iterative picking!")
                                break
                        
                        # calculating seismic moment Mo and corresponding moment
                        # magnitude Mw after Hanks and Kanamori (1979) from reliable
                        # picks/waveforms 
                        for key in picks:
                             if picks[key]['P']['weight'] < 4 and picks[key]['P']['w0'] is not None:
                                 selwf = wfdat.select(station=key)
                                 w0 = picks[key]['P']['w0']
                                 sourcepara = M0Mw(selwf, None, None,  None, w0, 5, \
                                 parameter.getParam('rho'), parameter.getParam('vp'))
                        
                    else:
                        print("autoPyLoT: No NLLoc-location file available! Stop iteration!")
            ##########################################################
            # write phase files for various location routines
            # HYPO71
            hypo71file = '%s/%s/autoPyLoT_HYPO71.pha' % (datapath, parameter.getParam('eventID'))
            writephases(picks, 'HYPO71', hypo71file)
           
            endsplash = '''------------------------------------------\n'
                           -----Finished event %s!-----\n' 
                           ------------------------------------------'''.format \
                            (version=_getVersionString()) % parameter.getParam('eventID')
            print(endsplash)
            if locflag == 0:
                print("autoPyLoT was running in non-location mode!")
                   
    endsp = '''####################################\n
               ************************************\n
               *********autoPyLoT terminates*******\n
               The Python picking and Location Tool\n
               ************************************'''.format(version=_getVersionString())
    print(endsp)

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
