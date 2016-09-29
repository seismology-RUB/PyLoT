#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import glob
import os

from obspy import read_events

import pylot.core.loc.hsat as hsat
import pylot.core.loc.nll as nll
from pylot.core.analysis.magnitude import MomentMagnitude, RichterMagnitude
from pylot.core.io.data import Data
from pylot.core.io.inputs import AutoPickParameter
from pylot.core.pick.autopick import autopickevent, iteratepicker
from pylot.core.util.dataprocessing import restitute_data, read_metadata, \
    remove_underscores
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()


def autoPyLoT(inputfile):
    """
    Determine phase onsets automatically utilizing the automatic picking
    algorithms by Kueperkoch et al. 2010/2012.

    :param inputfile: path to the input file containing all parameter
    information for automatic picking (for formatting details, see.
    `~pylot.core.io.inputs.AutoPickParameter`
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

    evt = None
    # getting information on data structure

    if parameter.hasParam('datastructure'):
        datastructure = DATASTRUCTURE[parameter.get('datastructure')]()
        dsfields = {'root': parameter.get('rootpath'),
                    'dpath': parameter.get('datapath'),
                    'dbase': parameter.get('database')}

        exf = ['root', 'dpath', 'dbase']

        if parameter.hasParam('eventID'):
            dsfields['eventID'] = parameter.get('eventID')
            exf.append('eventID')

        datastructure.modifyFields(**dsfields)
        datastructure.setExpandFields(exf)

        # check if default location routine NLLoc is available
        if parameter.hasParam('nllocbin'):
            locflag = 1
            # get NLLoc-root path
            nllocroot = parameter.get('nllocroot')
            # get path to NLLoc executable
            nllocbin = parameter.get('nllocbin')
            nlloccall = '%s/NLLoc' % nllocbin
            # get name of phase file
            phasef = parameter.get('phasefile')
            phasefile = '%s/obs/%s' % (nllocroot, phasef)
            # get name of NLLoc-control file
            ctrf = parameter.get('ctrfile')
            ctrfile = '%s/run/%s' % (nllocroot, ctrf)
            # pattern of NLLoc ttimes from location grid
            ttpat = parameter.get('ttpatter')
            # pattern of NLLoc-output file
            nllocoutpatter = parameter.get('outpatter')
            maxnumit = 3  # maximum number of iterations for re-picking
        else:
            locflag = 0
            print("                 !!!              ")
            print("!!No location routine available, autoPyLoT is running in non-location mode!!")
            print("!!No source parameter estimation possible!!")
            print("                 !!!              ")

        datapath = datastructure.expandDataPath()
        if not parameter.hasParam('eventID'):
            # multiple event processing
            # read each event in database
            events = [events for events in glob.glob(os.path.join(datapath, '*')) if os.path.isdir(events)]
        else:
            # single event processing
            events = glob.glob(os.path.join(datapath, parameter.get('eventID')))
        for event in events:
            data.setWFData(glob.glob(os.path.join(datapath, event, '*')))
            evID = os.path.split(event)[-1]
            print('Working on event %s' % event)
            print(data)
            wfdat = data.getWFData()  # all available streams
            wfdat = remove_underscores(wfdat)
            metadata =  read_metadata(parameter.get('invdir'))
            corr_dat, rest_flag = restitute_data(wfdat.copy(), *metadata)
            ##########################################################
            # !automated picking starts here!
            picks = autopickevent(wfdat, parameter)
            ##########################################################
            # locating
            if locflag == 1:
                # write phases to NLLoc-phase file
                nll.export(picks, phasefile)

                # For locating the event the NLLoc-control file has to be modified!
                nllocout = '%s_%s' % (evID, nllocoutpatter)
                # create comment line for NLLoc-control file
                nll.modify_inputs(ctrf, nllocroot, nllocout, phasef,
                                  ttpat)

                # locate the event
                nll.locate(ctrfile)

                # !iterative picking if traces remained unpicked or occupied with bad picks!
                # get theoretical onset times for picks with weights >= 4
                # in order to reprocess them using smaller time windows around theoretical onset
                # get stations with bad onsets
                badpicks = []
                for key in picks:
                    if picks[key]['P']['weight'] >= 4 or picks[key]['S']['weight'] >= 4:
                        badpicks.append([key, picks[key]['P']['mpp']])

                # TODO keep code DRY (Don't Repeat Yourself) the following part is written twice
                # suggestion: delete block and modify the later similar block to work properly

                if len(badpicks) == 0:
                    print("autoPyLoT: No bad onsets found, thus no iterative picking necessary!")
                    # get NLLoc-location file
                    locsearch = '%s/loc/%s.????????.??????.grid?.loc.hyp' % (nllocroot, nllocout)
                    if len(glob.glob(locsearch)) > 0:
                        # get latest NLLoc-location file if several are available
                        nllocfile = max(glob.glob(locsearch), key=os.path.getctime)
                        evt = read_events(nllocfile)[0]
                        # calculating seismic moment Mo and moment magnitude Mw
                        moment_mag = MomentMagnitude(corr_dat, evt, parameter.get('vp'),
                                                     parameter.get('Qp'),
                                                     parameter.get('rho'), True, 0)
                        # update pick with moment property values (w0, fc, Mo)
                        for station, props in moment_mag.moment_props.items():
                            picks[station]['P'].update(props)
                        evt = moment_mag.updated_event()
                        local_mag = RichterMagnitude(corr_dat, evt,
                                                     parameter.get('sstop'), True, 0)
                        for station, amplitude in local_mag.amplitudes.items():
                            picks[station]['S']['Ao'] = amplitude.generic_amplitude
                        evt = local_mag.updated_event()
                    else:
                        print("autoPyLoT: No NLLoc-location file available!")
                        print("No source parameter estimation possible!")
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
                            nll.export(picks, phasefile)
                            # remove actual NLLoc-location file to keep only the last
                            os.remove(nllocfile)
                            # locate the event
                            nll.locate(ctrfile)
                            print("autoPyLoT: Iteration No. %d finished." % nlloccounter)
                            # get updated NLLoc-location file
                            nllocfile = max(glob.glob(locsearch), key=os.path.getctime)
                            # check for bad picks
                            badpicks = []
                            for key in picks:
                                if picks[key]['P']['weight'] >= 4 or picks[key]['S']['weight'] >= 4:
                                    badpicks.append([key, picks[key]['P']['mpp']])
                            print("autoPyLoT: After iteration No. %d: %d bad onsets found ..." % (nlloccounter, \
                                                                                                  len(badpicks)))
                            if len(badpicks) == 0:
                                print("autoPyLoT: No more bad onsets found, stop iterative picking!")
                                nlloccounter = maxnumit
                        evt = read_events(nllocfile)[0]
                        # calculating seismic moment Mo and moment magnitude Mw
                        moment_mag = MomentMagnitude(corr_dat, evt, parameter.get('vp'),
                                                     parameter.get('Qp'),
                                                     parameter.get('rho'), True, 0)
                        # update pick with moment property values (w0, fc, Mo)
                        for station, props in moment_mag.moment_props.items():
                            picks[station]['P'].update(props)
                        evt = moment_mag.updated_event()
                        local_mag = RichterMagnitude(corr_dat, evt,
                                                     parameter.get('sstop'), True, 0)
                        for station, amplitude in local_mag.amplitudes.items():
                            picks[station]['S']['Ao'] = amplitude.generic_amplitude
                        evt = local_mag.updated_event()
                        net_mw = moment_mag.net_magnitude()
                        print("Network moment magnitude: %4.1f" % net_mw.mag)
                    else:
                        print("autoPyLoT: No NLLoc-location file available! Stop iteration!")
            ##########################################################
            # write phase files for various location routines
            # HYPO71
            hypo71file = '%s/autoPyLoT_HYPO71.pha' % event
            hsat.export(picks, hypo71file)
            data.applyEVTData(picks)
            if evt is not None:
                data.applyEVTData(evt, 'event')
            fnqml = '%s/autoPyLoT' % event
            data.exportEvent(fnqml)

            endsplash = '''------------------------------------------\n'
                           -----Finished event %s!-----\n'
                           ------------------------------------------'''.format \
                            (version=_getVersionString()) % evID
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
    from pylot.core.util.defaults import AUTOMATIC_DEFAULTS
    # parse arguments
    parser = argparse.ArgumentParser(
        description='''autoPyLoT automatically picks phase onset times using higher order statistics,
                       autoregressive prediction and AIC''')

    parser.add_argument('-i', '-I', '--inputfile', type=str,
                        action='store',
                        help='''full path to the file containing the input
                        parameters for autoPyLoT''',
                        default=AUTOMATIC_DEFAULTS
                        )
    parser.add_argument('-v', '-V', '--version', action='version',
                        version='autoPyLoT ' + __version__,
                        help='show version information and exit')

    cla = parser.parse_args()

    autoPyLoT(str(cla.inputfile))
