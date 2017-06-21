#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import glob
import os
import datetime
from obspy import read_events
import pylot.core.loc.hyposat as hyposat
import pylot.core.loc.hypo71 as hypo71
import pylot.core.loc.velest as velest
import pylot.core.loc.hypodd as hypodd
import pylot.core.loc.focmec as focmec
import pylot.core.loc.hash as hash
import pylot.core.loc.nll as nll
#from PySide.QtGui import QWidget, QInputDialog
from pylot.core.analysis.magnitude import MomentMagnitude, LocalMagnitude
from pylot.core.io.data import Data
from pylot.core.io.inputs import AutoPickParameter
from pylot.core.pick.autopick import autopickevent, iteratepicker
from pylot.core.util.dataprocessing import restitute_data, read_metadata, \
    remove_underscores
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()


def autoPyLoT(input_dict=None, parameter=None, inputfile=None, fnames=None, eventid=None, savepath=None, station='all', iplot=0):
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
                S. Wehling-Benatelli (Ruhr-Universitaet Bochum)\n
                L. Kueperkoch (BESTEC GmbH, Landau i. d. Pfalz)\n
                K. Olbert (Christian-Albrechts Universitaet zu Kiel)\n
                ***********************************'''.format(version=_getVersionString())
    print(splash)

    locflag = 1
    if input_dict and isinstance(input_dict, dict):
        if input_dict.has_key('parameter'):
            parameter = input_dict['parameter']
        if input_dict.has_key('fig_dict'):
            fig_dict = input_dict['fig_dict']
        if input_dict.has_key('station'):
            station = input_dict['station']
        if input_dict.has_key('fnames'):
            fnames = input_dict['fnames']
        if input_dict.has_key('iplot'):
            iplot = input_dict['iplot']
        if input_dict.has_key('locflag'):
            locflag = input_dict['locflag']

    if not parameter:
        if inputfile:
            parameter = AutoPickParameter(inputfile)
            iplot = parameter['iplot']
        else:
            print('No parameters set and no input file given. Choose either of both.')
            return
    else:
        if not type(parameter) == AutoPickParameter:
            print('Wrong input type for parameter: {}'.format(type(parameter)))
            return
        if inputfile:
            print('Parameters set and input file given. Choose either of both.')
            return
            
    data = Data()

    evt = None

    # reading parameter file
    if parameter.hasParam('datastructure'):
        # getting information on data structure
        datastructure = DATASTRUCTURE[parameter.get('datastructure')]()
        dsfields = {'root': parameter.get('rootpath'),
                    'dpath': parameter.get('datapath'),
                    'dbase': parameter.get('database')}

        exf = ['root', 'dpath', 'dbase']
        
        if parameter.hasParam('eventID') and fnames == 'None':
            dsfields['eventID'] = parameter.get('eventID')
            exf.append('eventID')

        datastructure.modifyFields(**dsfields)
        datastructure.setExpandFields(exf)

        # check if default location routine NLLoc is available
        if parameter['nllocbin'] and locflag:
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

        if not input_dict:
            # started in production mode
            datapath = datastructure.expandDataPath()
            if fnames == 'None' and not parameter['eventID']:
                # multiple event processing
                # read each event in database
                events = [events for events in glob.glob(os.path.join(datapath, '*')) if os.path.isdir(events)]
            elif fnames == 'None' and parameter['eventID']:
                # single event processing
                events = glob.glob(os.path.join(datapath, parameter.get('eventID')))
            else:
                # autoPyLoT was initialized from GUI
                events = []
                events.append(eventid)
                evID = os.path.split(eventid)[-1]
                locflag = 2
        else:
            # started in tune mode
            datapath = os.path.join(parameter['rootpath'],
                                    parameter['datapath'])
            events = []
            events.append(os.path.join(datapath,
                                       parameter['database'],
                                       parameter['eventID']))

        if not events:
            print('autoPyLoT: No events given. Return!')
            return
        
        for event in events:
            if fnames == 'None':
                data.setWFData(glob.glob(os.path.join(datapath, event, '*')))
                evID = os.path.split(event)[-1]
                # the following is necessary because within
                # multiple event processing no event ID is provided
                # in autopylot.in
                try:
                    parameter.get('eventID')
                except:
                    now = datetime.datetime.now()
                    eventID = '%d%02d%02d%02d%02d' % (now.year,
                                                      now.month,
                                                      now.day,
                                                      now.hour,
                                                      now.minute)
                    parameter.setParam(eventID=eventID)
            else:
                data.setWFData(fnames)
                
                event = events[0]
                #now = datetime.datetime.now()
                #evID = '%d%02d%02d%02d%02d' % (now.year,
                #                               now.month,
                #                               now.day,
                #                               now.hour,
                #                               now.minute)
                parameter.setParam(eventID=eventid)
            wfdat = data.getWFData()  # all available streams
            if not station == 'all':
                wfdat = wfdat.select(station=station)
                if not wfdat:
                    print('Could not find station {}. STOP!'.format(station))
                    return                      
            wfdat = remove_underscores(wfdat)
            metadata =  read_metadata(parameter.get('invdir'))
            print("Restitute data ...")
            corr_dat = restitute_data(wfdat.copy(), *metadata)
               
            print('Working on event %s. Stations: %s' % (event, station))
            print(wfdat)
            ##########################################################
            # !automated picking starts here!
            if input_dict:
                if input_dict.has_key('fig_dict'):
                    fig_dict = input_dict['fig_dict']
                    picks = autopickevent(wfdat, parameter, iplot=iplot, fig_dict=fig_dict)
            else:
                picks = autopickevent(wfdat, parameter, iplot=iplot)
            ##########################################################
            # locating
            if locflag > 0:
                # write phases to NLLoc-phase file
                nll.export(picks, phasefile, parameter)

                # For locating the event the NLLoc-control file has to be modified!
                nllocout = '%s_%s' % (evID, nllocoutpatter)
                # create comment line for NLLoc-control file
                nll.modify_inputs(ctrf, nllocroot, nllocout, phasef,
                                  ttpat)

                # locate the event
                nll.locate(ctrfile, inputfile)

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
                                                     parameter.get('rho'), True, \
                                                     iplot)
                        # update pick with moment property values (w0, fc, Mo)
                        for station, props in moment_mag.moment_props.items():
                            picks[station]['P'].update(props)
                        evt = moment_mag.updated_event()
                        local_mag = LocalMagnitude(corr_dat, evt,
                                                     parameter.get('sstop'), True,\
                                                     iplot)
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
                            if input_dict:
                                if input_dict.has_key('fig_dict'):
                                    fig_dict = input_dict['fig_dict']
                                    picks = iteratepicker(wfdat, nllocfile, picks, badpicks, parameter, fig_dict=fig_dict)
                            else:
                                picks = iteratepicker(wfdat, nllocfile, picks, badpicks, parameter)
                            # write phases to NLLoc-phase file
                            nll.export(picks, phasefile, parameter)
                            # remove actual NLLoc-location file to keep only the last
                            os.remove(nllocfile)
                            # locate the event
                            nll.locate(ctrfile, inputfile)
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
                        if locflag < 2:
                            # calculating seismic moment Mo and moment magnitude Mw
                            moment_mag = MomentMagnitude(corr_dat, evt, parameter.get('vp'),
                                                         parameter.get('Qp'),
                                                         parameter.get('rho'), True, \
                                                         iplot)
                            # update pick with moment property values (w0, fc, Mo)
                            for station, props in moment_mag.moment_props.items():
                                picks[station]['P'].update(props)
                            evt = moment_mag.updated_event()
                            local_mag = LocalMagnitude(corr_dat, evt,
                                                         parameter.get('sstop'), True, \
                                                         iplot)
                            for station, amplitude in local_mag.amplitudes.items():
                                picks[station]['S']['Ao'] = amplitude.generic_amplitude
                            evt = local_mag.updated_event()
                            net_mw = moment_mag.net_magnitude()
                            print("Network moment magnitude: %4.1f" % net_mw.mag)
                    else:
                        print("autoPyLoT: No NLLoc-location file available! Stop iteration!")
                        locflag = 9
            ##########################################################
            # write phase files for various location 
            # and fault mechanism calculation routines
            # ObsPy event object
            data.applyEVTData(picks)
            if evt is not None:
                data.applyEVTData(evt, 'event')
            fnqml = '%s/autoPyLoT' % event
            data.exportEvent(fnqml)
            # HYPO71
            hypo71file = '%s/autoPyLoT_HYPO71_phases' % event
            hypo71.export(picks, hypo71file, parameter)
            # HYPOSAT
            hyposatfile = '%s/autoPyLoT_HYPOSAT_phases' % event
            hyposat.export(picks, hyposatfile, parameter)
            if locflag == 1:
            	# VELEST
            	velestfile = '%s/autoPyLoT_VELEST_phases.cnv' % event
            	velest.export(picks, velestfile, parameter, evt)
            	# hypoDD
            	hypoddfile = '%s/autoPyLoT_hypoDD_phases.pha' % event
            	hypodd.export(picks, hypoddfile, parameter, evt)
            	# FOCMEC
            	focmecfile = '%s/autoPyLoT_FOCMEC.in' % event
            	focmec.export(picks, focmecfile, parameter, evt)
            	# HASH
            	hashfile = '%s/autoPyLoT_HASH' % event
            	hash.export(picks, hashfile, parameter, evt)

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
    return picks


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='''autoPyLoT automatically picks phase onset times using higher order statistics,
                       autoregressive prediction and AIC followed by locating the seismic events using 
                       NLLoc''')

    #parser.add_argument('-d', '-D', '--input_dict', type=str,
    #                    action='store',
    #                    help='''optional, dictionary containing processing parameters''')
    #parser.add_argument('-p', '-P', '--parameter', type=str,
    #                    action='store',
    #                    help='''parameter file, default=None''')
    parser.add_argument('-i', '-I', '--inputfile', type=str,
                        action='store',
                        help='''full path to the file containing the input
                        parameters for autoPyLoT''') 
    parser.add_argument('-f', '-F', '--fnames', type=str,
                        action='store',
                        help='''optional, list of data file names''')
    parser.add_argument('-e', '-E', '--eventid', type=str,
                        action='store',
                        help='''optional, event path incl. event ID''')
    parser.add_argument('-s', '-S', '--spath', type=str,
                        action='store',
                        help='''optional, save path for autoPyLoT output''')
    #parser.add_argument('-v', '-V', '--version', action='version',
    #                    version='autoPyLoT ' + __version__,
    #                    help='show version information and exit')

    cla = parser.parse_args()
    
    picks = autoPyLoT(inputfile=str(cla.inputfile), fnames=str(cla.fnames), 
                      eventid=str(cla.eventid), savepath=str(cla.spath))
