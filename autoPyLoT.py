#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import datetime
import glob
import os

import pylot.core.loc.focmec as focmec
import pylot.core.loc.hash as hash
import pylot.core.loc.hypo71 as hypo71
import pylot.core.loc.hypodd as hypodd
import pylot.core.loc.hyposat as hyposat
import pylot.core.loc.nll as nll
import pylot.core.loc.velest as velest
from obspy import read_events
from obspy.core.event import ResourceIdentifier
# from PySide.QtGui import QWidget, QInputDialog
from pylot.core.analysis.magnitude import MomentMagnitude, LocalMagnitude
from pylot.core.io.data import Data
from pylot.core.io.inputs import PylotParameter
from pylot.core.pick.autopick import autopickevent, iteratepicker
from pylot.core.util.dataprocessing import restitute_data, read_metadata
from pylot.core.util.defaults import SEPARATOR
from pylot.core.util.event import Event
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.util.utils import real_None, remove_underscores, trim_station_components, check4gaps, check4doubled, \
    check4rotated
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()


def autoPyLoT(input_dict=None, parameter=None, inputfile=None, fnames=None, eventid=None, savepath=None,
              savexml=True, station='all', iplot=0, ncores=0):
    """
    Determine phase onsets automatically utilizing the automatic picking
    algorithms by Kueperkoch et al. 2010/2012.
    :param input_dict:
    :type input_dict:
    :param parameter: PylotParameter object containing parameters used for automatic picking
    :type parameter: pylot.core.io.inputs.PylotParameter
    :param inputfile: path to the input file containing all parameter information for automatic picking
    (for formatting details, see. `~pylot.core.io.inputs.PylotParameter`
    :type inputfile: str
    :param fnames: list of data file names or None when called from GUI
    :type fnames: str
    :param eventid: event path incl. event ID (path to waveform files)
    :type eventid: str
    :param savepath: save path for autoPyLoT output, if None/"None" output will be saved in event folder
    :type savepath: str
    :param savexml: export results in XML file if True
    :type savexml: bool
    :param station: list of station names or 'all' to pick all stations
    :type station: str
    :param iplot: logical variable for plotting: 0=none, 1=partial, 2=all
    :type iplot: int
    :param ncores: number of cores used for parallel processing. Default (0) uses all available cores
    :type ncores: int
    :return: dictionary containing picks
    :rtype: dict
    """

    if ncores == 1:
        sp_info = 'autoPyLoT is running serial on 1 cores.'
    else:
        if ncores == 0:
            ncores_readable = 'all available'
        else:
            ncores_readable = ncores
        sp_info = 'autoPyLoT is running in parallel on {} cores.'.format(ncores_readable)

    splash = '''************************************\n
                *********autoPyLoT starting*********\n
                The Python picking and Location Tool\n
                Version {version} 2017\n
                \n
                Authors:\n
                L. Kueperkoch (BESTEC GmbH, Landau i. d. Pfalz)\n
                M. Paffrath (Ruhr-Universitaet Bochum)\n
                S. Wehling-Benatelli (Ruhr-Universitaet Bochum)\n
                
                {sp}
                ***********************************'''.format(version=_getVersionString(),
                                                              sp=sp_info)
    print(splash)

    parameter = real_None(parameter)
    inputfile = real_None(inputfile)
    eventid = real_None(eventid)

    fig_dict = None
    fig_dict_wadatijack = None

    locflag = 1
    if input_dict and isinstance(input_dict, dict):
        if 'parameter' in input_dict:
            parameter = input_dict['parameter']
        if 'fig_dict' in input_dict:
            fig_dict = input_dict['fig_dict']
        if 'fig_dict_wadatijack' in input_dict:
            fig_dict_wadatijack = input_dict['fig_dict_wadatijack']
        if 'station' in input_dict:
            station = input_dict['station']
        if 'fnames' in input_dict:
            fnames = input_dict['fnames']
        if 'eventid' in input_dict:
            eventid = input_dict['eventid']
        if 'iplot' in input_dict:
            iplot = input_dict['iplot']
        if 'locflag' in input_dict:
            locflag = input_dict['locflag']
        if 'savexml' in input_dict:
            savexml = input_dict['savexml']

    if not parameter:
        if inputfile:
            parameter = PylotParameter(inputfile)
            #iplot = parameter['iplot']
        else:
            infile = os.path.join(os.path.expanduser('~'), '.pylot', 'pylot.in')
            print('Using default input file {}'.format(infile))
            parameter = PylotParameter(infile)
    else:
        if not type(parameter) == PylotParameter:
            print('Wrong input type for parameter: {}'.format(type(parameter)))
            return
        if inputfile:
            print('Parameters set and input file given. Choose either of both.')
            return

    evt = None

    # reading parameter file
    if parameter.hasParam('datastructure'):
        # getting information on data structure
        datastructure = DATASTRUCTURE[parameter.get('datastructure')]()
        dsfields = {'root': parameter.get('rootpath'),
                    'dpath': parameter.get('datapath'),
                    'dbase': parameter.get('database')}

        exf = ['root', 'dpath', 'dbase']

        if parameter['eventID'] is not '*' and fnames == 'None':
            dsfields['eventID'] = parameter['eventID']
            exf.append('eventID')

        datastructure.modifyFields(**dsfields)
        datastructure.setExpandFields(exf)

        # check if default location routine NLLoc is available
        if real_None(parameter['nllocbin']) and locflag:
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
            if fnames == 'None' and parameter['eventID'] is '*':
                # multiple event processing
                # read each event in database
                events = [events for events in glob.glob(os.path.join(datapath, '*')) if os.path.isdir(events)]
            elif fnames == 'None' and parameter['eventID'] is not '*' and not type(parameter['eventID']) == list:
                # single event processing
                events = glob.glob(os.path.join(datapath, parameter['eventID']))
            elif fnames == 'None' and type(parameter['eventID']) == list:
                # multiple event processing
                events = []
                for eventID in parameter['eventID']:
                    events.append(os.path.join(datapath, eventID))
            else:
                # autoPyLoT was initialized from GUI
                events = []
                events.append(eventid)
                evID = os.path.split(eventid)[-1]
                locflag = 2
        else:
            # started in tune or interactive mode
            datapath = os.path.join(parameter['rootpath'],
                                    parameter['datapath'])
            events = []
            for eventID in eventid:
                events.append(os.path.join(datapath,
                                           parameter['database'],
                                           eventID))

        if not events:
            print('autoPyLoT: No events given. Return!')
            return

        # transform system path separator to '/'
        for index, eventpath in enumerate(events):
            eventpath = eventpath.replace(SEPARATOR, '/')
            events[index] = eventpath

        allpicks = {}
        glocflag = locflag
        for eventpath in events:
            evID = os.path.split(eventpath)[-1]
            fext = '.xml'
            filename = os.path.join(eventpath, 'PyLoT_' + evID + fext)
            try:
                data = Data(evtdata=filename)
                data.get_evt_data().path = eventpath
                print('Reading event data from filename {}...'.format(filename))
            except Exception as e:
                print('Could not read event from file {}: {}'.format(filename, e))
                data = Data()
                pylot_event = Event(eventpath)  # event should be path to event directory
                data.setEvtData(pylot_event)
            if fnames == 'None':
                data.setWFData(glob.glob(os.path.join(datapath, eventpath, '*')))
                # the following is necessary because within
                # multiple event processing no event ID is provided
                # in autopylot.in
                try:
                    parameter.get('eventID')
                except Exception:
                    now = datetime.datetime.now()
                    eventID = '%d%02d%02d%02d%02d' % (now.year,
                                                      now.month,
                                                      now.day,
                                                      now.hour,
                                                      now.minute)
                    parameter.setParam(eventID=eventID)
            else:
                data.setWFData(fnames)

                eventpath = events[0]
                # now = datetime.datetime.now()
                # evID = '%d%02d%02d%02d%02d' % (now.year,
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
            # trim components for each station to avoid problems with different trace starttimes for one station
            wfdat = check4gaps(wfdat)
            wfdat = check4doubled(wfdat)
            wfdat = trim_station_components(wfdat, trim_start=True, trim_end=False)
            metadata = read_metadata(parameter.get('invdir'))
            # rotate stations to ZNE
            wfdat = check4rotated(wfdat, metadata)
            corr_dat = None
            if locflag:
                print("Restitute data ...")
                corr_dat = restitute_data(wfdat.copy(), *metadata, ncores=ncores)
            if not corr_dat and locflag:
                locflag = 2
            print('Working on event %s. Stations: %s' % (eventpath, station))
            print(wfdat)
            ##########################################################
            # !automated picking starts here!
            fdwj = None
            if fig_dict_wadatijack:
                fdwj = fig_dict_wadatijack[evID]
            picks = autopickevent(wfdat, parameter, iplot=iplot, fig_dict=fig_dict,
                                  fig_dict_wadatijack=fdwj,
                                  ncores=ncores, metadata=metadata, origin=data.get_evt_data().origins)
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
                        # calculate seismic moment Mo and moment magnitude Mw
                        moment_mag = MomentMagnitude(corr_dat, evt, parameter.get('vp'),
                                                     parameter.get('Qp'),
                                                     parameter.get('rho'), True,
                                                     iplot)
                        # update pick with moment property values (w0, fc, Mo)
                        for stats, props in moment_mag.moment_props.items():
                            picks[stats]['P'].update(props)
                        evt = moment_mag.updated_event()
                        net_mw = moment_mag.net_magnitude()
                        print("Network moment magnitude: %4.1f" % net_mw.mag)
                        # calculate local (Richter) magntiude
                        WAscaling = parameter.get('WAscaling')
                        magscaling = parameter.get('magscaling')
                        local_mag = LocalMagnitude(corr_dat, evt,
                                                   parameter.get('sstop'),
                                                   WAscaling, True, iplot)
                        for stats, amplitude in local_mag.amplitudes.items():
                            picks[stats]['S']['Ao'] = amplitude.generic_amplitude
                        print("Local station magnitudes scaled with:")
                        print("log(Ao) + %f * log(r) + %f * r + %f" % (WAscaling[0],
                                                                       WAscaling[1],
                                                                       WAscaling[2]))
                        evt = local_mag.updated_event(magscaling)
                        net_ml = local_mag.net_magnitude(magscaling)
                        if net_ml:
                            print("Network local magnitude: %4.1f" % net_ml.mag)
                            print("Network local magnitude scaled with:")
                            print("%f * Ml + %f" % (magscaling[0], magscaling[1]))
                    else:
                        print("autoPyLoT: No NLLoc-location file available!")
                        print("No source parameter estimation possible!")
                        locflag = 9
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
                                if 'fig_dict' in input_dict:
                                    fig_dict = input_dict['fig_dict']
                                    picks = iteratepicker(wfdat, nllocfile, picks, badpicks, parameter,
                                                          fig_dict=fig_dict)
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
                            print("autoPyLoT: After iteration No. %d: %d bad onsets found ..." % (nlloccounter,
                                                                                                  len(badpicks)))
                            if len(badpicks) == 0:
                                print("autoPyLoT: No more bad onsets found, stop iterative picking!")
                                nlloccounter = maxnumit
                        evt = read_events(nllocfile)[0]
                        if locflag < 2:
                            # calculate seismic moment Mo and moment magnitude Mw
                            moment_mag = MomentMagnitude(corr_dat, evt, parameter.get('vp'),
                                                         parameter.get('Qp'),
                                                         parameter.get('rho'), True,
                                                         iplot)
                            # update pick with moment property values (w0, fc, Mo)
                            for stats, props in moment_mag.moment_props.items():
                                if stats in picks:
                                    picks[stats]['P'].update(props)
                            evt = moment_mag.updated_event()
                            net_mw = moment_mag.net_magnitude()
                            print("Network moment magnitude: %4.1f" % net_mw.mag)
                            # calculate local (Richter) magntiude
                            WAscaling = parameter.get('WAscaling')
                            magscaling = parameter.get('magscaling')
                            local_mag = LocalMagnitude(corr_dat, evt,
                                                       parameter.get('sstop'),
                                                       WAscaling, True, iplot)
                            for stats, amplitude in local_mag.amplitudes.items():
                                if stats in picks:
                                    picks[stats]['S']['Ao'] = amplitude.generic_amplitude
                            print("Local station magnitudes scaled with:")
                            print("log(Ao) + %f * log(r) + %f * r + %f" % (WAscaling[0],
                                                                           WAscaling[1],
                                                                           WAscaling[2]))
                            evt = local_mag.updated_event(magscaling)
                            net_ml = local_mag.net_magnitude(magscaling)
                            if net_ml:
                                print("Network local magnitude: %4.1f" % net_ml.mag)
                                print("Network local magnitude scaled with:")
                                print("%f * Ml + %f" % (magscaling[0], magscaling[1]))
                    else:
                        print("autoPyLoT: No NLLoc-location file available! Stop iteration!")
                        locflag = 9
            ##########################################################
            # write phase files for various location 
            # and fault mechanism calculation routines
            # ObsPy event object
            if evt is not None:
                event_id = eventpath.split('/')[-1]
                evt.resource_id = ResourceIdentifier('smi:local/' + event_id)
                data.applyEVTData(evt, 'event')
            data.applyEVTData(picks)
            if savexml:
                if savepath == 'None' or savepath is None:
                    saveEvtPath = eventpath
                else:
                    saveEvtPath = savepath
                fnqml = '%s/PyLoT_%s' % (saveEvtPath, evID)
                data.exportEvent(fnqml, fnext='.xml', fcheck=['auto', 'magnitude', 'origin'])
            if locflag == 1:
                # HYPO71
                hypo71file = '%s/PyLoT_%s_HYPO71_phases' % (eventpath, evID)
                hypo71.export(picks, hypo71file, parameter)
                # HYPOSAT
                hyposatfile = '%s/PyLoT_%s_HYPOSAT_phases' % (eventpath, evID)
                hyposat.export(picks, hyposatfile, parameter)
                # VELEST
                velestfile = '%s/PyLoT_%s_VELEST_phases.cnv' % (eventpath, evID)
                velest.export(picks, velestfile, evt, parameter)
                # hypoDD
                hypoddfile = '%s/PyLoT_%s_hypoDD_phases.pha' % (eventpath, evID)
                hypodd.export(picks, hypoddfile, parameter, evt)
                # FOCMEC
                focmecfile = '%s/PyLoT_%s_FOCMEC.in' % (eventpath, evID)
                focmec.export(picks, focmecfile, parameter, evt)
                # HASH
                hashfile = '%s/PyLoT_%s_HASH' % (eventpath, evID)
                hash.export(picks, hashfile, parameter, evt)

            endsplash = '''------------------------------------------\n'
                           -----Finished event %s!-----\n'
                           ------------------------------------------'''.format \
                           (version=_getVersionString()) % evID
            print(endsplash)
            locflag = glocflag
            if locflag == 0:
                print("autoPyLoT was running in non-location mode!")

            # save picks for current event ID to dictionary with ALL picks
            allpicks[evID] = picks

    endsp = '''####################################\n
               ************************************\n
               *********autoPyLoT terminates*******\n
               The Python picking and Location Tool\n
               ************************************'''.format(version=_getVersionString())
    print(endsp)
    return allpicks


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='''autoPyLoT automatically picks phase onset times using higher order statistics,
                       autoregressive prediction and AIC followed by locating the seismic events using 
                       NLLoc''')

    parser.add_argument('-i', '-I', '--inputfile', type=str,
                        action='store',
                        help='''full path to the file containing the input
                        parameters for autoPyLoT''')
    parser.add_argument('-p', '-P', '--iplot', type=int, 
                        action='store',
                        help='''optional, logical variable for plotting: 0=none, 1=partial, 2=all''') 
    parser.add_argument('-f', '-F', '--fnames', type=str,
                        action='store',
                        help='''optional, list of data file names''')
    parser.add_argument('-e', '--eventid', type=str,
                        action='store',
                        help='''optional, event path incl. event ID''')
    parser.add_argument('-s', '-S', '--spath', type=str,
                        action='store',
                        help='''optional, save path for autoPyLoT output''')
    parser.add_argument('-c', '-C', '--ncores', type=int,
                        action='store', default=0,
                        help='''optional, number of CPU cores used for parallel processing (default: all available(=0))''')

    cla = parser.parse_args()

    picks = autoPyLoT(inputfile=str(cla.inputfile), fnames=str(cla.fnames),
                      eventid=str(cla.eventid), savepath=str(cla.spath),
                      ncores=cla.ncores, iplot=int(cla.iplot))
