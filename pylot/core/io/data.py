#!/usr/bin/env pyth n
# -*- coding: utf-8 -*-

import copy
import os
from obspy import read_events
from obspy.core import read, Stream, UTCDateTime
from obspy.core.event import Event as ObsPyEvent
from obspy.io.sac import SacIOError

from PySide.QtGui import QMessageBox

import pylot.core.loc.velest as velest
import pylot.core.loc.focmec as focmec
import pylot.core.loc.hypodd as hypodd
from pylot.core.io.phases import readPILOTEvent, picks_from_picksdict, \
    picksdict_from_pilot, merge_picks, PylotParameter
from pylot.core.util.errors import FormatError, OverwriteError
from pylot.core.util.event import Event
from pylot.core.util.obspyDMT_interface import qml_from_obspyDMT
from pylot.core.util.utils import fnConstructor, full_range, check4rotated, \
    check4gapsAndMerge, trim_station_components

class Data(object):
    """
    Data container with attributes wfdata holding ~obspy.core.stream.

    :type parent: PySide.QtGui.QWidget object, optional
    :param parent: A PySide.QtGui.QWidget object utilized when
    called by a GUI to display a PySide.QtGui.QMessageBox instead of printing
    to standard out.
    :type evtdata: ~obspy.core.event.Event object, optional
    :param evtdata ~obspy.core.event.Event object containing all derived or
    loaded event. Container object holding, e.g. phase arrivals, etc.
    """

    def __init__(self, parent=None, evtdata=None):
        self._parent = parent
        if self.getParent():
            self.comp = parent.getComponent()
        else:
            self.comp = 'Z'
            self.wfdata = Stream()
        self._new = False
        if isinstance(evtdata, ObsPyEvent) or isinstance(evtdata, Event):
            pass
        elif isinstance(evtdata, dict):
            evt = readPILOTEvent(**evtdata)
            evtdata = evt
        elif type(evtdata) in [str, unicode]:
            try:
                cat = read_events(evtdata)
                if len(cat) is not 1:
                    raise ValueError('ambiguous event information for file: '
                                     '{file}'.format(file=evtdata))
                evtdata = cat[0]
            except TypeError as e:
                if 'Unknown format for file' in e.message:
                    if 'PHASES' in evtdata:
                        picks = picksdict_from_pilot(evtdata)
                        evtdata = ObsPyEvent()
                        evtdata.picks = picks_from_picksdict(picks)
                    elif 'LOC' in evtdata:
                        raise NotImplementedError('PILOT location information '
                                                  'read support not yet '
                                                  'implemeted.')
                    elif 'event.pkl' in evtdata:
                        evtdata = qml_from_obspyDMT(evtdata)
                    else:
                        raise e
                else:
                    raise e
        else:  # create an empty Event object
            self.setNew()
            evtdata = ObsPyEvent()
            evtdata.picks = []
        self.evtdata = evtdata
        self.wforiginal = None
        self.cuttimes = None
        self.dirty = False
        self.processed = None

    def __str__(self):
        return str(self.wfdata)

    def __add__(self, other):
        assert isinstance(other, Data), "operands must be of same type 'Data'"
        rs_id = self.get_evt_data().get('resource_id')
        rs_id_other = other.get_evt_data().get('resource_id')
        if other.isNew() and not self.isNew():
            picks_to_add = other.get_evt_data().picks
            old_picks = self.get_evt_data().picks
            wf_ids_old = [pick.waveform_id for pick in old_picks]
            for new_pick in picks_to_add:
                wf_id = new_pick.waveform_id
                if wf_id in wf_ids_old:
                    for old_pick in old_picks:
                        comparison = [old_pick.waveform_id == new_pick.waveform_id,
                                      old_pick.phase_hint == new_pick.phase_hint,
                                      old_pick.method_id == new_pick.method_id]
                        if all(comparison):
                            del(old_pick)
                old_picks.append(new_pick)
        elif not other.isNew() and self.isNew():
            new = other + self
            self.evtdata = new.get_evt_data()
        elif self.isNew() and other.isNew():
            pass
        elif rs_id == rs_id_other:
            other.setNew()
            return self + other
        else:
            raise ValueError("both Data objects have differing "
              "unique Event identifiers")
        return self

    def getPicksStr(self):
        """
        Return picks in event data
        :return: picks seperated by newlines
        :rtype: str
        """
        picks_str = ''
        for pick in self.get_evt_data().picks:
            picks_str += str(pick) + '\n'
        return picks_str

    def getParent(self):
        """
        Get PySide.QtGui.QWidget parent object
        """
        return self._parent

    def isNew(self):
        return self._new

    def setNew(self):
        self._new = True

    def getCutTimes(self):
        """
        Returns earliest start and latest end of all waveform data
        :return: minimum start time and maximum end time as a tuple
        :rtype: (UTCDateTime, UTCDateTime)
        """
        if self.cuttimes is None:
            self.updateCutTimes()
        return self.cuttimes

    def updateCutTimes(self):
        """
        Update cuttimes to contain earliest start and latest end time
        of all waveform data
        :rtype: None
        """
        self.cuttimes = full_range(self.getWFData())

    def getEventFileName(self):
        ID = self.getID()
        # handle forbidden filenames especially on windows systems
        return fnConstructor(str(ID))

    def checkEvent(self, event, fcheck, forceOverwrite=False):
        """
        Check information in supplied event and own event and replace own
        information with supplied information if own information not exiisting
        or forced by forceOverwrite
        :param event: Event that supplies information for comparison
        :type event: pylot.core.util.event.Event
        :param fcheck: check and delete existing information
        can be a str or a list of strings of ['manual', 'auto', 'origin', 'magnitude']
        :type fcheck: str, [str]
        :param forceOverwrite: Set to true to force overwrite own information. If false,
        supplied information from event is only used if there is no own information in that
        category (given in fcheck: manual, auto, origin, magnitude)
        :type forceOverwrite: bool
        :return:
        :rtype: None
        """
        if 'origin' in fcheck:
            self.replaceOrigin(event, forceOverwrite)
        if 'magnitude' in fcheck:
            self.replaceMagnitude(event, forceOverwrite)
        if 'auto' in fcheck:
            self.replacePicks(event, 'auto')
        if 'manual' in fcheck:
            self.replacePicks(event, 'manual')

    def replaceOrigin(self, event, forceOverwrite=False):
        """
        Replace own origin with the one supplied in event if own origin is not
        existing or forced by forceOverwrite = True
        :param event: Event that supplies information for comparison
        :type event: pylot.core.util.event.Event
        :param forceOverwrite: always replace own information with supplied one if true
        :type forceOverwrite: bool
        :return:
        :rtype: None
        """
        if self.get_evt_data().origins or forceOverwrite:
            if event.origins:
                print("Found origin, replace it by new origin.")
            event.origins = self.get_evt_data().origins

    def replaceMagnitude(self, event, forceOverwrite=False):
        """
        Replace own magnitude with the one supplied in event if own magnitude is not
        existing or forced by forceOverwrite = True
        :param event: Event that supplies information for comparison
        :type event: pylot.core.util.event.Event
        :param forceOverwrite: always replace own information with supplied one if true
        :type forceOverwrite: bool
        :return:
        :rtype: None
        """
        if self.get_evt_data().magnitudes or forceOverwrite:
            if event.magnitudes:
                print("Found magnitude, replace it by new magnitude")
            event.magnitudes = self.get_evt_data().magnitudes

    def replacePicks(self, event, picktype):
        """
        Replace own picks with the one in event
        :param event: Event that supplies information for comparison
        :type event: pylot.core.util.event.Event
        :param picktype: 'auto' or 'manual' picks
        :type picktype: str
        :return:
        :rtype: None
        """
        checkflag = 0
        picks = event.picks
        # remove existing picks
        for j, pick in reversed(list(enumerate(picks))):
            try:
                if picktype in str(pick.method_id.id):
                    picks.pop(j)
                    checkflag = 1
            except AttributeError as e:
                msg = '{}'.format(e)
                print(e)
                checkflag = 0
        if checkflag:
            print("Found %s pick(s), remove them and append new picks to catalog." % picktype)

            # append new picks
            for pick in self.get_evt_data().picks:
                if picktype in str(pick.method_id.id):
                    picks.append(pick)
                
    def exportEvent(self, fnout, fnext='.xml', fcheck='auto', upperErrors=None):
        """
        Export event to file
        :param fnout: basename of file
        :param fnext: file extensions xml, cnv, obs, focmec, or/and pha
        :param fcheck: check and delete existing information
        can be a str or a list of strings of ['manual', 'auto', 'origin', 'magnitude']
        """
        from pylot.core.util.defaults import OUTPUTFORMATS
        
        if not type(fcheck) == list:
            fcheck = [fcheck]

        try:
            evtformat = OUTPUTFORMATS[fnext]
        except KeyError as e:
            errmsg = '{0}; selected file extension {1} not ' \
                     'supported'.format(e, fnext)
            raise FormatError(errmsg)

        if hasattr(self.get_evt_data(), 'notes'):
            with open(os.path.join(os.path.dirname(fnout), 'notes.txt'), 'w') as notes_file:
                notes_file.write(self.get_evt_data().notes)

        # check for already existing xml-file
        if fnext == '.xml':
            if os.path.isfile(fnout + fnext):
                print("xml-file already exists! Check content ...")
                cat = read_events(fnout + fnext)
                if len(cat) > 1:
                    raise IOError('Ambigious event information in file {}'.format(fnout + fnext))
                if len(cat) < 1:
                    raise IOError('No event information in file {}'.format(fnout + fnext))
                event = cat[0]
                if not event.resource_id == self.get_evt_data().resource_id:
                    QMessageBox.warning(self, 'Warning', 'Different resource IDs!')
                    return
                self.checkEvent(event, fcheck)
                self.setEvtData(event)
                
            self.get_evt_data().write(fnout + fnext, format=evtformat)

        # try exporting event
        else:
            evtdata_org = self.get_evt_data()
            picks = evtdata_org.picks
            eventpath = evtdata_org.path
            picks_copy = copy.deepcopy(picks)
            evtdata_copy = Event(eventpath)
            evtdata_copy.picks = picks_copy

            # check for stations picked automatically as well as manually
            # Prefer manual picks!
            for i in range(len(picks)):
                if picks[i].method_id == 'manual':
                    mstation = picks[i].waveform_id.station_code
                    mstation_ext = mstation + '_'
                    for k in range(len(picks_copy)):
                        if ((picks_copy[k].waveform_id.station_code == mstation) or
                            (picks_copy[k].waveform_id.station_code == mstation_ext)) and \
                                (picks_copy[k].method_id == 'auto'):
                            del picks_copy[k]
                            break
            lendiff = len(picks) - len(picks_copy)
            if lendiff is not 0:
                print("Manual as well as automatic picks available. Prefered the {} manual ones!".format(lendiff))

            if upperErrors:
                # check for pick uncertainties exceeding adjusted upper errors
                # Picks with larger uncertainties will not be saved in output file!
                for j in range(len(picks)):
                    for i in range(len(picks_copy)):
                        if picks_copy[i].phase_hint[0] == 'P':
                            if (picks_copy[i].time_errors['upper_uncertainty'] >= upperErrors[0]) or \
                                    (picks_copy[i].time_errors['uncertainty'] is None):
                                print("Uncertainty exceeds or equal adjusted upper time error!")
                                print("Adjusted uncertainty: {}".format(upperErrors[0]))
                                print("Pick uncertainty: {}".format(picks_copy[i].time_errors['uncertainty']))
                                print("{1} P-Pick of station {0} will not be saved in outputfile".format(
                                    picks_copy[i].waveform_id.station_code,
                                    picks_copy[i].method_id))
                                print("#")
                                del picks_copy[i]
                                break
                        if picks_copy[i].phase_hint[0] == 'S':
                            if (picks_copy[i].time_errors['upper_uncertainty'] >= upperErrors[1]) or \
                                    (picks_copy[i].time_errors['uncertainty'] is None):
                                print("Uncertainty exceeds or equal adjusted upper time error!")
                                print("Adjusted uncertainty: {}".format(upperErrors[1]))
                                print("Pick uncertainty: {}".format(picks_copy[i].time_errors['uncertainty']))
                                print("{1} S-Pick of station {0} will not be saved in outputfile".format(
                                    picks_copy[i].waveform_id.station_code,
                                    picks_copy[i].method_id))
                                print("#")
                                del picks_copy[i]
                                break

            if fnext == '.obs':
                try:
                    evtdata_copy.write(fnout + fnext, format=evtformat)
                    # write header afterwards
                    evid = str(evtdata_org.resource_id).split('/')[1]
                    header = '# EQEVENT:  Label: EQ%s  Loc:  X 0.00  Y 0.00  Z 10.00  OT 0.00 \n' % evid
                    nllocfile = open(fnout + fnext)
                    l = nllocfile.readlines()
                    # Adding A0/Generic Amplitude to .obs file
                    #l2 = []
                    #for li in l:
                    #    for amp in evtdata_org.amplitudes:
                    #        if amp.waveform_id.station_code == li[0:5].strip():
                    #            li = li[0:64] + '{:0.2e}'.format(amp.generic_amplitude) + li[73:-1] + '\n'
                    #            l2.append(li)
                    #l = l2
                    nllocfile.close()
                    l.insert(0, header)
                    nllocfile = open(fnout + fnext, 'w')
                    nllocfile.write("".join(l))
                    nllocfile.close()
                except KeyError as e:
                    raise KeyError('''{0} export format
                                     not implemented: {1}'''.format(evtformat, e))
            if fnext == '.cnv':
                try:
                    velest.export(picks_copy, fnout + fnext, eventinfo=self.get_evt_data())
                except KeyError as e:
                    raise KeyError('''{0} export format
                                     not implemented: {1}'''.format(evtformat, e))
            if fnext == '_focmec.in':
                try:
                    infile = os.path.join(os.path.expanduser('~'), '.pylot', 'pylot.in')
                    print('Using default input file {}'.format(infile))
                    parameter = PylotParameter(infile)
                    focmec.export(picks_copy, fnout + fnext, parameter, eventinfo=self.get_evt_data())
                except KeyError as e:
                    raise KeyError('''{0} export format
                                     not implemented: {1}'''.format(evtformat, e))
            if fnext == '.pha':
                try:
                    infile = os.path.join(os.path.expanduser('~'), '.pylot', 'pylot.in')
                    print('Using default input file {}'.format(infile))
                    parameter = PylotParameter(infile)
                    hypodd.export(picks_copy, fnout + fnext, parameter, eventinfo=self.get_evt_data())
                except KeyError as e:
                    raise KeyError('''{0} export format
                                     not implemented: {1}'''.format(evtformat, e))

    def getComp(self):
        """
        Get component (ZNE)
        """
        return self.comp

    def getID(self):
        """
        Get unique resource id
        """
        try:
            return self.evtdata.get('resource_id').id
        except:
            return None

    def filterWFData(self, kwargs):
        """
        Filter waveform data
        :param kwargs: arguments to pass through to filter function
        """
        data = self.getWFData()
        data.detrend('linear')
        data.taper(0.02, type='cosine')
        data.filter(**kwargs)
        self.dirty = True

    def setWFData(self, fnames, fnames_syn=None, checkRotated=False, metadata=None, tstart=0, tstop=0):
        """
        Clear current waveform data and set given waveform data
        :param fnames: waveform data names to append
        :type fnames: list
        """
        self.wfdata = Stream()
        self.wforiginal = None
        self.wfsyn = Stream()
        if tstart == tstop:
            tstart = tstop = None
        self.tstart = tstart
        self.tstop = tstop

        # if obspy_dmt:
        #     wfdir = 'raw'
        #     self.processed = False
        #     for fname in fnames:
        #         if fname.endswith('processed'):
        #             wfdir = 'processed'
        #             self.processed = True
        #             break
        #     for fpath in fnames:
        #         if fpath.endswith(wfdir):
        #             wffnames = [os.path.join(fpath, fname) for fname in os.listdir(fpath)]
        #         if 'syngine' in fpath.split('/')[-1]:
        #             wffnames_syn = [os.path.join(fpath, fname) for fname in os.listdir(fpath)]
        # else:
        #     wffnames = fnames
        if fnames is not None:
            self.appendWFData(fnames)
            if fnames_syn is not None:
                self.appendWFData(fnames_syn, synthetic=True)
        else:
            return False

        # various pre-processing steps:
        # remove possible underscores in station names
        # self.wfdata = remove_underscores(self.wfdata)
        # check for gaps and merge
        self.wfdata = check4gapsAndMerge(self.wfdata)
        # check for stations with rotated components
        if checkRotated and metadata is not None:
            self.wfdata = check4rotated(self.wfdata, metadata, verbosity=0)
        # trim station components to same start value
        trim_station_components(self.wfdata, trim_start=True, trim_end=False)

        # make a copy of original data
        self.wforiginal = self.getWFData().copy()
        self.dirty = False
        return True

    def appendWFData(self, fnames, synthetic=False):
        """
        Read waveform data from fnames and append it to current wf data
        :param fnames: waveform data to append
        :type fnames: list
        """
        assert isinstance(fnames, list), "input parameter 'fnames' is " \
                                         "supposed to be of type 'list' " \
                                         "but is actually" \
                                         " {0}".format(type(fnames))
        if self.dirty:
            self.resetWFData()

        real_or_syn_data = {True: self.wfsyn,
                            False: self.wfdata}

        warnmsg = ''
        for fname in set(fnames):
            try:
                real_or_syn_data[synthetic] += read(fname, starttime=self.tstart, endtime=self.tstop)
            except TypeError:
                try:
                    real_or_syn_data[synthetic] += read(fname, format='GSE2', starttime=self.tstart, endtime=self.tstop)
                except Exception as e:
                    try:
                        real_or_syn_data[synthetic] += read(fname, format='SEGY', starttime=self.tstart, endtime=self.tstop)
                    except Exception as e:
                        warnmsg += '{0}\n{1}\n'.format(fname, e)
            except SacIOError as se:
                warnmsg += '{0}\n{1}\n'.format(fname, se)
        if warnmsg:
            warnmsg = 'WARNING in appendWFData: unable to read waveform data\n' + warnmsg
            print(warnmsg)

    def getWFData(self):
        return self.wfdata

    def getOriginalWFData(self):
        return self.wforiginal

    def getSynWFData(self):
        return self.wfsyn

    def resetWFData(self):
        """
        Set waveform data to original waveform data
        """
        if self.getOriginalWFData():
            self.wfdata = self.getOriginalWFData().copy()
        else:
            self.wfdata = Stream()
        self.dirty = False

    def resetPicks(self):
        """
        Clear all picks from event
        """
        self.get_evt_data().picks = []

    def get_evt_data(self):
        return self.evtdata

    def setEvtData(self, event):
        self.evtdata = event

    def applyEVTData(self, data, typ='pick', authority_id='rub'):
        """
        Either takes an `obspy.core.event.Event` object and applies all new
        information on the event to the actual data if typ is 'event or
        creates ObsPy pick objects and append it to the picks list from the
        PyLoT dictionary contain all picks if type is pick
        :param data: data to apply, either picks or complete event
        :type data:
        :param typ: which event data to apply, 'pick' or 'event'
        :type typ: str
        :param authority_id: (currently unused)
        :type: str
        :raise OverwriteError:
        """

        def applyPicks(picks):
            """
            Creates ObsPy pick objects and append it to the picks list from the
            PyLoT dictionary contain all picks.
            :param picks:
            :raise OverwriteError: raises an OverwriteError if the picks list is
             not empty. The GUI will then ask for a decision.
            """
            # firstonset = find_firstonset(picks)
            # check for automatic picks
            print("Writing phases to ObsPy-quakeml file")
            for key in picks:
                if not picks[key].get('P'):
                    continue
                if picks[key]['P']['picker'] == 'auto':
                    print("Existing auto-picks will be overwritten in pick-dictionary!")
                    picks = picks_from_picksdict(picks)
                    break
                else:
                    if self.get_evt_data().picks:
                        raise OverwriteError('Existing picks would be overwritten!')
                    else:
                        picks = picks_from_picksdict(picks)
                        break
            self.get_evt_data().picks = picks
            # if 'smi:local' in self.getID() and firstonset:
            #     fonset_str = firstonset.strftime('%Y_%m_%d_%H_%M_%S')
            #     ID = ResourceIdentifier('event/' + fonset_str)
            #     ID.convertIDToQuakeMLURI(authority_id=authority_id)
            #     self.get_evt_data().resource_id = ID

        def applyEvent(event):
            """
            takes an `obspy.core.event.Event` object and applies all new
            information on the event to the actual data
            :param event:
            """
            if event is None:
                print("applyEvent: Received None")
                return
            if self.isNew():
                self.setEvtData(event)
            else:
                # prevent overwriting original pick information
                event_old = self.get_evt_data()
                if not event_old.resource_id == event.resource_id:
                    print("WARNING: Missmatch in event resource id's: {} and {}".format(
                        event_old.resource_id,
                        event.resource_id))
                else:
                    picks = copy.deepcopy(event_old.picks)
                    event = merge_picks(event, picks)
                # apply event information from location
                event_old.update(event)

        applydata = {'pick': applyPicks,
                     'event': applyEvent}

        applydata[typ](data)
        self._new = False


class GenericDataStructure(object):
    """
    GenericDataBase type holds all information about the current data-
    base working on.
    """

    def __init__(self, **kwargs):

        self.allowedFields = []
        self.expandFields = ['root']
        self.dsFields = {}

        self.modifyFields(**kwargs)

    def modifyFields(self, **kwargs):

        """

        :param kwargs:
        """
        assert isinstance(kwargs, dict), 'dictionary type object expected'

        if not self.extraAllowed():
            kwargs = self.updateNotAllowed(kwargs)

        for key, value in kwargs.items():
            key = str(key).lower()
            if value is not None:
                if type(value) not in (str, int, float):
                    for n, val in enumerate(value):
                        value[n] = str(val)
                else:
                    value = str(value)
            try:
                self.setFieldValue(key, value)
            except KeyError as e:
                errmsg = ''
                errmsg += 'WARNING:\n'
                errmsg += 'unable to set values for datastructure fields\n'
                errmsg += '%s; desired value was: %s\n' % (e, value)
                print(errmsg)

    def isField(self, key):
        """

        :param key:
        :return:
        """
        return key in self.getFields().keys()

    def getFieldValue(self, key):
        """

        :param key:
        :return:
        """
        if self.isField(key):
            return self.getFields()[key]
        else:
            return

    def setFieldValue(self, key, value):
        """

        :param key:
        :param value:
        :raise KeyError:
        """
        if not self.extraAllowed() and key not in self.getAllowed():
            raise KeyError
        else:
            if not self.isField(key):
                print('creating new field "%s"' % key)
            self.getFields()[key] = value

    def getFields(self):
        """


        :return:
        """
        return self.dsFields

    def getExpandFields(self):
        """


        :return:
        """
        return self.expandFields

    def setExpandFields(self, keys):
        """

        :param keys:
        """
        expandFields = []
        for key in keys:
            if self.isField(key):
                expandFields.append(key)
        self.expandFields = expandFields

    def getAllowed(self):
        """


        :return:
        """
        return self.allowedFields

    def extraAllowed(self):
        """


        :return:
        """
        return not self.allowedFields

    def updateNotAllowed(self, kwargs):
        """

        :param kwargs:
        :return:
        """
        for key in kwargs:
            if key not in self.getAllowed():
                kwargs.__delitem__(key)
        return kwargs

    def hasSuffix(self):
        """


        :return:
        """
        try:
            self.getFieldValue('suffix')
        except KeyError:
            return False
        else:
            if self.getFieldValue('suffix'):
                return True
        return False

    def expandDataPath(self):
        """


        :return:
        """
        expandList = []
        for item in self.getExpandFields():
            expandList.append(self.getFieldValue(item))
        if self.hasSuffix():
            expandList.append('*%s' % self.getFieldValue('suffix'))
        return os.path.join(*expandList)

    def getCatalogName(self):
        """


        :return:
        """
        return os.path.join(self.getFieldValue('root'), 'catalog.qml')


class PilotDataStructure(GenericDataStructure):
    """
    Object containing the data access information for the old PILOT data
    structure.
    """

    def __init__(self, **fields):
        if not fields:
            fields = {'database': '2006.01',
                      'root': '/data/Egelados/EVENT_DATA/LOCAL'}

        GenericDataStructure.__init__(self, **fields)

        self.setExpandFields(['root', 'database'])


class ObspyDMTdataStructure(GenericDataStructure):
    """
    Object containing the data access information for the old PILOT data
    structure.
    """

    def __init__(self, **fields):
        if not fields:
            fields = {'database': '',
                      'root': ''}

        GenericDataStructure.__init__(self, **fields)

        self.setExpandFields(['root', 'database'])


class SeiscompDataStructure(GenericDataStructure):
    """
    Dictionary containing the data access information for an SDS data archive:

    :param str dataType: Desired data type. Default: ``'waveform'``
    :param sdate, edate: Either date string or an instance of
         :class:`obspy.core.utcdatetime.UTCDateTime. Default: ``None``
    :type sdate, edate: str or UTCDateTime or None
    """

    def __init__(self, rootpath='/data/SDS', dataformat='MSEED',
                 filesuffix=None, **kwargs):
        super(GenericDataStructure, self).__init__()

        edate = UTCDateTime()
        halfyear = UTCDateTime('1970-07-01')
        sdate = UTCDateTime(edate - halfyear)
        del halfyear

        year = ''
        if not edate.year == sdate.year:
            nyears = edate.year - sdate.year
            for yr in range(nyears):
                year += '{0:04d},'.format(sdate.year + yr)
            year = '{' + year[:-1] + '}'
        else:
            year = '{0:04d}'.format(sdate.year)

        # SDS fields' default values
        # definitions from
        # http://www.seiscomp3.org/wiki/doc/applications/slarchive/SDS

        self.dsFields = {'root': '/data/SDS', 'YEAR': year, 'NET': '??',
                         'STA': '????', 'CHAN': 'HH?', 'TYPE': 'D', 'LOC': '',
                         'DAY': '{0:03d}'.format(sdate.julday)
                         }
        self.modifiyFields(**kwargs)

    def modifiyFields(self, **kwargs):
        """

        :param kwargs:
        """
        if kwargs and isinstance(kwargs, dict):
            for key, value in kwargs.iteritems():
                key = str(key)
                if type(value) not in (str, int, float):
                    for n, val in enumerate(value):
                        value[n] = str(val)
                else:
                    value = str(value)
                try:
                    self.setFieldValue(key, value)
                except KeyError as e:
                    errmsg = ''
                    errmsg += 'WARNING:\n'
                    errmsg += 'unable to set values for SDS fields\n'
                    errmsg += '%s; desired value was: %s\n' % (e, value)
                    print(errmsg)

    def setFieldValue(self, key, value):
        """

        :param key:
        :param value:
        """
        if self.isField(key):
            self.getFields()[key] = value
        else:
            print('Warning: trying to set value of non-existent field '
                  '{field}'.format(field=key))

    def expandDataPath(self):
        """


        :return:
        """
        fullChan = '{0}.{1}'.format(self.getFields()['CHAN'], self.getType())
        dataPath = os.path.join(self.getFields()['SDSdir'],
                                self.getFields()['YEAR'],
                                self.getFields()['NET'],
                                self.getFields()['STA'],
                                fullChan,
                                '*{0}'.format(self.getFields()['DAY']))
        return dataPath
