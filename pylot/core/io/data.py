#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import logging
import os
import fnmatch
from dataclasses import dataclass, field
from typing import List

from PySide2.QtWidgets import QMessageBox
from obspy import read, read_events, Stream, Catalog, UTCDateTime
from obspy.core.event import Event as ObsPyEvent
from obspy.io.sac import SacIOError


import pylot.core.loc.focmec as focmec
import pylot.core.loc.hypodd as hypodd
import pylot.core.loc.velest as velest
from pylot.core.io.phases import readPILOTEvent, picks_from_picksdict, \
    picksdict_from_pilot, merge_picks, PylotParameter
from pylot.core.util.errors import FormatError, OverwriteError
from pylot.core.util.event import Event
from pylot.core.util.obspyDMT_interface import qml_from_obspyDMT
from pylot.core.util.utils import fnConstructor, full_range, check4rotated, \
    check_for_gaps_and_merge, trim_station_components, check_for_nan


class Data(object):
    """
    Data container with attributes wfdata holding ~obspy.core.stream.

    :type parent: PySide2.QtWidgets.QWidget object, optional
    :param parent: A PySide2.QtWidgets.QWidget object utilized when
    called by a GUI to display a PySide2.QtWidgets.QMessageBox instead of printing
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
        elif isinstance(evtdata, str):
            try:
                cat = read_events(evtdata)
                if len(cat) != 1:
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
                                                  'implemented.')
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
                            del (old_pick)
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

    def checkEvent(self, event, fcheck, forceOverwrite=False):
        """
        Check information in supplied event and own event and replace with own
        information if no other information are given or forced by forceOverwrite
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
        Replace picks in event with own picks
        :param event: Event that supplies information for comparison
        :type event: pylot.core.util.event.Event
        :param picktype: 'auto' or 'manual' picks
        :type picktype: str
        :return:
        :rtype: None
        """
        checkflag = 1
        picks = event.picks
        # remove existing picks
        for j, pick in reversed(list(enumerate(picks))):
            try:
                if picktype in str(pick.method_id.id):
                    picks.pop(j)
                    checkflag = 2
            except AttributeError as e:
                msg = '{}'.format(e)
                print(e)
                checkflag = 0
        if checkflag > 0:
            if checkflag == 1:
                print("Write new %s picks to catalog." % picktype)
            if checkflag == 2:
                print("Found %s pick(s), remove them and append new picks to catalog." % picktype)

            # append new picks
            for pick in self.get_evt_data().picks:
                if picktype in str(pick.method_id.id):
                    picks.append(pick)

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

    def setWFData(self, fnames, fnames_alt=None, checkRotated=False, metadata=None, tstart=0, tstop=0):
        """
        Clear current waveform data and set given waveform data
        :param fnames: waveform data names to append
        :param fnames_alt: alternative data to show (e.g. synthetic/processed)
        :type fnames: list
        """
        def check_fname_exists(filenames: list) -> list:
            if filenames:
                filenames = [fn for fn in filenames if os.path.isfile(fn)]
            return filenames

        self.wfdata = Stream()
        self.wforiginal = None
        self.wf_alt = Stream()
        if tstart == tstop:
            tstart = tstop = None
        self.tstart = tstart
        self.tstop = tstop

        # remove directories
        fnames = check_fname_exists(fnames)
        fnames_alt = check_fname_exists(fnames_alt)

        if fnames is not None:
            self.appendWFData(fnames)
            if fnames_alt is not None:
                self.appendWFData(fnames_alt, alternative=True)
        else:
            return False

        # check for gaps and merge
        self.wfdata, _ = check_for_gaps_and_merge(self.wfdata)
        # check for nans
        check_for_nan(self.wfdata)
        # check for stations with rotated components
        if checkRotated and metadata is not None:
            self.wfdata = check4rotated(self.wfdata, metadata, verbosity=0)
        # trim station components to same start value
        trim_station_components(self.wfdata, trim_start=True, trim_end=False)

        # make a copy of original data
        self.wforiginal = self.getWFData().copy()
        self.dirty = False
        return True

    def appendWFData(self, fnames, alternative=False):
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

        orig_or_alternative_data = {True: self.wf_alt,
                                    False: self.wfdata}

        warnmsg = ''
        for fname in set(fnames):
            try:
                orig_or_alternative_data[alternative] += read(fname, starttime=self.tstart, endtime=self.tstop)
            except TypeError:
                try:
                    orig_or_alternative_data[alternative] += read(fname, format='GSE2', starttime=self.tstart, endtime=self.tstop)
                except Exception as e:
                    try:
                        orig_or_alternative_data[alternative] += read(fname, format='SEGY', starttime=self.tstart,
                                                                      endtime=self.tstop)
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

    def getAltWFdata(self):
        return self.wf_alt

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

    def applyEVTData(self, data, typ='pick'):
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

@dataclass
class SeismicEventData:
    event_id: str = ""
    catalog: Catalog = field(default_factory=Catalog)

    def find_event_files(self, directory: str, extensions: List[str]) -> List[str]:
        """
        Browse the directory to find event files with specified extensions.

        Parameters:
        directory (str): The directory path to search for event files.
        extensions (List[str]): List of file extensions to search for.

        Returns:
        List[str]: List of file paths that match the given extensions.

        Example:
        >>> sed = SeismicEventData()
        >>> sed.find_event_files('test_directory', ['.xml', '.quakeml']) # doctest: +SKIP
        ['test_directory/event1.xml', 'test_directory/event2.quakeml']
        """
        matches = []
        for root, _, files in os.walk(directory):
            for ext in extensions:
                for filename in fnmatch.filter(files, f'*{ext}'):
                    matches.append(os.path.join(root, filename))
        return matches

    def read_event_from_directory(self, directory: str, extensions: List[str], format: str) -> None:
        """
        Read a seismic event from the first found file in the directory with specified format.

        Parameters:
        directory (str): The directory path to search for event files.
        extensions (List[str]): List of file extensions to search for.
        format (str): The format to read the event file.

        Example:
        >>> sed = SeismicEventData()
        >>> sed.read_event_from_directory('test_directory', ['.xml', '.quakeml'], 'QUAKEML') # doctest: +SKIP
        """
        event_files = self.find_event_files(directory, extensions)
        if event_files:
            self.read_event(event_files[0], format)
        else:
            raise FileNotFoundError(f"No event files found in directory {directory} with extensions {extensions}.")

    def read_event(self, file_path: str, format: str) -> None:
        """
        Read a seismic event from a file with specified format.

        Parameters:
        file_path (str): The path to the event file.
        format (str): The format to read the event file.

        Example:
        >>> sed = SeismicEventData()
        >>> sed.read_event('test_directory/event1.xml', 'QUAKEML') # doctest: +SKIP
        """
        if os.path.exists(file_path):
            self.catalog = read_events(file_path, format=format)
            self.event_id = self.catalog[0].resource_id.id.split('/')[-1] if self.catalog else ""
        else:
            raise FileNotFoundError(f"File {file_path} does not exist.")

    def write_event(self, file_path: str, format: str) -> None:
        """
        Write the seismic event to a file with specified format.

        Parameters:
        file_path (str): The path to the output file.
        format (str): The format to write the event file.

        Example:
        >>> sed = SeismicEventData(event_id='12345')
        >>> sed.write_event('output_directory/event1.xml', 'QUAKEML') # doctest: +SKIP
        """
        self.catalog.write(file_path, format=format)

@dataclass
class WaveformData:
    stream: Stream = field(default_factory=Stream)

    def find_waveform_files(self, directory: str, extensions: List[str]) -> List[str]:
        """
        Browse the directory to find waveform files with specified extensions.

        Parameters:
        directory (str): The directory path to search for waveform files.
        extensions (List[str]): List of file extensions to search for.

        Returns:
        List[str]: List of file paths that match the given extensions.

        Example:
        >>> wd = WaveformData()
        >>> wd.find_waveform_files('test_directory', ['.mseed']) # doctest: +SKIP
        ['test_directory/waveform1.mseed']
        """
        matches = []
        for root, _, files in os.walk(directory):
            for ext in extensions:
                for filename in fnmatch.filter(files, f'*{ext}'):
                    matches.append(os.path.join(root, filename))
        return matches

    def read_waveform_from_directory(self, directory: str, extensions: List[str], format: str) -> None:
        """
        Read waveform data from the first found file in the directory with specified format.

        Parameters:
        directory (str): The directory path to search for waveform files.
        extensions (List[str]): List of file extensions to search for.
        format (str): The format to read the waveform file.

        Example:
        >>> wd = WaveformData()
        >>> wd.read_waveform_from_directory('test_directory', ['.mseed'], 'MSEED') # doctest: +SKIP
        """
        waveform_files = self.find_waveform_files(directory, extensions)
        if waveform_files:
            self.read_waveform(waveform_files[0], format)
        else:
            raise FileNotFoundError(f"No waveform files found in directory {directory} with extensions {extensions}.")

    def read_waveform(self, file_path: str, format: str) -> None:
        """
        Read waveform data from a file with specified format.

        Parameters:
        file_path (str): The path to the waveform file.
        format (str): The format to read the waveform file.

        Example:
        >>> wd = WaveformData()
        >>> wd.read_waveform('test_directory/waveform1.mseed', 'MSEED') # doctest: +SKIP
        """
        if os.path.exists(file_path):
            self.stream = read(file_path, format=format)
        else:
            raise FileNotFoundError(f"File {file_path} does not exist.")

    def write_waveform(self, file_path: str, format: str) -> None:
        """
        Write the waveform data to a file with specified format.

        Parameters:
        file_path (str): The path to the output file.
        format (str): The format to write the waveform file.

        Example:
        >>> wd = WaveformData()
        >>> wd.write_waveform('output_directory/waveform1.mseed', 'MSEED') # doctest: +SKIP
        """
        self.stream.write(file_path, format=format)

# Example usage:
# seismic_event = SeismicEventData()
# seismic_event.read_event_from_directory("path_to_directory", extensions=[".xml", ".quakeml"], format="QUAKEML")
# seismic_event.write_event("output_event_file.xml", format="QUAKEML")

# waveform_data = WaveformData()
# waveform_data.read_waveform_from_directory("path_to_directory", extensions=[".mseed"], format="MSEED")
# waveform_data.write_waveform("output_waveform_file.mseed", format="MSEED")


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
