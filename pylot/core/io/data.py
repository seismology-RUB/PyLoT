#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import os
from obspy import read_events
from obspy.core import read, Stream, UTCDateTime
from obspy.core.event import Event

from pylot.core.io.phases import readPILOTEvent, picks_from_picksdict, \
    picksdict_from_pilot, merge_picks
from pylot.core.util.errors import FormatError, OverwriteError
from pylot.core.util.utils import fnConstructor, full_range


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
        if isinstance(evtdata, Event):
            pass
        elif isinstance(evtdata, dict):
            evt = readPILOTEvent(**evtdata)
            evtdata = evt
        elif isinstance(evtdata, basestring):
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
                        evtdata = Event()
                        evtdata.picks = picks_from_picksdict(picks)
                    elif 'LOC' in evtdata:
                        raise NotImplementedError('PILOT location information '
                                                  'read support not yet '
                                                  'implemeted.')
                    else:
                        raise e
                else:
                    raise e
        else:  # create an empty Event object
            self.setNew()
            evtdata = Event()
            evtdata.picks = []
        self.evtdata = evtdata
        self.wforiginal = None
        self.cuttimes = None
        self.dirty = False

    def __str__(self):
        return str(self.wfdata)

    def __add__(self, other):
        assert isinstance(other, Data), "operands must be of same type 'Data'"
        if other.isNew() and not self.isNew():
            picks_to_add = other.get_evt_data().picks
            old_picks = self.get_evt_data().picks
            for pick in picks_to_add:
                if pick not in old_picks:
                    old_picks.append(pick)
        elif not other.isNew() and self.isNew():
            new = other + self
            self.evtdata = new.get_evt_data()
        elif self.isNew() and other.isNew():
            pass
        elif self.get_evt_data().get('id') == other.get_evt_data().get('id'):
            other.setNew()
            return self + other
        else:
            raise ValueError("both Data objects have differing "
                             "unique Event identifiers")
        return self

    def getPicksStr(self):
        picks_str = ''
        for pick in self.get_evt_data().picks:
            picks_str += str(pick) + '\n'
        return picks_str

    def getParent(self):
        """


        :return:
        """
        return self._parent

    def isNew(self):
        """


        :return:
        """
        return self._new

    def setNew(self):
        self._new = True

    def getCutTimes(self):
        """


        :return:
        """
        if self.cuttimes is None:
            self.updateCutTimes()
        return self.cuttimes

    def updateCutTimes(self):
        """


        """
        self.cuttimes = full_range(self.getWFData())

    def getEventFileName(self):
        """


        :return:
        """
        ID = self.getID()
        # handle forbidden filenames especially on windows systems
        return fnConstructor(str(ID))

    def exportEvent(self, fnout, fnext='.xml'):

        """

        :param fnout:
        :param fnext:
        :raise KeyError:
        """
        from pylot.core.util.defaults import OUTPUTFORMATS

        try:
            evtformat = OUTPUTFORMATS[fnext]
        except KeyError as e:
            errmsg = '{0}; selected file extension {1} not ' \
                     'supported'.format(e, fnext)
            raise FormatError(errmsg)

        # try exporting event via ObsPy
        try:
            self.get_evt_data().write(fnout + fnext, format=evtformat)
        except KeyError as e:
            raise KeyError('''{0} export format
                              not implemented: {1}'''.format(evtformat, e))

    def getComp(self):
        """


        :return:
        """
        return self.comp

    def getID(self):
        """


        :return:
        """
        try:
            return self.evtdata.get('resource_id').id
        except:
            return None

    def filterWFData(self, kwargs):
        """

        :param kwargs:
        """
        self.getWFData().filter(**kwargs)
        self.dirty = True

    def setWFData(self, fnames):
        """

        :param fnames:
        """
        self.wfdata = Stream()
        self.wforiginal = None
        if fnames is not None:
            self.appendWFData(fnames)
        else:
            return False
        self.wforiginal = self.getWFData().copy()
        self.dirty = False
        return True

    def appendWFData(self, fnames):
        """

        :param fnames:
        """
        assert isinstance(fnames, list), "input parameter 'fnames' is " \
                                         "supposed to be of type 'list' " \
                                         "but is actually" \
                                         " {0}".format(type(fnames))
        if self.dirty:
            self.resetWFData()

        warnmsg = ''
        for fname in fnames:
            try:
                self.wfdata += read(fname)
            except TypeError:
                try:
                    self.wfdata += read(fname, format='GSE2')
                except Exception as e:
                    warnmsg += '{0}\n{1}\n'.format(fname, e)
        if warnmsg:
            warnmsg = 'WARNING: unable to read\n' + warnmsg
            print(warnmsg)

    def getWFData(self):
        """


        :return:
        """
        return self.wfdata

    def getOriginalWFData(self):
        """


        :return:
        """
        return self.wforiginal

    def resetWFData(self):
        """


        """
        self.wfdata = self.getOriginalWFData().copy()
        self.dirty = False

    def resetPicks(self):
        """


        """
        self.get_evt_data().picks = []

    def get_evt_data(self):
        """


        :return:
        """
        return self.evtdata

    def setEvtData(self, event):
        self.evtdata = event

    def applyEVTData(self, data, type='pick', authority_id='rub'):

        """

        :param data:
        :param type:
        :param authority_id:
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

            #firstonset = find_firstonset(picks)
            if self.get_evt_data().picks:
                raise OverwriteError('Actual picks would be overwritten!')
            else:
                picks = picks_from_picksdict(picks)
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
            if not self.isNew():
                self.setEvtData(event)
            else:
                # prevent overwriting original pick information
                picks =  copy.deepcopy(self.get_evt_data().picks)
                event = merge_picks(event, picks)
                # apply event information from location
                self.get_evt_data().update(event)

        applydata = {'pick': applyPicks,
                     'event': applyEvent}

        applydata[type](data)


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
