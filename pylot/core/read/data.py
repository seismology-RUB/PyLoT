#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from obspy.core import read, Stream, UTCDateTime
from obspy import readEvents, read_inventory
from obspy.core.event import Event, ResourceIdentifier, Pick, WaveformStreamID

from pylot.core.read.io import readPILOTEvent
from pylot.core.util.utils import fnConstructor, getGlobalTimes
from pylot.core.util.errors import FormatError


class Data(object):
    '''
    Data container with attributes wfdata holding ~obspy.core.stream.

    :type parent: PySide.QtGui.QWidget object, optional
    :param parent: A PySide.QtGui.QWidget object utilized when
    called by a GUI to display a PySide.QtGui.QMessageBox instead of printing
    to standard out.
    :type wfdata: ~obspy.core.stream.Stream object, optional
    :param wfdata: ~obspy.core.stream.Stream object containing all available
    waveform data for the actual event
    :type evtdata: ~obspy.core.event.Event object, optional
    :param evtdata ~obspy.core.event.Event object containing all derived or
    loaded event. Container object holding, e.g. phase arrivals, etc.
    '''

    def __init__(self, parent=None, evtdata=None):
        self._parent = parent
        if self.getParent():
            self.comp = parent.getComponent()
        else:
            self.comp = 'Z'
            self.wfdata = Stream()
        self.newevent = False
        if evtdata is not None and isinstance(evtdata, Event):
            self.evtdata = evtdata
        elif evtdata is not None and not isinstance(evtdata, dict):
            cat = readEvents(evtdata)
            self.evtdata = cat[0]
        elif evtdata is not None:
            cat = readPILOTEvent(**evtdata)
        else:  # create an empty Event object
            self.newevent = True
            self.evtdata = Event()
            self.getEvtData().picks = []
        self.wforiginal = None
        self.cuttimes = None
        self.dirty = False

    def __str__(self):
        return str(self.wfdata)

    def getParent(self):
        return self._parent

    def isNew(self):
        return self.newevent

    def getCutTimes(self):
        if self.cuttimes is None:
            self.updateCutTimes()
        return self.cuttimes

    def updateCutTimes(self):
        self.cuttimes = getGlobalTimes(self.getWFData())

    def exportEvent(self, fnout=None, evtformat='QUAKEML'):

        from pylot.core.util.defaults import OUTPUTFORMATS

        if evtformat.strip() not in OUTPUTFORMATS.values():
            errmsg = 'selected format {0} not available'.format(evtformat)
            raise FormatError(errmsg)

        if fnout is None:
            ID = self.getID()
            # handle forbidden filenames especially on windows systems
            fnout = fnConstructor(str(ID))
        else:
            fnout = fnConstructor(str(fnout))

        evtformat = evtformat.upper().strip()

        # establish catalog object (event object has no write method)
        cat = Catalog()
        cat.append(self.getEvtData())
        # try exporting event via ObsPy
        try:
            cat.write(fnout + evtformat.lower(), format=evtformat)
        except KeyError, e:
            raise KeyError('''{0} export format
                              not implemented: {1}'''.format(evtformat, e))

    def getComp(self):
        return self.comp

    def getID(self):
        try:
            return self.evtdata.get('resource_id').id
        except:
            return None

    def filterWFData(self, kwargs):
        self.getWFData().filter(**kwargs)
        self.dirty = True

    def setWFData(self, fnames):
        self.wfdata = Stream()
        self.wforiginal = None
        if fnames is not None:
            self.appendWFData(fnames)
        self.wforiginal = self.getWFData().copy()
        self.dirty = False

    def appendWFData(self, fnames):
        assert isinstance(fnames, list), "input parameter 'fnames' is " \
                                         "supposed to be of type 'list' " \
                                         "but is actually {0}".format(type(
                                                                        fnames))
        if self.dirty:
            self.resetWFData()

        warnmsg = ''
        for fname in fnames:
            try:
                self.wfdata += read(fname)
            except TypeError:
                try:
                    self.wfdata += read(fname, format='GSE2')
                except Exception, e:
                    warnmsg += '{0}\n{1}\n'.format(fname, e)
        if warnmsg:
            warnmsg = 'WARNING: unable to read\n' + warnmsg
            print warnmsg

    def getWFData(self):
        return self.wfdata

    def getOriginalWFData(self):
        return self.wforiginal

    def resetWFData(self):
        self.wfdata = self.getOriginalWFData().copy()
        self.dirty = False

    def restituteWFData(self, fninventory):
        st = self.getWFData()
        inv = read_inventory(fninventory)
        st.attach_response(inv)
        pre_filt = (0.005, 0.006, 30.0, 35.0) # set in autoPyLoT.in
        st.remove_response(output='VEL', pre_filt=pre_filt)

    def getEvtData(self):
        return self.evtdata

    def applyEVTData(self, data, type='pick', authority_id='rub'):

        def applyPicks(picks):
            firstonset = None
            for station, onsets in picks.items():
                print 'Reading picks on station %s' % station
                for label, phase in onsets.items():
                    onset = phase['mpp']
                    epp = phase['epp']
                    lpp = phase['lpp']
                    error = phase['spe']
                    pick = Pick()
                    pick.time = onset
                    pick.time_errors.lower_uncertainty = onset - epp
                    pick.time_errors.upper_uncertainty = lpp - onset
                    pick.time_errors.uncertainty = error
                    pick.phase_hint = label
                    pick.waveform_id = WaveformStreamID(station_code=station)
                    self.getEvtData().picks.append(pick)
                    try:
                        polarity = phase['fm']
                    except KeyError, e:
                        print 'No polarity information found for %s' % phase
                    if firstonset is None or firstonset > onset:
                        firstonset = onset

            if 'smi:local' in self.getID():
                fonset_str = firstonset.strftime('%Y_%m_%d_%H_%M_%S')
                ID = ResourceIdentifier('event/' + fonset_str)
                ID.convertIDToQuakeMLURI(authority_id=authority_id)
                self.getEvtData().resource_id = ID

        def applyArrivals(arrivals):
            pass

        def applyEvent(event):
            pass

        applydata = {'pick': applyPicks,
                     'arrival': applyArrivals,
                     'event': applyEvent}

        applydata[type](data)


class GenericDataStructure(object):
    '''
    GenericDataBase type holds all information about the current data-
    base working on.
    '''

    def __init__(self, **kwargs):

        self.allowedFields = []
        self.expandFields = ['root']
        self.dsFields = {}

        self.modifyFields(**kwargs)

    def modifyFields(self, **kwargs):

        assert isinstance(kwargs, dict), 'dictionary type object expected'

        if not self.extraAllowed():
            kwargs = self.updateNotAllowed(kwargs)

        for key, value in kwargs.iteritems():
            key = str(key).lower()
            if value is not None:
                if type(value) not in (str, int, float):
                    for n, val in enumerate(value):
                        value[n] = str(val)
                else:
                    value = str(value)
            try:
                self.setFieldValue(key, value)
            except KeyError, e:
                errmsg = ''
                errmsg += 'WARNING:\n'
                errmsg += 'unable to set values for datastructure fields\n'
                errmsg += '%s; desired value was: %s\n' % (e, value)
                print errmsg

    def isField(self, key):
        return key in self.getFields().keys()

    def getFieldValue(self, key):
        if self.isField(key):
            return self.getFields()[key]
        else:
            return

    def setFieldValue(self, key, value):
        if not self.extraAllowed() and key not in self.getAllowed():
            raise KeyError
        else:
            if not self.isField(key):
                print 'creating new field "%s"' % key
            self.getFields()[key] = value

    def getFields(self):
        return self.dsFields

    def getExpandFields(self):
        return self.expandFields

    def setExpandFields(self, keys):
        expandFields = []
        for key in keys:
            if self.isField(key):
                expandFields.append(key)
        self.expandFields = expandFields

    def getAllowed(self):
        return self.allowedFields

    def extraAllowed(self):
        return not self.allowedFields

    def updateNotAllowed(self, kwargs):
        for key in kwargs:
            if key not in self.getAllowed():
                kwargs.__delitem__(key)
        return kwargs

    def hasSuffix(self):
        try:
            self.getFieldValue('suffix')
        except KeyError:
            return False
        else:
            if self.getFieldValue('suffix'):
                return True
        return False

    def expandDataPath(self):
        expandList = []
        for item in self.getExpandFields():
            expandList.append(self.getFieldValue(item))
        if self.hasSuffix():
            expandList.append('*%s' % self.getFieldValue('suffix'))
        return os.path.join(*expandList)

    def getCatalogName(self):
        return os.path.join(self.getFieldValue('root'), 'catalog.qml')


class PilotDataStructure(GenericDataStructure):
    '''
    Object containing the data access information for the old PILOT data
    structure.
    '''

    def __init__(self, **fields):

        if not fields:
            fields = {'database':'2006.01',
                      'root':'/data/Egelados/EVENT_DATA/LOCAL'}

        GenericDataStructure.__init__(self, **fields)

        self.setExpandFields(['root','database'])


class SeiscompDataStructure(GenericDataStructure):
    '''
    Dictionary containing the data access information for an SDS data archive:

    :param str dataType: Desired data type. Default: ``'waveform'``
    :param sdate, edate: Either date string or an instance of
         :class:`obspy.core.utcdatetime.UTCDateTime. Default: ``None``
    :type sdate, edate: str or UTCDateTime or None
    '''

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
        if kwargs and isinstance(kwargs, dict):
            for key, value in kwargs.iteritems():
                key = str(key)
                if type(value) not in (str, int, float):
                    for n, val in enumerate(value):
                        value[n] = str(val)
                else:
                    value = str(value)
                try:
                    self.setField(key, value)
                except KeyError, e:
                    errmsg = ''
                    errmsg += 'WARNING:\n'
                    errmsg += 'unable to set values for SDS fields\n'
                    errmsg += '%s; desired value was: %s\n' % (e, value)
                    print errmsg

    def setFieldValue(self, key, value):
        if self.isField(key):
            self.getFields()[key] = value
        else:
            print('Warning: trying to set value of non-existent field '
                  '{field}'.format(field=key))

    def getFields(self):
        return self.__sdsFields

    def getName(self):
        return self.__name

    def expandDataPath(self):
        fullChan = '{0}.{1}'.format(self.getFields()['CHAN'], self.getType())
        dataPath = os.path.join(self.getFields()['SDSdir'],
                                self.getFields()['YEAR'],
                                self.getFields()['NET'],
                                self.getFields()['STA'],
                                fullChan,
                                '*{0}'.format(self.getFields()['DAY']))
        return dataPath
