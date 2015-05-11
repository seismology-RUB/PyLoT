#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from obspy.core import (read, Stream, UTCDateTime)
from obspy import readEvents
from obspy.core.event import (Event, Catalog)

from pylot.core.read import readPILOTEvent
from pylot.core.util import fnConstructor, FormatError, \
    getGlobalTimes


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
        self.wforiginal = None
        self.cuttimes = None
        self.dirty = False

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
            fnout = fnConstructor(ID)
        else:
            fnout = fnConstructor(fnout)

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
            return 'smi:bug/pylot/1234'

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

    def getEvtData(self):
        return self.evtdata

#TODO: write superclass DataStructure instead of three fully featured individual classes

class GenericDataStructure(object):
    '''
    GenericDataBase type holds all information about the current data-
    base working on.
    '''

    def __init__(self, rootpath='/data', dataformat=None, filesuffix=None,
                 **kwargs):
        self.__allowExtraFields = True
        self.__expandFields = ['root']
        self.__dsFields = {'root': rootpath}

        self.modifyFields(**kwargs)

    def modifyFields(self, **kwargs):

        assert isinstance(kwargs, dict), 'dictionary type object expected'

        for key, value in kwargs.iteritems():
            key = str(key).lower()
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
        if self.isField(key):
            self.getFields()[key] = value
        else:
            print('Warning: trying to set value of non-existent field '
                  '{field}'.format(field=key))

    def getFields(self):
        return self.__dsFields

    def getExpandFields(self):
        return self.__expandFields

    def expandDataPath(self):
        expandList = []
        for item in self.getExpandFields():
            expandList.append(self.getFieldValue(item))
        return os.path.join(*self.getExpandFields())


class PilotDataStructure(object):
    '''
    Object containing the data access information for the old PILOT data
    structure.
    '''

    def __init__(self, dataformat='GSE2', fsuffix='gse',
                 root='/data/Egelados/EVENT_DATA/LOCAL/', database='2006.01',
                 **kwargs):
        self.__name = 'PDS'
        self.dataFormat = dataformat
        self.dataType = 'waveform'
        self.__pdsFields = {'root': root,
                            'database': database,
                            'suffix': fsuffix
        }

        self.modifyFields(**kwargs)

    def modifyFields(self, **kwargs):
        if kwargs and isinstance(kwargs, dict):
            for key, value in kwargs.iteritems():
                key = str(key)
                if type(value) not in (str, int, float):
                    for n, val in enumerate(value):
                        value[n] = str(val)
                else:
                    value = str(value)
                try:
                    if key in self.getFields().keys():
                        self.getFields()[key] = value
                    else:
                        raise KeyError('unknown PDS wildcard: %s.' % key)
                except KeyError, e:
                    errmsg = ''
                    errmsg += 'WARNING:\n'
                    errmsg += 'unable to set values for PDS fields\n'
                    errmsg += '%s; desired value was: %s\n' % (e, value)
                    print errmsg

    def getType(self):
        return self.dataType

    def getFields(self):
        return self.__pdsFields

    def getName(self):
        return self.__name

    def expandDataPath(self):
        datapath = os.path.join(self.getFields()['ROOT'],
                                self.getFields()['DATABASE'])
        return datapath


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
        super(GenericDataStructure, self).__init__(rootpath, dataformat,
                                                   filesuffix, **kwargs)


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

        self.__sdsFields = {'SDSdir': '/data/SDS',  # base directory
                            'YEAR': year,  # 4 digits
                            'NET': '??',  # up to 8 characters
                            'STA': '????',  # up to 8 characters
                            'CHAN': 'HH?',  # up to 8 characters
                            'TYPE': self.getType(),  # 1 character
                            'LOC': '',  # up to 8 characters
                            'DAY': '{0:03d}'.format(sdate.julday)  # 3 digits
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

    def getType(self):
        return self.__typeOptions[self.dataType]

    def isField(self, key):
        return key in self.getFields().keys()

    def getFieldValue(self, key):
        if self.isField(key):
            return self.getFields()[key]
        else:
            return

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
