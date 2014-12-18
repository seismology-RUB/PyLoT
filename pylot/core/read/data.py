#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from PySide.QtGui import QMessageBox
from obspy.core import (read, Stream)
from obspy import readEvents
from obspy.core.event import (Event, Catalog)
from pylot.core.util import fnConstructor
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
        try:
            if parent:
                self.wfdata = read(parent.fnames)
            else:
                self.wfdata = read()
        except IOError, e:
            msg = 'An I/O error occured while loading data!'
            inform = 'Variable wfdata will be empty.'
            details = '{0}'.format(e)
            if parent is not None:
                warnio = QMessageBox(parent=parent)
                warnio.setText(msg)
                warnio.setDetailedText(details)
                warnio.setInformativeText(inform)
                warnio.setStandarButtons(QMessageBox.Ok)
                warnio.setIcon(QMessageBox.Warning)
            else:
                print msg, '\n', details
            self.wfdata = Stream()
        else:
            self.wfdata = Stream()
        self.newevent = False
        if evtdata is not None and isinstance(evtdata, Event):
            self.evtdata = evtdata
        elif evtdata is not None and not evtdata.endswith('.mat'):
            cat = readEvents(evtdata)
            self.evtdata = cat[0]
        elif evtdata is not None:
            cat = readMatPhases(evtdata)
        else:  # create an empty Event object
            self.newevent = True
            self.evtdata = Event()
        self.orig = self.wfdata.copy()

    def isNew(self):
        return self.newevent

    def readMatPhases(self, fname):
        pass

    def exportEvent(self, fnout=None, evtformat='QUAKEML'):

        from pylot.core.util.defaults import OUTPUTFORMATS

        if evtformat.strip() not in OUTPUTFORMATS.values():
            errmsg = 'selected format {0} not available'.format(evtformat)
            raise FormatError(errmsg)

        if fnout is None:
            ID = self.evtdata.getEventID()
        # handle forbidden filenames especially on windows systems
        fnout = fnConstructor(ID)

        evtformat = evtformat.upper().strip()

        # establish catalog object (event object has no write method)
        cat = Catalog()
        cat.append(self.event)
        # try exporting event via ObsPy
        try:
            cat.write(fnout + evtformat.lower(), format=evtformat)
        except KeyError, e:
            raise KeyError('''{0} export format 
                              not implemented: {1}'''.format(evtformat, e))

    def plotData(self, widget):
        time_ax = np.arange(0, len(self.wfdata[0].data)/self.wfdata[0].stats.sampling_rate, self.wfdata[0].stats.delta)
        widget.axes.plot(time_ax, self.wfdata[0].data)

    def getID(self):
        try:
            return self.evtdata.get('resource_id').id
        except:
            return 'smi:bug/pylot/1234'


class GenericDataStructure(object):
    '''
    GenericDataBase type holds all information about the current data-
    base working on.
    '''
    def __init__(self, stexp=None, folderdepth=4, **kwargs):
        structExpression = []
        depth = 0
        while stexp is not os.path.sep:
            try:
                [stexp, tlexp] = os.path.split(stexp)
            except AttributeError:
                rootExpression = None
                structExpression = None
                break
            structExpression.append(tlexp)
            depth += 1
            if depth is folderdepth:
                rootExpression = stexp
                break
        structExpression.reverse()

        self.folderDepth = folderdepth
        self.dataBaseDict = {}


class SeiscompDataStructure(object):
    '''
    Dictionary containing the data access information for an SDS data archive:

    :param str dataType: Desired data type. Default: ``'waveform'``
    :param sdate, edate: Either date string or an instance of
         :class:`obspy.core.utcdatetime.UTCDateTime. Default: ``None``
    :type sdate, edate: str or UTCDateTime or None
    '''

    # Data type options
    __typeOptions = {'waveform': 'D',  # Waveform data
                          'detect': 'E',  # Detection data
                          'log': 'L',  # Log data
                          'timing': 'T',  # Timing data
                          'calib': 'C',  # Calibration data
                          'resp': 'R',  # Response data
                          'opaque': 'O'  # Opaque data
                          }

    def __init__(self, dataType='waveform', sdate=None, edate=None, **kwargs):
        # imports
        from obspy.core import UTCDateTime

        def checkDate(date):
            if not isinstance(date, UTCDateTime):
                return True
            return False

        try:
            if checkDate(sdate):
                sdate = UTCDateTime(sdate)
            if checkDate(edate):
                edate = UTCDateTime(edate)
        except TypeError:
            edate = UTCDateTime()
            halfyear = UTCDateTime('1970-07-01')
            sdate = UTCDateTime(edate - halfyear)
            del halfyear

        year = ''
        if not edate.year == sdate.year:
            nyears = edate.year - sdate.year
            for yr in range(nyears):
                year += '{0:04d},'.format(sdate.year+yr)
            year = '{'+year[:-1]+'}'
        else:
            year = '{0:04d}'.format(sdate.year)

        if dataType in self.__typeOptions.keys():
            self.dataType = dataType
        else:
            self.dataType = 'waveform'  # default value for dataType
            print '''Warning: Selected datatype ('%s') not available.\n
                     Using 'waveform' instead!'''.format(dataType)

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
                    if key in self.getSDSFields().keys():
                        self.getSDSFields()[key] = value
                    else:
                        raise KeyError('unknown SDS wildcard: %s.' % key)
                except KeyError, e:
                    errmsg = ''
                    errmsg += 'WARNING:\n'
                    errmsg += 'unable to set values for SDS fields\n'
                    errmsg += '%s; desired value was: %s\n' % (e, value)
                    print errmsg

    def getType(self):
        return self.__typeOptions[self.dataType]

    def getSDSFields(self):
        return self.__sdsFields

    def expandDataPath(self):
        fullChan = '{0}.{1}'.format(self.getSDSFields()['CHAN'], self.getType())
        dataPath = os.path.join(self.getSDSFields()['SDSdir'],
                                self.getSDSFields()['YEAR'],
                                self.getSDSFields()['NET'],
                                self.getSDSFields()['STA'],
                                fullChan,
                                '*{0}'.format(self.getSDSFields()['DAY'])
                                )
        return dataPath
