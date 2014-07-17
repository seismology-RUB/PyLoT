#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from obspy.core import (read, Stream)
from obspy.core.event import Event


class Data(object):
    '''
    Data container class providing ObSpy-Stream and -Event instances as
    variables.
    
    :type parent: PySide.QtGui.QWidget object, optional
    :param parent: A PySide.QtGui.QWidget class utilized when
    called by a GUI to display a PySide.QtGui.QMessageBox instead of printing
    to standard out.
    '''

    def __init__(self, parent=None, wfdata=None, evtdata=None):
        if wfdata is not None and isinstance(wfdata, Stream):
            self.wfdata = wfdata
        elif wfdata is not None:
            try:
                self.wfdata = read(wfdata)
            except IOError, e:
                msg = 'An I/O error occured while loading data!'
                inform = 'Variable wfdata will be empty.'
                details = '{0}'.format(e)
                if parent is not None:
                    from PySide.QtGui import QMessageBox
                    warnio = QMessageBox(parent=parent)
                    warnio.setText(msg)
                    warnio.setDetailedText(details)
                    warnio.setStandarButtons(QMessageBox.Ok)
                    warnio.setIcon(QMessageBox.Warning)
                else:
                    print msg, '\n', details
                self.wfdata = Stream()
        else:
            self.wfdata = Stream()
        if evtdata is not None and isinstance(evtdata, Event):
            self.evtdata = evtdata
        else:
            self.evtdata = Event()


class GenericDataStructure(object):
    '''
    GenericDataBase type holds all information about the current data-
    base working on.
    '''
    def __init__(self, stexp=None, folderdepth=4, **kwargs):
        dbRegExp = {}
        structExpression = []
        depth = 0
        while stexp is not os.path.sep:
            [stexp, tlexp] = os.path.split(stexp)
            structExpression.append(tlexp)
            depth += 1
            if depth is folderdepth:
                rootExpression = stexp
                break
        structExpression.reverse()

        self.folderDepth = folderdepth
        self.dataBaseDict = kwargs


class SeiscompDataStructure(object):
    '''
    Dictionary containing the data acces information for an SDS data archive:

    :param str dataType: Desired data type. Default: ``'waveform'``
    :param sdate, edate: Either date string or an instance of
         :class:`obspy.core.utcdatetime.UTCDateTime. Default: ``None``
    :type sdate, edate: str or UTCDateTime or None
    '''

    def __init__(self, dataType='waveform', sdate=None, edate=None, **kwargs):
        # imports
        from obspy.core import UTCDateTime
        import numpy as np

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
            sdate = edate - np.pi*1e7/2

        year = ''
        if not edate.year == sdate.year:
            nyears = edate.year - sdate.year
            for yr in range(nyears):
                year += '{0:04d},'.format(sdate.year+yr)
            year = '{'+year[:-1]+'}'
        else:
            year = '{0:04d},'.format(sdate.year)

        # Data type options
        self.__typeOptions = {'waveform': 'D',  # Waveform data
                              'detect': 'E',  # Detection data
                              'log': 'L',  # Log data
                              'timing': 'T',  # Timing data
                              'calib': 'C',  # Calibration data
                              'resp': 'R',  # Response data
                              'opaque': 'O'  # Opaque data
                              }

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
                            'TYPE': self._getType(),  # 1 character
                            'LOC': '',  # up to 8 characters
                            'DAY': '{0:03d}'.format(sdate.julday)  # 3 digits
                            }
        self.modifiyFields(**kwargs)

    def modifiyFields(self, **kwargs):
        if kwargs and isinstance(kwargs, dict):
            for key, value in kwargs.iteritems():
                key = str(key)
                value = str(value)
                try:
                    if key in self.__sdsFields.keys():
                        self.__sdsFields[key] = str(value)
                    else:
                        raise KeyError('unknown SDS wildcard: %s.' % key)
                except KeyError, e:
                    errmsg = ''
                    errmsg += 'WARNING:\n'
                    errmsg += 'unable to set values for SDS fields\n'
                    errmsg += '%s; desired value was: %s\n' % (e, value)
                    print errmsg

    def _getType(self):
        return self.__typeOptions[self.dataType]

    def expandDataPath(self):
        fullChan = '{0}.{1}'.format(self.__sdsFields['CHAN'], self._getType())
        dataPath = os.path.join(self.__sdsFields['SDSdir'],
                                self.__sdsFields['YEAR'],
                                self.__sdsFields['NET'],
                                self.__sdsFields['STA'],
                                fullChan,
                                '*{0}'.format(self.__sdsFields['DAY'])
                                )
        return dataPath
