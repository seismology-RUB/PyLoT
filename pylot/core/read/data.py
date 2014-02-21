#!/usr/bin/env python
#
# Provide user the opportunity to read arbitrary organized database
# types. This means e.g. seiscomp data structure (SDS) or event based
# EGELADOS structure.
#

import os


class GenericDataBase(object):
    '''
    GenericDataBase type holds all information about the current data-
    base working on.
    '''
    def __init__(self, stexp=None, **kwargs):
        dbRegExp = {}

        structExpression = os.path.split(stexp)
        self.dataBaseDict = kwargs


class SeiscompDataStructure(object):

    def __init__(self, dataType='waveform', **kwargs):
        '''
        Object initialization method:

        :type dataType: str
        :param dataType: Desired data type. Default: ``'waveform'``
        '''
        # imports
        from obspy.core import UTCDateTime

        now = UTCDateTime()
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

        # SDS fields' default values
        # definitions from
        # http://www.seiscomp3.org/wiki/doc/applications/slarchive/SDS
        self.__sdsFields = {'SDSdir': '/data/SDS',  # base directory
                            'YEAR': '{0:04d}'.format(now.year),  # 4 digits
                            'NET': '??',  # up to 8 characters
                            'STA': '????',  # up to 8 characters
                            'CHAN': 'HH?',  # up to 8 characters
                            'TYPE': self._getType(),  # 1 character
                            'LOC': '',  # up to 8 characters
                            'DAY': '{0:03d}'.format(now.julday)  # 3 digit doy
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
