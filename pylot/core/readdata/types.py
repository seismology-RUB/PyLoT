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


class SeiscompDataStructure(GenericDataBase):

    def __init__(self, **kwargs):
        # Data type options
        self.__typeOptions = {'waveform': 'D',  # Waveform data
                              'detect': 'E',  # Detection data
                              'log': 'L',  # Log data
                              'timing': 'T',  # Timing data
                              'calib': 'C',  # Calibration data
                              'resp': 'R',  # Response data
                              'opaque': 'O'  # Opaque data
                              }
        # SDS fields' default values
        # definitions from
        # http://www.seiscomp3.org/wiki/doc/applications/slarchive/SDS
        self.__sdsFields = {'SDSdir': ['data', 'SDS'],  # base directory list
                            'YEAR': '1970',  # 4 digits
                            'NET': 'XX',  # up to 8 characters
                            'STA': 'XXXX',  # up to 8 characters
                            'CHAN': 'HHZ',  # up to 8 characters
                            'TYPE': self.getType('waveform'),  # 1 character
                            'LOC': '',  # up to 8 characters
                            'DAY': '001'  # 3 digit day of year
                            }
        pass

    def getType(self, typeName=None):
        pass
