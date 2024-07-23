#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from dataclasses import dataclass, field
from typing import List, Union

from obspy import UTCDateTime

from pylot.core.io.event import EventData
from pylot.core.io.waveformdata import WaveformData

@dataclass
class Data:
    event_data: EventData = field(default_factory=EventData)
    waveform_data: WaveformData = field(default_factory=WaveformData)
    _parent: Union[None, 'QtWidgets.QWidget'] = None

    def __init__(self, parent=None, evtdata=None):
        self._parent = parent
        self.event_data = EventData(evtdata)
        self.waveform_data = WaveformData()

    def __str__(self):
        return str(self.waveform_data.wfdata)

    def __add__(self, other):
        if not isinstance(other, Data):
            raise TypeError("Operands must be of type 'Data'")
        if self.event_data.is_new() and other.event_data.is_new():
            pass
        elif other.event_data.is_new():
            new_picks = other.event_data.evtdata.picks
            old_picks = self.event_data.evtdata.picks
            old_picks.extend([pick for pick in new_picks if pick not in old_picks])
        elif self.event_data.is_new():
            return other + self
        elif self.event_data.get_id() == other.event_data.get_id():
            other.event_data.set_new()
            return self + other
        else:
            raise ValueError("Both Data objects have differing unique Event identifiers")
        return self

    def get_parent(self):
        return self._parent

    def filter_wf_data(self, **kwargs):
        self.waveform_data.wfdata.detrend('linear')
        self.waveform_data.wfdata.taper(0.02, type='cosine')
        self.waveform_data.wfdata.filter(**kwargs)
        self.waveform_data.dirty = True

    def set_wf_data(self, fnames: List[str], fnames_alt: List[str] = None, check_rotated=False, metadata=None, tstart=0, tstop=0):
        return self.waveform_data.set_wf_data(fnames, fnames_alt, check_rotated, metadata, tstart, tstop)

    def reset_wf_data(self):
        self.waveform_data.reset_wf_data()

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
