#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from obspy import UTCDateTime
from obspy.core.event import Event as ObsPyEvent
from obspy.core.event import Origin, ResourceIdentifier
from pylot.core.io.phases import picks_from_picksdict
from pylot.core.util.obspyDMT_interface import qml_from_obspyDMT


class Event(ObsPyEvent):
    '''
    Pickable class derived from ~obspy.core.event.Event containing information on a single event.
    '''

    def __init__(self, path):
        """
        Initialize event by event directory
        :param path: path to event directory
        :type path: str
        """
        self.pylot_id = path.split('/')[-1]
        # initialize super class
        super(Event, self).__init__(resource_id=ResourceIdentifier('smi:local/' + self.pylot_id))
        self.path = path
        self.database = path.split('/')[-2]
        self.datapath = path.split('/')[-3]
        self.rootpath = '/' + os.path.join(*path.split('/')[:-3])
        self.pylot_autopicks = {}
        self.pylot_picks = {}
        self.notes = ''
        self._testEvent = False
        self._refEvent = False
        self.get_notes()
        self.get_obspy_event_info()

    def get_notes_path(self):
        """
        Notes files is freely editable by the user and can contain notes regarding the event
        :return: path to notes file
        :rtype: str
        """
        notesfile = os.path.join(self.path, 'notes.txt')
        return notesfile

    def get_obspy_event_info(self):
        infile_pickle = os.path.join(self.path, 'info/event.pkl')
        if not os.path.isfile(infile_pickle):
            return
        try:
            event_dmt = qml_from_obspyDMT(infile_pickle)
        except Exception as e:
            print('Could not get obspy event info: {}'.format(e))
            return
        self.magnitudes = event_dmt.magnitudes
        self.origins = event_dmt.origins

    def get_notes(self):
        """
        set self.note attribute to content of notes file
        :return:
        :rtype: None
        """
        notesfile = self.get_notes_path()
        if os.path.isfile(notesfile):
            with open(notesfile) as infile:
                path = str(infile.readlines()[0].split('\n')[0])
                text = '[eventInfo: ' + path + ']'
                self.addNotes(text)
                try:
                    datetime = UTCDateTime(path.split('/')[-1])
                    origin = Origin(resource_id=self.resource_id, time=datetime, latitude=0, longitude=0, depth=0)
                    self.origins.append(origin)
                except:
                    pass

    def addNotes(self, notes):
        """
        Set new notes string
        :param notes: notes to save in Event object
        :type notes: str
        :return:
        :rtype: None
        """
        self.notes = str(notes)

    def clearNotes(self):
        """
        Clear event notes
        :return:
        :rtype: None
        """
        self.notes = None

    def isRefEvent(self):
        """
        Return reference event flag
        :return: True if event is refence event
        :rtype: bool
        """
        return self._refEvent

    def isTestEvent(self):
        """
        Return test event flag
        :return: True if event is test event
        :rtype: bool
        """
        return self._testEvent

    def setRefEvent(self, bool):
        """
        Set reference event flag
        :param bool: new reference event flag
        :type bool: bool
        :return:
        :rtype: None
        """
        self._refEvent = bool
        if bool: self._testEvent = False

    def setTestEvent(self, bool):
        """
        Set test event flag
        :param bool: new test event flag
        :type bool: bool
        :return:
        :rtype: None
        """
        self._testEvent = bool
        if bool: self._refEvent = False

    def clearObsPyPicks(self, picktype):
        """
        Remove picks of a certain type from event
        :param picktype: type of picks to remove, 'auto' or 'manual'
        :type picktype: str
        :return:
        :rtype: None
        """
        for index, pick in reversed(list(enumerate(self.picks))):
            if picktype in str(pick.method_id):
                self.picks.pop(index)

    def addPicks(self, picks):
        """
        add pylot picks and overwrite existing ones
        :param picks: picks to add to event in pick dictionary
        :type picks: dict
        :return:
        :rtype: None
        """
        for station in picks:
            self.pylot_picks[station] = picks[station]
        # add ObsPy picks (clear old manual and copy all new manual from pylot)
        self.clearObsPyPicks('manual')
        self.picks += picks_from_picksdict(self.pylot_picks)

    def addAutopicks(self, autopicks):
        """
        Add automatic picks to event
        :param autopicks: automatic picks to add to event
        :return:
        :rtype: None
        """
        for station in autopicks:
            self.pylot_autopicks[station] = autopicks[station]
        # add ObsPy picks (clear old auto and copy all new auto from pylot)
        self.clearObsPyPicks('auto')
        self.picks += picks_from_picksdict(self.pylot_autopicks)

    def setPick(self, station, pick):
        """
        Set pick for a station
        :param station: station name
        :type station: str
        :param pick:
        :type pick: dict
        :return:
        :rtype:
        """
        if pick:
            self.pylot_picks[station] = pick
        else:
            try:
                if station in self.pylot_picks:
                    self.pylot_picks.pop(station)
            except Exception as e:
                print('Could not remove pick {} from station {}: {}'.format(pick, station, e))
        self.clearObsPyPicks('manual')
        self.picks += picks_from_picksdict(self.pylot_picks)

    def setPicks(self, picks):
        """
        Set pylot picks and delete and overwrite all existing
        :param picks: new picks
        :type picks: dict
        :return:
        :rtype: None
        """
        self.pylot_picks = picks
        self.clearObsPyPicks('manual')
        self.picks += picks_from_picksdict(self.pylot_picks)

    def getPick(self, station):
        """
        Get pick at station
        :param station: station name
        :type station: str
        :return: pick dictionary of station
        :rtype: dict
        """
        if station in self.pylot_picks.keys():
            return self.pylot_picks[station]

    def getPicks(self):
        """
        Return pylot picks
        :return:
        :rtype: dict
        """
        return self.pylot_picks

    def setAutopick(self, station, pick):
        """
        Set pick at station
        :param station: station name
        :type station: str
        :param pick:
        :type pick: dict
        :return:
        :rtype: None
        """
        if pick:
            self.pylot_autopicks[station] = pick
        else:
            try:
                if station in self.pylot_autopicks:
                    self.pylot_autopicks.pop(station)
            except Exception as e:
                print('Could not remove pick {} from station {}: {}'.format(pick, station, e))
        self.clearObsPyPicks('auto')
        self.picks += picks_from_picksdict(self.pylot_autopicks)

    def setAutopicks(self, picks):
        """
        Set pylot picks and delete and overwrite all existing
        :param picks:  new picks
        :type picks: dict
        :return:
        :rtype: None
        """
        self.pylot_autopicks = picks
        self.clearObsPyPicks('auto')
        self.picks += picks_from_picksdict(self.pylot_autopicks)

    def getAutopick(self, station):
        """
        Return autopick at station
        :param station: station name
        :type station: str
        :return: pick dictionary
        :rtype: dict
        """
        if station in self.pylot_autopicks.keys():
            return self.pylot_autopicks[station]

    def getAutopicks(self):
        """
        Get autopicks of event
        :return: dict containing automatic picks
        :rtype: dict
        """
        return self.pylot_autopicks

    def save(self, filename):
        """
        Save PyLoT Event to a file.
        Can be loaded by using event.load(filename).
        Uses pickling to save event object to file
        :param filename: filename to save project under
        :type filename: str
        :return:
        :rtype: None
        """
        try:
            import cPickle
        except ImportError:
            import _pickle as cPickle

        try:
            outfile = open(filename, 'wb')
            cPickle.dump(self, outfile, -1)
        except Exception as e:
            print('Could not pickle PyLoT event. Reason: {}'.format(e))

    @staticmethod
    def load(filename):
        """
        Load project from filename
        :param filename: to load event file
        :type filename: str
        :return: event loaded from file
        :rtype: Event
        """
        try:
            import cPickle
        except ImportError:
            import _pickle as cPickle
        infile = open(filename, 'rb')
        event = cPickle.load(infile)
        print('Loaded %s' % filename)
        return event
