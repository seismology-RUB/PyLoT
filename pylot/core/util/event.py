#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from obspy import UTCDateTime
from obspy.core.event import Event as ObsPyEvent
from obspy.core.event import Origin, ResourceIdentifier
from pylot.core.io.phases import picks_from_picksdict


class Event(ObsPyEvent):
    '''
    Pickable class derived from ~obspy.core.event.Event containing information on a single event.
    '''

    def __init__(self, path):
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

    def get_notes_path(self):
        notesfile = os.path.join(self.path, 'notes.txt')
        return notesfile

    def get_notes(self):
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
        self.notes = str(notes)

    def clearNotes(self):
        self.notes = None

    def isRefEvent(self):
        return self._refEvent

    def isTestEvent(self):
        return self._testEvent

    def setRefEvent(self, bool):
        self._refEvent = bool
        if bool: self._testEvent = False

    def setTestEvent(self, bool):
        self._testEvent = bool
        if bool: self._refEvent = False

    def clearObsPyPicks(self, picktype):
        for index, pick in reversed(list(enumerate(self.picks))):
            if picktype in str(pick.method_id):
                self.picks.pop(index)

    def addPicks(self, picks):
        '''
        add pylot picks and overwrite existing ones
        '''
        for station in picks:
            self.pylot_picks[station] = picks[station]
        # add ObsPy picks (clear old manual and copy all new manual from pylot)
        self.clearObsPyPicks('manual')
        self.picks += picks_from_picksdict(self.pylot_picks)

    def addAutopicks(self, autopicks):
        for station in autopicks:
            self.pylot_autopicks[station] = autopicks[station]
        # add ObsPy picks (clear old auto and copy all new auto from pylot)
        self.clearObsPyPicks('auto')
        self.picks += picks_from_picksdict(self.pylot_autopicks)

    def setPick(self, station, pick):
        if pick:
            self.pylot_picks[station] = pick
        else:
            try:
                self.pylot_picks.pop(station)
            except Exception as e:
                print('Could not remove pick {} from station {}: {}'.format(pick, station, e))
        self.clearObsPyPicks('manual')
        self.picks += picks_from_picksdict(self.pylot_picks)

    def setPicks(self, picks):
        '''
        set pylot picks and delete and overwrite all existing
        '''
        self.pylot_picks = picks
        self.clearObsPyPicks('manual')
        self.picks += picks_from_picksdict(self.pylot_picks)

    def getPick(self, station):
        if station in self.pylot_picks.keys():
            return self.pylot_picks[station]

    def getPicks(self):
        return self.pylot_picks

    def setAutopick(self, station, pick):
        if pick:
            self.pylot_autopicks[station] = pick
        else:
            try:
                self.pylot_autopicks.pop(station)
            except Exception as e:
                print('Could not remove pick {} from station {}: {}'.format(pick, station, e))
        self.clearObsPyPicks('auto')
        self.picks += picks_from_picksdict(self.pylot_autopicks)

    def setAutopicks(self, picks):
        '''
        set pylot picks and delete and overwrite all existing
        '''
        self.pylot_autopicks = picks
        self.clearObsPyPicks('auto')
        self.picks += picks_from_picksdict(self.pylot_autopicks)

    def getAutopick(self, station):
        if station in self.pylot_autopicks.keys():
            return self.pylot_autopicks[station]

    def getAutopicks(self):
        return self.pylot_autopicks

    def save(self, filename):
        '''
        Save PyLoT Event to a file. 
        Can be loaded by using event.load(filename).
        '''
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
        '''
        Load project from filename.
        '''
        try:
            import cPickle
        except ImportError:
            import _pickle as cPickle
        infile = open(filename, 'rb')
        event = cPickle.load(infile)
        print('Loaded %s' % filename)
        return event
