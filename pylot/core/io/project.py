import os

from pylot.core.util.event import Event


class Project(object):
    '''
    Pickable class containing information of a PyLoT project, like event lists and file locations.
    '''

    # TODO: remove rootpath
    def __init__(self):
        self.eventlist = []
        self.location = None
        self.rootpath = None
        self.datapath = None
        self.dirty = False
        self.parameter = None
        self._table = None

    def add_eventlist(self, eventlist):
        '''
        Add events from an eventlist containing paths to event directories.
        Will skip existing paths.
        '''
        if len(eventlist) == 0:
            return
        for item in eventlist:
            event = Event(item)
            event.rootpath = self.parameter['rootpath']
            event.database = self.parameter['database']
            event.datapath = self.parameter['datapath']
            if not event.path in self.getPaths():
                self.eventlist.append(event)
                self.setDirty()
            else:
                print('Skipping event with path {}. Already part of project.'.format(event.path))
        self.eventlist.sort(key=lambda x: x.pylot_id)
        self.search_eventfile_info()

    def remove_event(self, event):
        self.eventlist.remove(event)

    def remove_event_by_id(self, eventID):
        for event in self.eventlist:
            if eventID in str(event.resource_id):
                self.remove_event(event)
                break

    def read_eventfile_info(self, filename, separator=','):
        '''
        Try to read event information from file (:param:filename) comparing specific event datetimes.
        File structure (each row): event, date, time, magnitude, latitude, longitude, depth
        separated by :param:separator each.
        '''
        with open(filename, 'r') as infile:
            for line in infile.readlines():
                eventID, date, time, mag, lat, lon, depth = line.split(separator)[:7]
                # skip first line
                try:
                    day, month, year = date.split('/')
                except:
                    continue
                year = int(year)
                # hardcoded, if year only consists of 2 digits (e.g. 16 instead of 2016)
                if year < 100:
                    year += 2000
                datetime = '{}-{}-{}T{}'.format(year, month, day, time)
                try:
                    datetime = UTCDateTime(datetime)
                except Exception as e:
                    print(e, datetime, filename)
                    continue
                for event in self.eventlist:
                    if eventID in str(event.resource_id) or eventID in event.origins:
                        if event.origins:
                            origin = event.origins[0]  # should have only one origin
                            if origin.time == datetime:
                                origin.latitude = float(lat)
                                origin.longitude = float(lon)
                                origin.depth = float(depth)
                            else:
                                continue
                        elif not event.origins:
                            origin = Origin(resource_id=event.resource_id,
                                            time=datetime, latitude=float(lat),
                                            longitude=float(lon), depth=float(depth))
                            event.origins.append(origin)
                        event.magnitudes.append(Magnitude(resource_id=event.resource_id,
                                                          mag=float(mag),
                                                          mag_type='M'))
                        break

    def search_eventfile_info(self):
        '''
        Search all datapaths in rootpath for filenames with given file extension fext
        and try to read event info from it
        '''
        datapaths = []
        fext = '.csv'
        for event in self.eventlist:
            if not event.datapath in datapaths:
                datapaths.append(event.datapath)
        for datapath in datapaths:
            # datapath = os.path.join(self.rootpath, datapath)
            if os.path.isdir(datapath):
                for filename in os.listdir(datapath):
                    filename = os.path.join(datapath, filename)
                    if os.path.isfile(filename) and filename.endswith(fext):
                        try:
                            self.read_eventfile_info(filename)
                        except Exception as e:
                            print('Failed on reading eventfile info from file {}: {}'.format(filename, e))
            else:
                print("Directory %s does not exist!" % datapath)

    def getPaths(self):
        '''
        Returns paths (eventlist) of all events saved in the project.
        '''
        paths = []
        for event in self.eventlist:
            paths.append(event.path)
        return paths

    def setDirty(self, value=True):
        self.dirty = value

    def getEventFromPath(self, path):
        '''
        Search for an event in the project by event path.
        '''
        for event in self.eventlist:
            if event.path == path:
                return event

    def save(self, filename=None):
        '''
        Save PyLoT Project to a file.
        Can be loaded by using project.load(filename).
        '''
        try:
            import pickle
        except ImportError:
            import _pickle as pickle

        if filename:
            self.location = filename
        else:
            filename = self.location

        table = self._table  # MP: see below
        try:
            outfile = open(filename, 'wb')
            self._table = []  # MP: Workaround as long as table cannot be saved as part of project
            pickle.dump(self, outfile, protocol=pickle.HIGHEST_PROTOCOL)
            self.setDirty(False)
            self._table = table  # MP: see above
            return True
        except Exception as e:
            print('Could not pickle PyLoT project. Reason: {}'.format(e))
            self.setDirty()
            self._table = table  # MP: see above
            return False

    @staticmethod
    def load(filename):
        '''
        Load project from filename.
        '''
        import pickle
        infile = open(filename, 'rb')
        project = pickle.load(infile)
        infile.close()
        project.location = filename
        print('Loaded %s' % filename)
        return project
