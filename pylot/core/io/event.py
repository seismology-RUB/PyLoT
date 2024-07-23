import copy
from dataclasses import dataclass
from typing import Union

from obspy import read_events
from obspy.core.event import Event as ObsPyEvent


@dataclass
class EventData:
    evtdata: Union[ObsPyEvent, None] = None
    _new: bool = False

    def __init__(self, evtdata=None):
        self.set_event_data(evtdata)

    def set_event_data(self, evtdata):
        if isinstance(evtdata, ObsPyEvent):
            self.evtdata = evtdata
        elif isinstance(evtdata, dict):
            self.evtdata = self.read_pilot_event(**evtdata)
        elif isinstance(evtdata, str):
            self.evtdata = self.read_event_file(evtdata)
        else:
            self.set_new()
            self.evtdata = ObsPyEvent(picks=[])

    def read_event_file(self, evtdata: str) -> ObsPyEvent:
        try:
            cat = read_events(evtdata)
            if len(cat) != 1:
                raise ValueError(f'ambiguous event information for file: {evtdata}')
            return cat[0]
        except TypeError as e:
            self.handle_event_file_error(e, evtdata)

    def handle_event_file_error(self, e: TypeError, evtdata: str):
        if 'Unknown format for file' in str(e):
            if 'PHASES' in evtdata:
                picks = self.picksdict_from_pilot(evtdata)
                evtdata = ObsPyEvent(picks=self.picks_from_picksdict(picks))
            elif 'LOC' in evtdata:
                raise NotImplementedError('PILOT location information read support not yet implemented.')
            elif 'event.pkl' in evtdata:
                evtdata = self.qml_from_obspy_dmt(evtdata)
            else:
                raise e
        else:
            raise e

    def set_new(self):
        self._new = True

    def is_new(self) -> bool:
        return self._new

    def get_picks_str(self) -> str:
        return '\n'.join(str(pick) for pick in self.evtdata.picks)

    def replace_origin(self, event: ObsPyEvent, force_overwrite: bool = False):
        if self.evtdata.origins or force_overwrite:
            event.origins = self.evtdata.origins

    def replace_magnitude(self, event: ObsPyEvent, force_overwrite: bool = False):
        if self.evtdata.magnitudes or force_overwrite:
            event.magnitudes = self.evtdata.magnitudes

    def replace_picks(self, event: ObsPyEvent, picktype: str):
        checkflag = 1
        picks = event.picks
        for j, pick in reversed(list(enumerate(picks))):
            if picktype in str(pick.method_id.id):
                picks.pop(j)
                checkflag = 2

        if checkflag > 0:
            for pick in self.evtdata.picks:
                if picktype in str(pick.method_id.id):
                    picks.append(pick)

    def get_id(self) -> Union[str, None]:
        try:
            return self.evtdata.resource_id.id
        except:
            return None

    def apply_event_data(self, data, typ='pick'):
        if typ == 'pick':
            self.apply_picks(data)
        elif typ == 'event':
            self.apply_event(data)

    def apply_picks(self, picks):
        self.evtdata.picks = picks

    def apply_event(self, event: ObsPyEvent):
        if self.is_new():
            self.evtdata = event
        else:
            old_event = self.evtdata
            if old_event.resource_id == event.resource_id:
                picks = copy.deepcopy(old_event.picks)
                event = self.merge_picks(event, picks)
                old_event.update(event)
            else:
                print(f"WARNING: Mismatch in event resource id's: {old_event.resource_id} and {event.resource_id}")
