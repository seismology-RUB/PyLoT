from obspy import UTCDateTime
from obspy.core import event as ope

from pylot.core.util.utils import getLogin, getHash


def create_amplitude(pickID, amp, unit, category, cinfo):
    '''

    :param pickID:
    :param amp:
    :param unit:
    :param category:
    :param cinfo:
    :return:
    '''
    amplitude = ope.Amplitude()
    amplitude.creation_info = cinfo
    amplitude.generic_amplitude = amp
    amplitude.unit = ope.AmplitudeUnit(unit)
    amplitude.type = ope.AmplitudeCategory(category)
    amplitude.pick_id = pickID
    return amplitude


def create_arrival(pickresID, cinfo, phase, azimuth=None, dist=None):
    '''
    create_arrival - function to create an Obspy Arrival

    :param pickresID: Resource identifier of the created pick
    :type pickresID: :class: `~obspy.core.event.ResourceIdentifier` object
    :param cinfo: An ObsPy :class: `~obspy.core.event.CreationInfo` object
    holding information on the creation of the returned object
    :type cinfo: :class: `~obspy.core.event.CreationInfo` object
    :param phase: name of the arrivals seismic phase
    :type phase: str
    :param azimuth: azimuth between source and receiver
    :type azimuth: float or int, optional
    :param dist: distance between source and receiver
    :type dist: float or int, optional
    :return: An ObsPy :class: `~obspy.core.event.Arrival` object
    '''
    arrival = ope.Arrival()
    arrival.creation_info = cinfo
    arrival.pick_id = pickresID
    arrival.phase = phase
    if azimuth is not None:
        arrival.azimuth = float(azimuth) if azimuth > -180 else azimuth + 360.
    else:
        arrival.azimuth = azimuth
    arrival.distance = dist
    return arrival


def create_creation_info(agency_id=None, creation_time=None, author=None):
    '''
    get creation info of obspy event
    :param agency_id:
    :param creation_time:
    :param author:
    :return:
    '''
    if author is None:
        author = getLogin()
    if creation_time is None:
        creation_time = UTCDateTime()
    return ope.CreationInfo(agency_id=agency_id, author=author,
                            creation_time=creation_time)


def create_event(origintime, cinfo, originloc=None, etype='earthquake',
                 resID=None, authority_id=None):
    '''
    create_event - funtion to create an ObsPy Event

    :param origintime: the events origintime
    :type origintime: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :param cinfo: An ObsPy :class: `~obspy.core.event.CreationInfo` object
        holding information on the creation of the returned object
    :type cinfo: :class: `~obspy.core.event.CreationInfo` object
    :param originloc: tuple containing the location of the origin
        (LAT, LON, DEP) affiliated with the event which is created
    :type originloc: tuple, list
    :param etype: Event type str object. converted via ObsPy to a valid event
        type string.
    :type etype: str
    :param resID: Resource identifier of the created event
    :type resID: :class: `~obspy.core.event.ResourceIdentifier` object, str
    :param authority_id: name of the institution carrying out the processing
    :type authority_id: str
    :return: An ObsPy :class: `~obspy.core.event.Event` object
    '''

    if originloc is not None:
        o = create_origin(origintime, cinfo,
                          originloc[0], originloc[1], originloc[2])
    else:
        o = None
    if not resID:
        resID = create_resourceID(origintime, etype, authority_id)
    elif isinstance(resID, str):
        resID = create_resourceID(origintime, etype, authority_id, resID)
    elif not isinstance(resID, ope.ResourceIdentifier):
        raise TypeError("unsupported type(resID) for resource identifier "
                        "generation: %s" % type(resID))
    event = ope.Event(resource_id=resID)
    event.creation_info = cinfo
    event.event_type = etype
    if o:
        event.origins = [o]
    return event


def create_magnitude(originID, cinfo):
    '''
    create_magnitude - function to create an ObsPy Magnitude object
    :param originID:
    :type originID:
    :param cinfo:
    :type cinfo:
    :return:
    '''
    magnitude = ope.Magnitude()
    magnitude.creation_info = cinfo
    magnitude.origin_id = originID
    return magnitude


def create_origin(origintime, cinfo, latitude, longitude, depth):
    '''
    create_origin - function to create an ObsPy Origin
    :param origintime: the origins time of occurence
    :type origintime: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :param cinfo:
    :type cinfo:
    :param latitude: latitude in decimal degree of the origins location
    :type latitude: float
    :param longitude: longitude in decimal degree of the origins location
    :type longitude: float
    :param depth: hypocentral depth of the origin
    :type depth: float
    :return: An ObsPy :class: `~obspy.core.event.Origin` object
    '''

    assert isinstance(origintime, UTCDateTime), "origintime has to be " \
                                                "a UTCDateTime object, but " \
                                                "actually is of type " \
                                                "'%s'" % type(origintime)

    origin = ope.Origin()
    origin.time = origintime
    origin.creation_info = cinfo
    origin.latitude = latitude
    origin.longitude = longitude
    origin.depth = depth
    return origin


def create_pick(origintime, picknum, picktime, eventnum, cinfo, phase, station,
                wfseedstr, authority_id):
    '''
    create_pick - function to create an ObsPy Pick

    :param origintime:
    :type origintime:
    :param picknum: number of the created pick
    :type picknum: int
    :param picktime:
    :type picktime:
    :param eventnum: human-readable event identifier
    :type eventnum: str
    :param cinfo: An ObsPy :class: `~obspy.core.event.CreationInfo` object
        holding information on the creation of the returned object
    :type cinfo: :class: `~obspy.core.event.CreationInfo` object
    :param phase: name of the arrivals seismic phase
    :type phase: str
    :param station: name of the station at which the seismic phase has been
        picked
    :type station: str
    :param wfseedstr: A SEED formatted string of the form
        network.station.location.channel in order to set a referenced waveform
    :type wfseedstr: str, SEED formatted
    :param authority_id: name of the institution carrying out the processing
    :type authority_id: str
    :return: An ObsPy :class: `~obspy.core.event.Pick` object
    '''
    pickID = eventnum + '_' + station.strip() + '/{0:03d}'.format(picknum)
    pickresID = create_resourceID(origintime, 'pick', authority_id, pickID)
    pick = ope.Pick()
    pick.resource_id = pickresID
    pick.time = picktime
    pick.creation_info = cinfo
    pick.phase_hint = phase
    pick.waveform_id = ope.ResourceIdentifier(id=wfseedstr, prefix='file:/')
    return pick


def create_resourceID(timetohash, restype, authority_id=None, hrstr=None):
    '''
    create unique resource id
    :param timetohash: event origin time to hash
    :type timetohash: class: `~obspy.core.utcdatetime.UTCDateTime` object
    :param restype: type of the resource, e.g. 'orig', 'earthquake' ...
    :type restype: str
    :param authority_id: name of the institution carrying out the processing
    :type authority_id: str, optional
    :param hrstr:
    :type hrstr:
    :return:
    '''
    assert isinstance(timetohash, UTCDateTime), "'timetohash' is not an ObsPy" \
                                                "UTCDateTime object"
    hid = getHash(timetohash)
    if hrstr is None:
        resID = ope.ResourceIdentifier(restype + '/' + hid[0:6])
    else:
        resID = ope.ResourceIdentifier(restype + '/' + hrstr)
    if authority_id is not None:
        resID.convertIDToQuakeMLURI(authority_id=authority_id)
    return resID
